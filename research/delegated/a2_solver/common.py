"""Shared cost ledger, calibration and tuning/final policy for A2 E4."""
from __future__ import annotations

import json
import math
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
SOLVER_DIR = Path(__file__).resolve().parent

C0 = 299792458.0
EPS0 = 8.8541878128e-12
SQRT2 = float(np.sqrt(2.0))

# Physical task constants (fixed by A2 protocol / TASK.md).
ALPHA_BOX = (0.03, 2.5)
X_BOX = (-0.7, 0.7)
POSE_LEVER = 1.5
GRID_H_N16 = 1.0 / 16.0

# rank ladder for prasc (capped at model n at runtime)
PRASC_RANK_LADDER = (8, 16, 24, 36, 64, 96, 128, 192, 256)
FIXED_RANKS = (16, 36, 64)
STAGE_FRACTIONS = (0.12, 0.28, 0.50, 1.0)
FREQ_STAGES = ((0,), (0, 1), (0, 1, 2), (0, 1, 2, 3))

# safety: final test seeds must never be touched from tuning
FINAL_SEED_RANGE = range(1001, 1021)
TUNING_SEED_RANGE = range(1, 11)


def path_for_seed(seed: int, allow_final: bool = False) -> None:
    if seed in FINAL_SEED_RANGE and not allow_final:
        raise RuntimeError(
            "final seed access attempted; final mode requires an explicit "
            "--frozen-config and an explicit final launch"
        )


def default_settings() -> dict[str, Any]:
    return {
        "optimizer": "L-BFGS-B",
        "maxls": 10,
        "ftol": 1e-8,
        "gtol": 1e-6,
        "prior_strength": 0.0,
        "prior_center": [0.5] * 9,
        "fixed_rank": None,
        "prasc_ranks": list(PRASC_RANK_LADDER),
        "tsom_append": False,
        "tsom_target_fraction": 0.25,
        "exact_fallback": True,
        "cert_residual_factor": 0.05,
        "cert_gradient_rel": 0.1,
        "init_max_units": 32,
        "control_rounds": True,
    }


def settings_hash(settings: dict[str, Any]) -> str:
    import hashlib

    blob = json.dumps(
        settings, sort_keys=True, separators=(",", ":"), default=_jsonable
    )
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()[:16]


def _jsonable(v: Any) -> Any:
    if isinstance(v, np.ndarray):
        return v.tolist()
    if isinstance(v, (np.integer, np.floating)):
        return float(v)
    if isinstance(v, tuple):
        return list(v)
    raise TypeError(f"not jsonable: {type(v)}")


def jsonable(v: Any) -> Any:
    return json.loads(json.dumps(v, default=_jsonable))


def pose_scale(lambda_min: float, lever: float = POSE_LEVER) -> np.ndarray:
    """Scale x=(tx,ty,theta) so the squared scaled norm equals
    (lever-arm metric)/lambda_min^2."""
    return np.array([1.0 / lambda_min, 1.0 / lambda_min, lever / lambda_min])


def lever_metric_error(x1: np.ndarray, x2: np.ndarray, lever: float = POSE_LEVER) -> float:
    d = np.asarray(x1, dtype=float) - np.asarray(x2, dtype=float)
    return float(np.sqrt(d[0] ** 2 + d[1] ** 2 + lever**2 * d[2] ** 2))


def frequency_row_slices() -> list[slice]:
    """Stable frequency-major row layout of physics forward output."""
    rows_per_freq = 3 * 2 * 12
    return [slice(i * rows_per_freq, (i + 1) * rows_per_freq) for i in range(4)]


def select_rows(y: np.ndarray, freq_ids: tuple[int, ...]) -> np.ndarray:
    slices = frequency_row_slices()
    parts = [np.asarray(y)[slices[fi]] for fi in sorted(freq_ids)]
    return np.concatenate(parts)


@dataclass
class CostLedger:
    """Equivalent-work accounting for one run.

    Unit policy (documented for the parent): one unit is one full-wave
    RHS solve column through an LU factorisation. Reduced RHS solves are
    converted using conservative time ratios from the tuning-only
    microbenchmark. LU factorisations, dense operator products, SVD/RRQR
    basis work and wall time are counted separately and are published raw;
    they are not silently converted into solve units.
    """

    full_rhs_solves: int = 0
    rhs_by_kind: dict[str, int] = field(default_factory=dict)
    full_factorizations: int = 0
    full_operator_products: int = 0
    basis_products_nr: int = 0  # dense n x r U-basis products in reduced path
    reduced_solves: int = 0
    reduced_factorizations: int = 0
    reduced_rank_sum: int = 0
    svd_count: int = 0
    rrqr_count: int = 0
    exact_acceptance_checks: int = 0
    wall_seconds: float = 0.0
    _units: float = 0.0
    calibration: Any = None

    def charge_rhs(self, n: int, kind: str = "other") -> None:
        self.full_rhs_solves += int(n)
        self.rhs_by_kind[kind] = self.rhs_by_kind.get(kind, 0) + int(n)
        self._units += float(n)

    def charge_model_work(self, work: dict[str, Any]) -> None:
        rhs = int(work.get("rhs_solves_total", 0))
        self.full_rhs_solves += rhs
        for k, v in (work.get("rhs_solves") or {}).items():
            self.rhs_by_kind[k] = self.rhs_by_kind.get(k, 0) + int(v)
        self._units += rhs
        self.full_factorizations += int(work.get("factorizations", 0))
        self.full_operator_products += int(work.get("operator_products", 0))
        self.wall_seconds += float(work.get("wall_seconds", 0.0))

    def charge_factorizations(self, n: int) -> None:
        self.full_factorizations += int(n)

    def charge_operator_products(self, n: int) -> None:
        self.full_operator_products += int(n)

    def charge_basis_products(self, n: int) -> None:
        self.basis_products_nr += int(n)

    def charge_reduced(
        self, solves: int = 0, factorizations: int = 0, rank: int | None = None
    ) -> None:
        """Charge reduced linear algebra converted by calibration."""
        self.reduced_solves += int(solves)
        self.reduced_factorizations += int(factorizations)
        if rank is not None:
            self.reduced_rank_sum += int(solves) * int(rank)
        cal = self.calibration
        if cal is None:
            # no calibration yet: treat a reduced solve as cheap, not free.
            self._units += 0.0
            return
        r = int(rank or cal.nominal_n)
        self._units += float(solves) * cal.reduced_solve_ratio(r)
        self._units += float(factorizations) * cal.reduced_qr_ratio(r)

    def charge_svd(self, n: int = 1) -> None:
        self.svd_count += int(n)

    def charge_rrqr(self, n: int = 1) -> None:
        self.rrqr_count += int(n)

    def charge_exact_acceptance(self, n: int = 1) -> None:
        self.exact_acceptance_checks += int(n)

    def add_wall(self, seconds: float) -> None:
        self.wall_seconds += float(seconds)

    @property
    def units(self) -> float:
        return float(self._units)

    def snapshot(self) -> dict[str, Any]:
        return {
            "units": self.units,
            "full_rhs_solves": self.full_rhs_solves,
            "rhs_by_kind": dict(self.rhs_by_kind),
            "full_factorizations": self.full_factorizations,
            "full_operator_products": self.full_operator_products,
            "basis_products_nr": self.basis_products_nr,
            "reduced_solves": self.reduced_solves,
            "reduced_factorizations": self.reduced_factorizations,
            "svd_count": self.svd_count,
            "rrqr_count": self.rrqr_count,
            "exact_acceptance_checks": self.exact_acceptance_checks,
            "wall_seconds": self.wall_seconds,
        }


@dataclass
class Calibration:
    """Time ratios fixed on tuning seeds; conservative factors included."""

    nominal_n: int = 256
    full_rhs_seconds: float = 0.001
    lu_factor_seconds: float = 0.001
    dense_matvec_seconds: float = 0.001
    reduced_qr_samples: dict[int, float] = field(default_factory=dict)
    reduced_solve_samples: dict[int, float] = field(default_factory=dict)
    svd_samples: dict[int, float] = field(default_factory=dict)  # keyed rows
    rrqr_samples: dict[int, float] = field(default_factory=dict)  # keyed k
    overhead_factor: float = 1.15
    raw_samples: dict[str, Any] = field(default_factory=dict)

    def lu_ratio(self) -> float:
        return (
            self.overhead_factor
            * self.lu_factor_seconds
            / max(self.full_rhs_seconds, 1e-12)
        )

    def matvec_ratio(self) -> float:
        return (
            self.overhead_factor
            * self.dense_matvec_seconds
            / max(self.full_rhs_seconds, 1e-12)
        )

    def reduced_qr_ratio(self, r: int) -> float:
        r = int(r)
        if not self.reduced_qr_samples:
            # flop estimate for an economic QR of an n x r matrix compared
            # with one full n x n triangular RHS solve (~n^2 flops)
            return 2.0 * r * r / float(self.nominal_n)
        t = _interp_table(self.reduced_qr_samples, r, self.nominal_n)
        return self.overhead_factor * t / max(self.full_rhs_seconds, 1e-12)

    def reduced_solve_ratio(self, r: int) -> float:
        r = int(r)
        if not self.reduced_solve_samples:
            return (r / float(self.nominal_n)) ** 2
        t = _interp_table(self.reduced_solve_samples, r, self.nominal_n)
        return self.overhead_factor * t / max(self.full_rhs_seconds, 1e-12)

    def svd_ratio(self, rows: int) -> float:
        rows = int(rows)
        keys = sorted(self.svd_samples)
        if not keys:
            return 0.0
        t = _interp_x(self.svd_samples, rows)
        return self.overhead_factor * t / max(self.full_rhs_seconds, 1e-12)

    def to_dict(self) -> dict[str, Any]:
        return {
            "nominal_n": self.nominal_n,
            "full_rhs_seconds": self.full_rhs_seconds,
            "lu_factor_seconds": self.lu_factor_seconds,
            "dense_matvec_seconds": self.dense_matvec_seconds,
            "reduced_qr_samples": {str(k): v for k, v in self.reduced_qr_samples.items()},
            "reduced_solve_samples": {
                str(k): v for k, v in self.reduced_solve_samples.items()
            },
            "svd_samples": {str(k): v for k, v in self.svd_samples.items()},
            "rrqr_samples": {str(k): v for k, v in self.rrqr_samples.items()},
            "overhead_factor": self.overhead_factor,
            "raw_samples": self.raw_samples,
        }

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "Calibration":
        return cls(
            nominal_n=int(d.get("nominal_n", 256)),
            full_rhs_seconds=float(d.get("full_rhs_seconds")),
            lu_factor_seconds=float(d.get("lu_factor_seconds")),
            dense_matvec_seconds=float(d.get("dense_matvec_seconds")),
            reduced_qr_samples=_as_int_dict(d.get("reduced_qr_samples", {})),
            reduced_solve_samples=_as_int_dict(d.get("reduced_solve_samples", {})),
            svd_samples=_as_int_dict(d.get("svd_samples", {})),
            rrqr_samples=_as_int_dict(d.get("rrqr_samples", {})),
            overhead_factor=float(d.get("overhead_factor", 1.15)),
            raw_samples=d.get("raw_samples", {}),
        )


def _as_int_dict(d: Any) -> dict[int, float]:
    return {int(k): float(v) for k, v in d.items()}


def _interp_table(table: dict[int, float], r: int, n: int) -> float:
    if r <= 0:
        return 0.0
    r = min(r, n)
    keys = sorted(table)
    if not keys:
        # fall back to a physical (r/n)^2 scaling of a full RHS
        return max(1e-12, (r / float(n)) ** 2)
    if r <= keys[0]:
        return table[keys[0]] * max(1.0, r / float(keys[0]))
    for lo, hi in zip(keys, keys[1:]):
        if lo <= r <= hi:
            a = (r - lo) / float(hi - lo)
            return math.exp((1 - a) * math.log(table[lo]) + a * math.log(table[hi]))
    return table[keys[-1]] * (r / float(keys[-1])) ** 2


def _interp_x(table: dict[int, float], x: int) -> float:
    keys = sorted(table)
    if not keys:
        return 0.0
    if x <= keys[0]:
        return table[keys[0]]
    for lo, hi in zip(keys, keys[1:]):
        if lo <= x <= hi:
            a = (x - lo) / float(hi - lo)
            return table[lo] + a * (table[hi] - table[lo])
    return table[keys[-1]]


class BudgetExhausted(Exception):
    """Raised before a full evaluation when the run work cap is reached."""


def calibration_path() -> Path:
    return SOLVER_DIR / "calibration_tune.json"


def frozen_timing_path(frozen_dir: str | Path) -> Path:
    return Path(frozen_dir) / "timing_frozen.json"


def write_json(path: str | Path, obj: Any) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(obj, indent=2, sort_keys=True, default=_jsonable) + "\n",
        encoding="utf-8",
    )


def read_json(path: str | Path) -> Any:
    return json.loads(Path(path).read_text(encoding="utf-8"))
