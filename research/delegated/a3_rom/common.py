"""Shared deterministic configuration and audited cost accounting for the
A3 physical-ROM reconstruction worker.

This worker develops a *fixed-chart, residual-minimising* full-wave reduced
model and a joint material/pose solver whose trial acceptance always uses the
exact physical objective.  It reuses the A2 ``a2_physics`` core read-only:
no source file there is modified and no module-level mutation is performed.
All runs here are DEVELOPMENT runs on seeds 4101-4104 only.
"""

from __future__ import annotations

import hashlib
import json
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
DIR = Path(__file__).resolve().parent
A2_PHYSICS = ROOT / "research" / "delegated" / "a2_physics"
if str(A2_PHYSICS) not in sys.path:
    sys.path.insert(0, str(A2_PHYSICS))

import physics as A2  # noqa: E402  read-only import

# ---------------------------------------------------------------------------
# Fixed development configuration
# ---------------------------------------------------------------------------

DEVELOPMENT_SEEDS = tuple(range(4101, 4105))
FORBIDDEN_FINAL_SEEDS = tuple(range(1001, 1021))
KMAX = 12.0
FREQ_IDS = (0, 1, 3)  # k = 3, 6, 12 rad/m under Config(kmax=12)
KS = tuple(
    float(k) for k in (KMAX * np.array([0.25, 0.5, 0.75, 1.0]))[list(FREQ_IDS)]
)
POSE_LEVER = 1.5
LAMBDA_MIN = 2.0 * np.pi / KMAX
POSE_SCALE = np.array(
    [1.0 / LAMBDA_MIN, 1.0 / LAMBDA_MIN, POSE_LEVER / LAMBDA_MIN]
)
SNR_DB = 30.0
BG_ALPHA = 0.08
ALPHA_BOX = (0.0, 2.5)
X_TRANS_BOX = (-0.7, 0.7)
X_ROT_BOX = (-0.7, 0.7)
ROWS_PER_FREQ = 72  # 3 poses x 2 illuminations x 12 receivers


def guard_seed(seed: int) -> None:
    if int(seed) in FORBIDDEN_FINAL_SEEDS:
        raise RuntimeError(
            f"seed {seed} is in the forbidden E4 final range; never use it"
        )
    if int(seed) not in DEVELOPMENT_SEEDS:
        raise RuntimeError(
            f"seed {seed} is not in the predeclared development range "
            f"{DEVELOPMENT_SEEDS}"
        )


def make_model(
    n: int = 16,
    aperture: str = "full",
    kmax: float = KMAX,
    material_basis: np.ndarray | None = None,
) -> A2.Model:
    if int(n) not in (8, 16):
        raise ValueError("only grid8 smoke and grid16 development are allowed")
    cfg = A2.Config(
        N=int(n),
        aperture=aperture,
        kmax=float(kmax),
        material_basis=material_basis,
    )
    return A2.Model(cfg)


def config_json(n: int = 16) -> dict[str, Any]:
    return {
        "grid": int(n),
        "aperture": "full",
        "kmax": KMAX,
        "frequency_k_rad_per_m": list(map(float, KS)),
        "freq_ids_in_model": list(map(int, FREQ_IDS)),
        "pose_lever_m": POSE_LEVER,
        "lambda_min_m": float(LAMBDA_MIN),
        "pose_scale": list(map(float, POSE_SCALE)),
        "material_q_default": 9,
        "noise_snr_db": SNR_DB,
        "snr_reference": "grid16_background_alpha_0.08",
    }


# ---------------------------------------------------------------------------
# JSON helpers
# ---------------------------------------------------------------------------


def _jsonable(v: Any) -> Any:
    if isinstance(v, np.ndarray):
        return v.tolist()
    if isinstance(v, (np.integer, np.floating)):
        return float(v)
    if isinstance(v, complex):
        return [float(v.real), float(v.imag)]
    if isinstance(v, tuple):
        return list(v)
    if isinstance(v, Path):
        return str(v)
    raise TypeError(f"not JSON serialisable: {type(v)!r}")


def jsonable(v: Any) -> Any:
    return json.loads(json.dumps(v, default=_jsonable))


def settings_hash(content: Any) -> str:
    b = json.dumps(jsonable(content), sort_keys=True, separators=(",", ":")).encode(
        "utf-8"
    )
    return hashlib.sha256(b).hexdigest()


def write_json(path: str | Path, obj: Any) -> None:
    Path(path).write_text(
        json.dumps(obj, indent=2, sort_keys=True, default=_jsonable) + "\n",
        encoding="utf-8",
    )


# ---------------------------------------------------------------------------
# Parameter scaling and physics objective helpers
# ---------------------------------------------------------------------------


def scaled_from_params(alpha: np.ndarray, x: np.ndarray) -> np.ndarray:
    return np.concatenate(
        [np.asarray(alpha, dtype=float), np.asarray(x, dtype=float) * POSE_SCALE]
    )


def params_from_scaled(z: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    z = np.asarray(z, dtype=float)
    return z[:9].copy(), z[9:] / POSE_SCALE


def z_bounds() -> list[tuple[float, float]]:
    bounds = [tuple(ALPHA_BOX)] * 9
    bounds.append((X_TRANS_BOX[0] * POSE_SCALE[0], X_TRANS_BOX[1] * POSE_SCALE[0]))
    bounds.append((X_TRANS_BOX[0] * POSE_SCALE[1], X_TRANS_BOX[1] * POSE_SCALE[1]))
    bounds.append((X_ROT_BOX[0] * POSE_SCALE[2], X_ROT_BOX[1] * POSE_SCALE[2]))
    return bounds


def select_rows(y: np.ndarray, freq_ids: tuple[int, ...]) -> np.ndarray:
    y = np.asarray(y, dtype=np.complex128).reshape(-1)
    ids = sorted(set(int(f) for f in freq_ids))
    return np.concatenate([y[int(f) * ROWS_PER_FREQ : (int(f) + 1) * ROWS_PER_FREQ] for f in ids])


def realify(z: np.ndarray, sigma: float = 1.0) -> np.ndarray:
    return A2.realify(np.asarray(z), float(sigma))


def exact_loss(res_complex: np.ndarray, sigma: float) -> float:
    return 0.5 * float(np.sum(np.abs(realify(res_complex, sigma)) ** 2))


def lever_metric_error(x1: np.ndarray, x2: np.ndarray) -> float:
    d = np.asarray(x1, dtype=float) - np.asarray(x2, dtype=float)
    return float(np.sqrt(d[0] ** 2 + d[1] ** 2 + POSE_LEVER**2 * d[2] ** 2))


# ---------------------------------------------------------------------------
# Work ledger
# ---------------------------------------------------------------------------


@dataclass
class WorkLedger:
    """Transparent accounting for full/reduced linear algebra.

    ``full_rhs_columns`` counts dense M^{-1} LU RHS columns and mirrors the
    A2 core's ``rhs_solves_total``.  Reduced triangular solves are reported
    separately, never converted to full-RHS-equivalent units.  LU, QR and SVD
    factorisations, matrix products (including their column widths), basis
    products and wall time are each explicit.
    """

    full_rhs_columns: int = 0
    full_rhs_by_kind: dict[str, int] = field(default_factory=dict)
    lu_factorizations: int = 0
    qr_factorizations: int = 0
    svd_factorizations: int = 0
    reduced_tri_rhs_columns: int = 0
    full_operator_products: int = 0
    basis_products: int = 0
    product_columns: int = 0
    matrix_dimensions: list[list[int]] = field(default_factory=list)
    exact_acceptance_evals: int = 0
    wall_seconds: float = 0.0

    def charge_full(
        self,
        rhs: int = 0,
        kind: str = "other",
        lu: int = 0,
        products: int = 0,
        product_columns: int = 0,
        wall: float = 0.0,
        dims: list[int] | tuple[int, ...] | None = None,
    ) -> None:
        self.full_rhs_columns += int(rhs)
        self.full_rhs_by_kind[kind] = self.full_rhs_by_kind.get(kind, 0) + int(rhs)
        self.lu_factorizations += int(lu)
        self.full_operator_products += int(products)
        self.product_columns += int(product_columns)
        self.wall_seconds += float(wall)
        if lu:
            self.matrix_dimensions.append([int(dims[0]) if dims else -1, int(dims[1]) if dims else -1])

    def charge_model_work(self, work: dict[str, Any]) -> None:
        self.full_rhs_columns += int(work.get("rhs_solves_total", 0))
        for k, v in (work.get("rhs_solves") or {}).items():
            self.full_rhs_by_kind[str(k)] = (
                self.full_rhs_by_kind.get(str(k), 0) + int(v)
            )
        self.lu_factorizations += int(work.get("factorizations", 0))
        self.full_operator_products += int(work.get("operator_products", 0))
        self.wall_seconds += float(work.get("wall_seconds", 0.0))

    def charge_reduced(
        self,
        tri_rhs: int = 0,
        qr: int = 0,
        svd: int = 0,
        products: int = 0,
        basis: int = 0,
        product_columns: int = 0,
        wall: float = 0.0,
        dims: list[int] | tuple[int, ...] | None = None,
    ) -> None:
        self.reduced_tri_rhs_columns += int(tri_rhs)
        self.qr_factorizations += int(qr)
        self.svd_factorizations += int(svd)
        self.full_operator_products += int(products)
        self.basis_products += int(basis)
        self.product_columns += int(product_columns)
        self.wall_seconds += float(wall)
        if qr or svd:
            self.matrix_dimensions.append(
                [int(dims[0]) if dims else -1, int(dims[1]) if dims else -1]
            )

    def charge_acceptance(self, n: int = 1) -> None:
        self.exact_acceptance_evals += int(n)

    def add_wall(self, seconds: float) -> None:
        self.wall_seconds += float(seconds)

    def snapshot(self) -> dict[str, Any]:
        return {
            "full_rhs_columns": self.full_rhs_columns,
            "full_rhs_by_kind": dict(self.full_rhs_by_kind),
            "lu_factorizations": self.lu_factorizations,
            "qr_factorizations": self.qr_factorizations,
            "svd_factorizations": self.svd_factorizations,
            "reduced_tri_rhs_columns": self.reduced_tri_rhs_columns,
            "full_operator_products": self.full_operator_products,
            "basis_products": self.basis_products,
            "product_columns": self.product_columns,
            "matrix_dimensions": self.matrix_dimensions,
            "exact_acceptance_evals": self.exact_acceptance_evals,
            "wall_seconds": self.wall_seconds,
        }


# ---------------------------------------------------------------------------
# Scenes / data (development seeds only)
# ---------------------------------------------------------------------------


def _noise_draw(rng: np.random.Generator, shape: tuple[int, ...], sigma: float):
    return (float(sigma) / np.sqrt(2.0)) * (
        rng.standard_normal(shape) + 1j * rng.standard_normal(shape)
    )


def nominal_sigma() -> float:
    """30 dB total-field noise scale fixed from grid16 alpha=0.08 reference."""
    model = make_model(n=16)
    fw = model.forward(
        np.full(9, BG_ALPHA), np.zeros(3), list(FREQ_IDS), jacobian=False
    )
    power = float(np.mean(np.abs(fw["total"]) ** 2))
    return float(np.sqrt(power / 10.0 ** (SNR_DB / 10.0)))


def make_alpha_true(seed: int) -> np.ndarray:
    guard_seed(seed)
    rng = np.random.default_rng(seed)
    alpha = np.full(9, BG_ALPHA)
    n_targets = int(rng.integers(2, 4))
    idx = rng.choice(9, size=n_targets, replace=False)
    strengths = rng.uniform(0.35, 0.90, size=n_targets)
    for i, s in zip(idx, strengths):
        alpha[int(i)] += float(s)
    return alpha


def pose_init(seed: int) -> np.ndarray:
    """Small shared-pose initialisation error (translation + rotation).

    Development nonlinear acceptance is checked at a small (<=~0.08 m,
    <=~0.05 rad) initial offset.  Larger-offset basin comparisons are a
    separate later experiment, not this development foundation.
    """
    guard_seed(seed)
    rng = np.random.default_rng(seed + 10_000)
    angle = rng.uniform(0.0, 2.0 * np.pi)
    radius = rng.uniform(0.04, 0.08)
    dtheta = rng.uniform(0.01, 0.05) * (1.0 if rng.random() < 0.5 else -1.0)
    return np.array(
        [radius * np.cos(angle), radius * np.sin(angle), dtheta], dtype=float
    )


def make_scene(seed: int, sigma: float | None = None) -> dict[str, Any]:
    guard_seed(seed)
    sigma = nominal_sigma() if sigma is None else float(sigma)
    model = make_model(n=16)
    alpha_true = make_alpha_true(seed)
    # Full four-frequency raw observations are retained so every stage can
    # address a model frequency index directly with select_rows().
    fw = model.forward(alpha_true, np.zeros(3), jacobian=False)
    rng = np.random.default_rng(seed)
    y = fw["total"] + _noise_draw(rng, fw["total"].shape, sigma)
    rng_val = np.random.default_rng(seed + 50_000)
    y_val = fw["total"] + _noise_draw(rng_val, fw["total"].shape, sigma)
    return {
        "seed": int(seed),
        "alpha_true": alpha_true.tolist(),
        "x_true": [0.0, 0.0, 0.0],
        "y": jsonable(y),
        "y_val": jsonable(y_val),
        "sigma": float(sigma),
        "y_true": jsonable(fw["total"]),
    }


def scene_arrays(scene: dict[str, Any]) -> dict[str, Any]:
    return {
        "alpha_true": np.asarray(scene["alpha_true"], dtype=float),
        "x_true": np.asarray(scene["x_true"], dtype=float),
        "y": complex_from_json(scene["y"]),
        "y_val": complex_from_json(scene["y_val"]),
        "y_true": complex_from_json(scene["y_true"]),
        "sigma": float(scene["sigma"]),
    }


def complex_from_json(v: Any) -> np.ndarray:
    a = np.asarray(v)
    if a.ndim == 2 and a.shape[1] == 2:
        return (a[:, 0] + 1j * a[:, 1]).astype(np.complex128)
    return a.astype(np.complex128)


def utcnow_iso() -> str:
    from datetime import datetime, timezone

    return datetime.now(timezone.utc).isoformat(timespec="seconds")
