"""Shared deterministic settings, scene/data generation and N20 adapter for
the exploratory higher-dimensional 2D imaging comparison (a2_highdim).

This worker module is deliberately self-contained and makes no research
claims.  The physical full-wave core is imported read-only from
``research/delegated/a2_physics/physics.py``.  Physics ``Config`` normally
restricts the grid to N in {8,16,32}; this exploratory task predeclares an
N=20 inverse grid, so this module applies a runtime-only extension of the
allowed-N tuple (``_ALLOWED_N += (20,)``) before any N=20 ``Config`` is
constructed.  No physics source file is modified and no physical formula is
changed.

Single-thread CPU only.  Do not access E4 final seeds 1001-1020.
"""

from __future__ import annotations

import hashlib
import json
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
DIR = Path(__file__).resolve().parent
PHYS_DIR = ROOT / "research" / "delegated" / "a2_physics"

if str(PHYS_DIR) not in sys.path:
    sys.path.insert(0, str(PHYS_DIR))

import physics  # noqa: E402  (read-only import after path setup)

# Runtime-only N20 extension of the validation tuple (see module docstring).
if 20 not in tuple(int(v) for v in physics._ALLOWED_N):  # type: ignore[attr-defined]
    physics._ALLOWED_N = tuple(physics._ALLOWED_N) + (20,)  # type: ignore[attr-defined]

# ---------------------------------------------------------------------------
# Fixed physical / experimental design constants
# ---------------------------------------------------------------------------

KMAX = 4.0 * np.pi
KS = KMAX * np.array([0.25, 0.5, 0.75, 1.0])
LAMBDA_MIN = 2.0 * np.pi / KMAX  # 0.5 m
POSE_LEVER = 1.5
POSE_SCALE = np.array(
    [1.0 / LAMBDA_MIN, 1.0 / LAMBDA_MIN, POSE_LEVER / LAMBDA_MIN]
)  # z = x * POSE_SCALE; grad_z = grad_x / POSE_SCALE

N_INV = 20
N_DATA = 32
N_ALPHA = 49
BASIS_CENTRES = tuple(float(v) for v in np.linspace(-0.35, 0.35, 7))
BASIS_WIDTH = 0.10
BG_ALPHA = 0.08

ALPHA_BOX = (0.0, 2.5)
X_BOX_TRANSLATION = (-0.7, 0.7)
X_BOX_ROTATION = (-0.7, 0.7)

SNR_DB = 30.0
SQRT2 = float(np.sqrt(2.0))

FINAL_E4_SEEDS = range(1001, 1021)
TUNING_SEEDS = (71, 72)
TEST_SEEDS = tuple(range(801, 807))
MISMATCH_SEED = 881

# Same low-to-high cumulative frequency schedule for every method.  Stage s
# uses frequencies KS[: s+1].  Stage ``maxiter`` is modest and frozen; the
# global per-run RHS budget is a hard cap for all methods.
STAGE_FREQS: tuple[tuple[int, ...], ...] = ((0,), (0, 1), (0, 1, 2), (0, 1, 2, 3))
STAGE_MAXITER_DEFAULT = (18, 18, 14, 18)
BUDGET_RHS_DEFAULT = 5000
MAXLS_DEFAULT = 12


def default_settings() -> dict[str, Any]:
    return {
        "version": 1,
        "label": "a2_highdim_exploratory",
        "budget_rhs_columns": BUDGET_RHS_DEFAULT,
        "stage_freqs": [list(s) for s in STAGE_FREQS],
        "stage_maxiter": list(STAGE_MAXITER_DEFAULT),
        "maxls": MAXLS_DEFAULT,
        "ftol": 1e-10,
        "gtol": 1e-7,
        "prior_strength": 0.0,
        "prior_center": [0.0] * N_ALPHA,
        "alpha_box": list(ALPHA_BOX),
        "x_box_translation_m": list(X_BOX_TRANSLATION),
        "x_box_rotation_rad": list(X_BOX_ROTATION),
        "alpha_init": BG_ALPHA,
        "pose_scale": list(map(float, POSE_SCALE)),
        "pose_lever_m": POSE_LEVER,
        "lambda_min_m": float(LAMBDA_MIN),
        "kmax": float(KMAX),
        "frequencies_rad_per_m": list(map(float, KS)),
        "snr_db": SNR_DB,
        "noise_reference": "nominal_n32_background_alpha_0_08",
        "inverse_grid": N_INV,
        "data_grid": N_DATA,
        "material_basis": {
            "type": "gaussian_unnormalised",
            "centres_min": -0.35,
            "centres_max": 0.35,
            "centres_per_axis": 7,
            "width_std_m": BASIS_WIDTH,
        },
    }


def canonical_bytes(obj: Any) -> bytes:
    return json.dumps(
        obj, sort_keys=True, separators=(",", ":"), default=_jsonable
    ).encode("utf-8")


def settings_hash(settings: dict[str, Any]) -> str:
    content = {k: v for k, v in settings.items() if k not in ("hash", "timestamp")}
    return hashlib.sha256(canonical_bytes(content)).hexdigest()


def _jsonable(v: Any) -> Any:
    if isinstance(v, np.ndarray):
        return v.tolist()
    if isinstance(v, (np.integer, np.floating)):
        return float(v)
    if isinstance(v, tuple):
        return list(v)
    if isinstance(v, Path):
        return str(v)
    if isinstance(v, complex):
        return [v.real, v.imag]
    raise TypeError(f"not JSON-serialisable: {type(v)}")


def jsonable(v: Any) -> Any:
    return json.loads(json.dumps(v, default=_jsonable))


def guard_seed(seed: int, allow_final: bool = False) -> None:
    if seed in FINAL_E4_SEEDS and not allow_final:
        raise RuntimeError(
            f"seed {seed} is in the E4 final range 1001-1020; this "
            "exploratory worker must never touch it"
        )


def write_frozen(settings: dict[str, Any], path: Path | None = None) -> dict[str, Any]:
    """Write the one frozen settings file with a hash and timestamp."""
    path = path or (DIR / "frozen.json")
    if path.exists():
        raise FileExistsError(
            f"frozen settings already exist at {path}; refusing to overwrite"
        )
    content = jsonable(settings)
    digest = settings_hash(content)
    frozen = {
        "settings": content,
        "sha256": digest,
        "frozen_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "frozen_by": "a2_highdim worker",
    }
    path.write_text(
        json.dumps(frozen, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return frozen


def load_frozen(path: Path | None = None) -> dict[str, Any]:
    path = path or (DIR / "frozen.json")
    data = json.loads(path.read_text(encoding="utf-8"))
    settings = dict(data["settings"])
    digest = settings_hash(settings)
    if digest != data["sha256"]:
        raise RuntimeError(
            f"frozen settings hash mismatch: file says {data['sha256']}, "
            f"recomputed {digest}"
        )
    return settings


# ---------------------------------------------------------------------------
# Models and grids
# ---------------------------------------------------------------------------


def gaussian_basis_for_grid(
    points: np.ndarray, width: float = BASIS_WIDTH
) -> np.ndarray:
    cols = []
    for cx in BASIS_CENTRES:
        for cy in BASIS_CENTRES:
            c = np.array([float(cx), float(cy)], dtype=float)
            d = points - c[None, :]
            cols.append(np.exp(-np.einsum("ij,ij->i", d, d) / (2.0 * width**2)))
    return np.stack(cols, axis=1)


_model_cache: dict[str, physics.Model] = {}


def make_model(N: int, basis: str = "gaussian49") -> physics.Model:
    """Cached N20/N32 forward model with the same analytic 49-column basis.

    ``basis="pixel32"`` constructs an N=32 identity basis (used only for the
    outside-49-basis mismatch truth where the true spatial map is defined
    pixel-wise).
    """
    key = f"{N}:{basis}"
    if key in _model_cache:
        return _model_cache[key]
    if basis == "gaussian49":
        points, _ = physics.make_grid(N)
        Phi = gaussian_basis_for_grid(points)
        cfg = physics.Config(N=N, material_basis=Phi)
        mdl = physics.Model(cfg)
    elif basis == "pixel32":
        if N != 32:
            raise ValueError("pixel basis is only predeclared for N=32")
        n = N * N
        mdl = physics.Model(physics.Config(N=N, material_basis=np.eye(n)))
    else:
        raise ValueError(f"unknown basis {basis!r}")
    _model_cache[key] = mdl
    return mdl


def make_grid_points(N: int) -> tuple[np.ndarray, float]:
    return physics.make_grid(int(N))


def spatial_chi(map_points: np.ndarray, alpha: np.ndarray) -> np.ndarray:
    Phi = gaussian_basis_for_grid(map_points)
    return Phi @ np.asarray(alpha, dtype=float).reshape(-1)


# ---------------------------------------------------------------------------
# Reference noise, scenes and pose initialisation
# ---------------------------------------------------------------------------

_sigma_cache: float | None = None


def nominal_sigma() -> float:
    """30 dB total-field noise scale fixed once from the declared nominal
    N32 reference material (uniform alpha=0.08 background)."""
    global _sigma_cache
    if _sigma_cache is None:
        model = make_model(N_DATA, "gaussian49")
        alpha_ref = np.full(N_ALPHA, BG_ALPHA)
        fw = model.forward(alpha_ref, np.zeros(3), jacobian=False)
        power = float(np.mean(np.abs(fw["total"]) ** 2))
        _sigma_cache = float(np.sqrt(power / 10.0 ** (SNR_DB / 10.0)))
    return _sigma_cache


def make_alpha_true(seed: int) -> np.ndarray:
    """Nonuniform positive coefficient scene: declared background alpha 0.08
    plus 2-3 smooth Gaussian target coefficients.  Deterministic per seed;
    no E4 final seeds are ever used."""
    guard_seed(seed)
    rng = np.random.default_rng(seed)
    alpha = np.full(N_ALPHA, BG_ALPHA)
    n_targets = int(rng.integers(2, 4))  # 2 or 3
    idx = rng.choice(N_ALPHA, size=n_targets, replace=False)
    strengths = rng.uniform(0.35, 0.90, size=n_targets)
    for i, s in zip(idx, strengths):
        alpha[int(i)] += float(s)
    return alpha


def _noise_draw(rng: np.random.Generator, shape: tuple[int, ...], sigma: float):
    return (sigma / SQRT2) * (
        rng.standard_normal(shape) + 1j * rng.standard_normal(shape)
    )


def make_scene(seed: int) -> dict[str, Any]:
    """Coherent total-field observations plus an independent noisy replicate.
    Data are generated with the N32 model on the analytic 49-Gaussian basis.
    The same y/y_val are shared by every method and both pose-error radii."""
    guard_seed(seed)
    sigma = nominal_sigma()
    model = make_model(N_DATA, "gaussian49")
    alpha_true = make_alpha_true(seed)
    fw = model.forward(alpha_true, np.zeros(3), jacobian=False)
    y_true = fw["total"]
    rng = np.random.default_rng(seed)
    noise = _noise_draw(rng, y_true.shape, sigma)
    y = y_true + noise
    rng_val = np.random.default_rng(10_000 + seed)
    y_val = y_true + _noise_draw(rng_val, y_true.shape, sigma)
    return {
        "seed": int(seed),
        "alpha_true": alpha_true.tolist(),
        "x_true": [0.0, 0.0, 0.0],
        "y_true": _complex_to_json(y_true),
        "y": _complex_to_json(y),
        "y_val": _complex_to_json(y_val),
        "sigma": float(sigma),
    }


def scene_arrays(scene: dict[str, Any]) -> dict[str, Any]:
    return {
        "alpha_true": np.asarray(scene["alpha_true"], dtype=float),
        "x_true": np.asarray(scene["x_true"], dtype=float),
        "y_true": np.asarray(scene["y_true"], dtype=complex),
        "y": np.asarray(scene["y"], dtype=complex),
        "y_val": np.asarray(scene["y_val"], dtype=complex),
        "sigma": float(scene["sigma"]),
    }


def _complex_to_json(z: np.ndarray) -> list[list[float]]:
    return [[float(v.real), float(v.imag)] for v in np.asarray(z).reshape(-1)]


def complex_from_json(v: Any) -> np.ndarray:
    """Decode the [re, im] JSON encoding back to a complex vector."""
    a = np.asarray(v)
    if a.ndim == 2 and a.shape[1] == 2:
        return (a[:, 0] + 1j * a[:, 1]).astype(np.complex128)
    if a.ndim == 0:
        raise ValueError("complex scalar JSON is not used by this worker")
    return a.astype(np.complex128)


def pose_error_radius(seed: int, radius_index: int) -> dict[str, Any]:
    """Predeclared 12-combination init table used by tuning/test runs.

    For each (seed,radius) pair k=(seed-base)*2+radius_index cycles a
    translation angle around the circle; the squared pose metric is split
    equally between translation and lever-arm rotation, alternating the
    rotation sign as in the A2 protocol (exploratory, not E4 primary).
    """
    if seed in TUNING_SEEDS:
        base = int(seed) - 71
        if radius_index not in (0,):
            raise ValueError("tuning seeds use only the .125 lambda radius")
        r = 0.125 * LAMBDA_MIN
        k = base * 2
    else:
        base = int(seed) - TEST_SEEDS[0]
        r = (0.125 if radius_index == 0 else 0.5) * LAMBDA_MIN
        k = base * 2 + radius_index
    trans = r / SQRT2
    dtheta = r / (SQRT2 * POSE_LEVER)
    angle = k * np.pi / 6.0
    sign = 1.0 if k % 2 == 0 else -1.0
    x0 = [trans * np.cos(angle), trans * np.sin(angle), sign * dtheta]
    return {
        "radius_m": float(r),
        "radius_label": f"{r / LAMBDA_MIN:.3f}lambda",
        "angle_rad": float(angle),
        "x0": list(map(float, x0)),
        "seed": int(seed),
        "radius_index": int(radius_index),
    }


def mismatch_init_x0() -> np.ndarray:
    """Predeclared .125-lambda error for mismatch seed 881: translation along
    +x and positive lever-arm rotation (same split as the init table)."""
    r = 0.125 * LAMBDA_MIN
    return np.array(
        [r / SQRT2, 0.0, r / (SQRT2 * POSE_LEVER)], dtype=float
    )


def lever_metric_error(x1: np.ndarray, x2: np.ndarray) -> float:
    d = np.asarray(x1, dtype=float) - np.asarray(x2, dtype=float)
    return float(
        np.sqrt(d[0] ** 2 + d[1] ** 2 + POSE_LEVER**2 * d[2] ** 2)
    )


def scaled_from_params(alpha: np.ndarray, x: np.ndarray) -> np.ndarray:
    return np.concatenate(
        [np.asarray(alpha, dtype=float), np.asarray(x, dtype=float) * POSE_SCALE]
    )


def params_from_scaled(z: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    z = np.asarray(z, dtype=float)
    return z[:N_ALPHA].copy(), z[N_ALPHA:] / POSE_SCALE


def z_bounds() -> list[tuple[float, float]]:
    b = [tuple(ALPHA_BOX)] * N_ALPHA
    b.append(
        (
            X_BOX_TRANSLATION[0] * POSE_SCALE[0],
            X_BOX_TRANSLATION[1] * POSE_SCALE[0],
        )
    )
    b.append(
        (
            X_BOX_TRANSLATION[0] * POSE_SCALE[1],
            X_BOX_TRANSLATION[1] * POSE_SCALE[1],
        )
    )
    b.append(
        (
            X_BOX_ROTATION[0] * POSE_SCALE[2],
            X_BOX_ROTATION[1] * POSE_SCALE[2],
        )
    )
    return b


@dataclass
class CostLedger:
    """One unit = one full-wave LU RHS solve column.  Factorisations,
    operator products and wall seconds are reported separately and never
    converted into solve units."""

    rhs_solves: int = 0
    rhs_by_kind: dict[str, int] = None  # type: ignore[assignment]
    factorizations: int = 0
    operator_products: int = 0
    wall_seconds: float = 0.0

    def __post_init__(self) -> None:
        if self.rhs_by_kind is None:
            self.rhs_by_kind = {}

    def charge(self, rhs: int = 0, kind: str = "other", factorizations: int = 0,
               operator_products: int = 0, wall: float = 0.0) -> None:
        self.rhs_solves += int(rhs)
        self.rhs_by_kind[kind] = self.rhs_by_kind.get(kind, 0) + int(rhs)
        self.factorizations += int(factorizations)
        self.operator_products += int(operator_products)
        self.wall_seconds += float(wall)

    def charge_model_work(self, work: dict[str, Any]) -> None:
        self.rhs_solves += int(work.get("rhs_solves_total", 0))
        for k, v in (work.get("rhs_solves") or {}).items():
            self.rhs_by_kind[k] = self.rhs_by_kind.get(k, 0) + int(v)
        self.factorizations += int(work.get("factorizations", 0))
        self.operator_products += int(work.get("operator_products", 0))
        self.wall_seconds += float(work.get("wall_seconds", 0.0))

    def add_wall(self, seconds: float) -> None:
        self.wall_seconds += float(seconds)

    @property
    def units(self) -> int:
        return self.rhs_solves

    def snapshot(self) -> dict[str, Any]:
        return {
            "units_rhs_columns": self.rhs_solves,
            "rhs_by_kind": dict(self.rhs_by_kind),
            "factorizations": self.factorizations,
            "operator_products": self.operator_products,
            "wall_seconds": self.wall_seconds,
        }
