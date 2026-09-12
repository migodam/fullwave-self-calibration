"""Physical and projected-state interfaces for the A12 HG graph package.

This is a standalone 2-D scalar VIE implementation on the specified 0.4 m
domain.  It intentionally keeps material rendering and HG current ports
separate.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import math
import numpy as np
from scipy.special import hankel1

C0 = 299_792_458.0
FREQ_GHZ = np.array([1.5, 2.25, 3.0])
FREF_GHZ = 2.25


@dataclass(frozen=True)
class Geometry:
    n: int = 32
    side: float = 0.4
    n_tx: int = 8
    n_rx: int = 32
    radius_tx: float = 0.34
    radius_rx: float = 0.30


def grid(g: Geometry) -> tuple[np.ndarray, float]:
    h = g.side / g.n
    x = -g.side / 2 + (np.arange(g.n) + .5) * h
    xx, yy = np.meshgrid(x, x, indexing="ij")
    return np.c_[xx.ravel(), yy.ravel()], h


def _green(src: np.ndarray, dst: np.ndarray, k: float) -> np.ndarray:
    d = dst[:, None] - src[None, :]
    r = np.sqrt(np.sum(d * d, axis=-1))
    out = np.zeros(r.shape, np.complex128)
    nz = r > 0
    out[nz] = .25j * hankel1(0, k * r[nz])
    return out


def _self(k: float, h: float) -> complex:
    a = h / math.sqrt(math.pi)
    return .5j * math.pi * a / k * hankel1(1, k * a) - 1 / k**2


def _ports_for_patch(points: np.ndarray, center: np.ndarray, sigma: float, count: int, mask: np.ndarray) -> np.ndarray:
    x = (points[:, 0] - center[0]) / sigma
    y = (points[:, 1] - center[1]) / sigma
    base = np.exp(-.5 * (x * x + y * y))
    def hermite(n, z):
        h0=np.ones_like(z)
        if n==0: return h0
        h1=z
        if n==1: return h1
        for k in range(2,n+1): h0,h1=h1,z*h1-(k-1)*h0
        return h1
    # Preserve the already-frozen H0/H1/H2/H3 column order for 1/3/6/10
    # ports, then append higher total-degree modes.  This keeps stored R96
    # coordinates stable while permitting rank-aware higher-order closure.
    pairs=[(0,0),(1,0),(0,1),(2,0),(0,2),(1,1),(3,0),(0,3),(2,1),(1,2)]
    degree=4
    while len(pairs)<count:
        pairs.extend((nx,degree-nx) for nx in range(degree+1)); degree+=1
    cols=[hermite(nx,x)*hermite(ny,y)*base for nx,ny in pairs[:count]]
    raw = np.column_stack(cols[:count])
    raw[~mask] = 0.0
    # QR is discrete Gram whitening under the grid's Euclidean quadrature;
    # h is a common factor and does not alter the orthonormal span.
    q, _ = np.linalg.qr(raw, mode="reduced")
    return q[:, :count]


def hg_basis(g: Geometry, ports_per_patch: int = 6) -> tuple[np.ndarray, list[np.ndarray]]:
    if ports_per_patch not in (1, 3, 6, 10, 15, 21, 28, 36):
        raise ValueError("unsupported HG total-degree port count")
    pts, _ = grid(g)
    patch_side = g.side / 4
    sigma = .36 * patch_side
    cols, index = [], []
    for ix in range(4):
        for iy in range(4):
            c = np.array([-g.side/2 + (ix+.5)*patch_side,
                          -g.side/2 + (iy+.5)*patch_side])
            lo = -g.side/2 + np.array([ix, iy]) * patch_side
            hi = lo + patch_side
            mask = ((pts[:, 0] >= lo[0]) & (pts[:, 0] < hi[0]) &
                    (pts[:, 1] >= lo[1]) & (pts[:, 1] < hi[1]))
            q = _ports_for_patch(pts, c, sigma, ports_per_patch, mask)
            start = len(cols)
            cols.extend(q[:, j] for j in range(ports_per_patch))
            index.append(np.arange(start, start + ports_per_patch))
    return np.column_stack(cols).astype(np.complex128), index


class Operators:
    """Cached dense physical operators and HG projection operators."""
    def __init__(self, geom: Geometry = Geometry(), ports_per_patch: int = 6, aperture: str = "full"):
        self.geom, self.points, self.h = geom, *grid(geom)
        self.Q, self.node_index = hg_basis(geom, ports_per_patch)
        self.ports_per_patch = ports_per_patch
        self.freq_hz = FREQ_GHZ * 1e9
        tx_a = np.linspace(0, 2*np.pi, geom.n_tx, endpoint=False)
        rx_a = (np.linspace(-np.pi/2, np.pi/2, geom.n_rx) if aperture == "half"
                else np.linspace(0, 2*np.pi, geom.n_rx, endpoint=False))
        self.tx = geom.radius_tx * np.c_[np.cos(tx_a), np.sin(tx_a)]
        self.rx = geom.radius_rx * np.c_[np.cos(rx_a), np.sin(rx_a)]
        self.D, self.E, self.S, self.direct = [], [], [], []
        for f in self.freq_hz:
            k = 2*np.pi*f/C0
            d = k*k*self.h*self.h*_green(self.points, self.points, k)
            np.fill_diagonal(d, k*k*_self(k, self.h))
            e = _green(self.tx, self.points, k)  # [N,Tx]
            s = k*k*self.h*self.h*_green(self.points, self.rx, k)  # [Rx,N]
            direct = _green(self.tx, self.rx, k)  # [Rx,Tx]
            self.D.append(d); self.E.append(e); self.S.append(s); self.direct.append(direct)
        # Material-independent contraction used by every scene.  Caching it
        # keeps per-scene data cost to chi weighting + small reduced algebra.
        self.DQ = [d @ self.Q for d in self.D]

    def chi_frequency(self, chi_ref: np.ndarray, fi: int) -> np.ndarray:
        return np.asarray(chi_ref, float) * (1 + .02j * FREF_GHZ / FREQ_GHZ[fi])

    def reduced_system(self, chi_ref: np.ndarray, fi: int) -> tuple[np.ndarray, np.ndarray]:
        chi = self.chi_frequency(chi_ref, fi)
        dq = self.DQ[fi]
        h = self.Q.conj().T @ (chi[:, None] * dq)
        f = self.Q.conj().T @ (chi[:, None] * self.E[fi])
        return h, f


def material_render(theta: list[dict] | np.ndarray, points: np.ndarray, soft_partition: bool = False) -> np.ndarray:
    """Render additive Gaussian material or the sharp OOD soft partition."""
    if soft_partition:
        # theta=[cx,cy,wx,wy,amplitude,softness] with smooth signed edges.
        cx, cy, wx, wy, a, s = np.asarray(theta, float)
        x, y = points[:, 0], points[:, 1]
        return a * .25 * (1 + np.tanh((wx-abs(x-cx))/s)) * (1 + np.tanh((wy-abs(y-cy))/s))
    out = np.zeros(points.shape[0], float)
    for q in theta:
        a, cx, cy, sx, sy, ang = (q[k] for k in ("a","cx","cy","sx","sy","angle"))
        c, s = np.cos(ang), np.sin(ang)
        d = points - np.array([cx, cy])
        u = c*d[:, 0] + s*d[:, 1]
        v = -s*d[:, 0] + c*d[:, 1]
        out += a * np.exp(-.5*((u/sx)**2 + (v/sy)**2))
    return out


def random_material(rng: np.random.Generator, points: np.ndarray, family: str = "gaussian", regime: str | None = None) -> tuple[np.ndarray, object]:
    regime = regime or str(rng.choice(["weak","medium","strong"],p=[.34,.34,.32]))
    scale = {"weak":1.0,"medium":4.0,"strong":20.0}[regime]
    if family == "sharp":
        t = np.array([rng.uniform(-.06,.06), rng.uniform(-.06,.06), rng.uniform(.045,.10), rng.uniform(.04,.09), scale*rng.uniform(.025,.050), rng.uniform(.006,.014)])
        return material_render(t, points, True), {"family":"sharp", "regime":regime, "scale":scale, "theta":t.tolist()}
    p = int(rng.integers(2, 6)); terms=[]
    for _ in range(p):
        terms.append({"a":float(scale*rng.uniform(.015,.045)), "cx":float(rng.uniform(-.12,.12)), "cy":float(rng.uniform(-.12,.12)), "sx":float(rng.uniform(.025,.06)), "sy":float(rng.uniform(.025,.06)), "angle":float(rng.uniform(-np.pi,np.pi))})
    return material_render(terms, points), {"family":"gaussian", "regime":regime, "scale":scale, "theta":terms}


def forward_full(ops: Operators, chi_ref: np.ndarray, return_current: bool = False) -> dict:
    currents=[]; fields=[]; scattered=[]
    for fi in range(len(FREQ_GHZ)):
        chi=ops.chi_frequency(chi_ref, fi)
        m=np.eye(ops.points.shape[0], dtype=np.complex128)-chi[:,None]*ops.D[fi]
        j=np.linalg.solve(m, chi[:,None]*ops.E[fi])
        currents.append(j.T)
        sca=ops.S[fi]@j; scattered.append(sca.T); fields.append((ops.direct[fi]+sca).T)
    out={"field":np.asarray(fields),"scattered":np.asarray(scattered)} # [F,Tx,Rx]
    if return_current: out["current"]=np.asarray(currents) # [F,Tx,N]
    return out


def forward_reduced(ops: Operators, chi_ref: np.ndarray) -> dict:
    bs=[]; fields=[]; scattered=[]; hs=[]; fs=[]
    for fi in range(len(FREQ_GHZ)):
        h, f=ops.reduced_system(chi_ref,fi)
        b=np.linalg.solve(np.eye(h.shape[0],dtype=complex)-h,f)
        chi=ops.chi_frequency(chi_ref,fi)
        c=ops.S[fi]@ops.Q
        sca=c@b; scattered.append(sca.T); fields.append((ops.direct[fi]+sca).T); bs.append(b.T); hs.append(h); fs.append(f.T)
    return {"b":np.asarray(bs),"field":np.asarray(fields),"scattered":np.asarray(scattered),"H":np.asarray(hs),"f":np.asarray(fs)}
