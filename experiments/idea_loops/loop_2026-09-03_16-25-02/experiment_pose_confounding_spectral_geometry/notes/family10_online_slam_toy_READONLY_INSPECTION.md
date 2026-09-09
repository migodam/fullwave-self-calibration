# Family 10 online-SLAM toy: read-only code inspection

Source inspected: `src/family10_online_slam_toy.py` (1,298 lines).
The source was **not modified** and the experiment was **not run**.
This note is an annotated, faithful condensation of the requested parts; line
numbers refer to the source file.

## 1. `CONFIG` (lines 89-188) and trial-count logic

```python
CONFIG = {
    "N": 16, "T": 6, "n_rx": 4,
    "m_c": 24,        # T * n_rx complex rows per frequency
    "m_real": 144,    # 2 * m_c * len(frequencies) after stacking
    "p": 24,          # smooth basis dimension
    "q_pose": 18,     # 3 * T pose parameters
    "arc_radius": 1.6,
    "arc_phi_deg": [-45.0, 45.0],
    "pose_theta_convention": "theta = atan2(-p_y, -p_x): body +x axis points toward the origin",
    "rx_offsets": [[-0.06,0.0],[0.06,0.0],[0.0,-0.06],[0.0,0.06]],
    "tx_offset": [0.0,0.0],
    "chi_blobs": {"amp1":0.3,"amp2":0.5,"sigma1":0.09,"sigma2":0.07,
                  "c1":[-0.15,-0.12],"c2":[0.18,0.14]},
    "smooth_basis": {"p":24,"x_centers_n":4,"y_centers_n":6,
                     "x_span":[-0.3,0.3],"y_span":[-0.3,0.3],
                     "sigma_b":0.16,"unit_columns":True},
    "k_b_rule": "k_b(f) = 2*pi*f with f in {1.0, 1.4, 1.8}",
    "frequencies": [1.0, 1.4, 1.8],
    "snr": 100.0,
    "alpha": 1.0,
    "noise_model": {
        "per_frequency_complex_covariance":
            "sigma_f^2 I, sigma_f^2 = ||A_pix||_F^2 / (m_c*snr)",
        "whitening": "scalar division of complex blocks and forward by 1/sigma_f",
        "realification":
            "A_R = sqrt(2)*[Re(A_w); Im(A_w)], B_R = sqrt(2)*[Re(B_w); Im(B_w)] "
            "via hh.whiten_realify(A_w, B_w, None); "
            "y_R = sqrt(2)*[Re(F_w); Im(F_w)]",
        "stacked_noise_covariance":
            "I_m_real after realification (independent standard Gaussians on each stacked row)",
        "noise_draw":
            "one default_rng(20260903); the same per-trial draw is applied to "
            "both Born and full-wave true data (matched comparison)",
    },
    "monte_carlo": {"n_trials": 200, "seed": 20260903,
                    "modes": ["born", "full_wave"]},
    "fit": {"method":"trf","x_scale":"jac","max_nfev":1200,
            "xtol":1e-10,"ftol":1e-10,"gtol":1e-10,
            "note": "local fits start at the true parameters; x_scale='jac' ..."},
    "n_confounded_directions": 3,
    "finite_difference_self_check": {"eps":1e-6,"seed":20260903,"mode":"born"},
    "coefficient_choice": "c0 = argmin_c ||S c - chi0||_2 (numpy lstsq); chi_true = S c0",
    "scope_note": "finite-dimensional toy only; no continuum/global/real-SLAM claim",
}
```

Actual n-trials used by the run is decided in `main()` (lines 1041+): the CLI
arguments `--n-trials`, `--max-nfev`, and `--seed` default to `None` and
overwrite `CONFIG["monte_carlo"]["n_trials"]`,
`CONFIG["fit"]["max_nfev"]`, and `CONFIG["monte_carlo"]["seed"]` when supplied.

## 2. Whitened/realified stacks (per frequency and across frequencies)

Constant and helper (lines 190, 216-220):

```python
_SQRT2 = float(np.sqrt(2.0))

def realify_vector(z):            # exact vector convention
    return _SQRT2 * np.concatenate([np.real(z), np.imag(z)])
```

The reused module realifies matrices identically (`helmholtz.py`,
`whiten_realify`, lines 518-531):

```python
# A_hat = W A, B_hat = W B; here W=None so A_hat=A_w, B_hat=B_w.
A_R = _SQRT2 * np.vstack([np.real(A_hat), np.imag(A_hat)])
B_R = _SQRT2 * np.vstack([np.real(B_hat), np.imag(B_hat)])
```

Core stack builder (`build_model_blocks`, lines 299-360):

```python
for f in cfg["frequencies"]:
    k_b = 2.0 * np.pi * float(f)
    if mode == "full_wave":
        A_c, B_c, F_c, _ = hh.build_AB(chi, poses, rx, tx, N, k_b)
    else:  # "born"
        F_c, A_c, B_c = born_forward_AB(chi, poses, rx, tx, N, k_b)

    sigma2 = float(np.linalg.norm(A_c, ord="fro") ** 2) / (float(m_c) * snr)
    sigma  = float(np.sqrt(sigma2))

    A_w = A_c / sigma
    B_w = B_c / sigma
    F_w = F_c / sigma

    A_pix_R, B_R = hh.whiten_realify(A_w, B_w, None)   # x sqrt(2)
    A_s = A_pix_R @ scene["S"]                          # map-coef Jacobian
    y_R = realify_vector(F_w)                           # x sqrt(2)

    As.append(A_s); Bs.append(B_R); ys.append(y_R); sigmas.append(sigma)

A_stack = np.vstack(As)          # (2*m_c*n_freq, p) = (144, 24)
B_stack = np.vstack(Bs)          # (144, 3*T)       = (144, 18)
y_stack = np.concatenate(ys)     # (144,)
```

The Born helper (`born_forward_AB`, lines 246-298) is only present because no
Born pose Jacobian exists in the reused modules. For each pose row it builds

```python
G_S  = (k_b**2) * (h**2) * hh.green_matrix(points, rx_t, k_b)
E_inc = hh.green_matrix(tx_t[None, :], points, k_b)[:, 0]
A[sl, :] = G_S * E_inc[None, :]
F[sl]    = G_S @ (chi * E_inc)
# per pose-motion component l in {0,1,2}:
DXGS_l = (k_b**2) * (h**2) * np.einsum("asd,ad->as", grad1, drx_t[l])
DXE_l  = np.einsum("sd,d->s", grad2, dtx_t[l])
B[sl, 3*t+l] = (DXGS_l * E_inc[None, :] + G_S * DXE_l[None, :]) @ chi
```

where `grad1 = green_grad_first(...)` is the receiver-side derivative and
`grad2 = green_grad_source(...)` the source-motion derivative used by
`drx_t`/`dtx_t`.

## 3. Residual/Jacobian passed to `least_squares`

`ResidualModel._forward_at(theta, free)` (lines 416-492) is the exact API used
by the `residual`/`jacobian` callbacks:

```python
p   = S.shape[1]; q = x0.size
c   = theta[:p]
dx  = theta[p:] if free else np.zeros(q)
chi = S @ c
poses = (x0 + dx).reshape(T, 3)

for f, sigma in zip(frequencies, self._sigmas):
    k_b = 2*pi*f
    if mode == "full_wave": A_c, B_c, F_c, _ = hh.build_AB(chi, poses, rx, tx, N, k_b)
    else:                   F_c, A_c, B_c = born_forward_AB(chi, poses, rx, tx, N, k_b)
    A_w, B_w, F_w = A_c/sigma, B_c/sigma, F_c/sigma
    A_pix_R, B_R = hh.whiten_realify(A_w, B_w, None)
    As.append(A_pix_R @ S); Bs.append(B_R); ms.append(realify_vector(F_w))

A = vstack(As); B = vstack(Bs); model = concatenate(ms)

if free:
    res = np.concatenate([self.target - model, np.sqrt(alpha) * dx])
    J_data  = np.hstack([-A, -B])
    J_prior = np.hstack([np.zeros((q, p)), np.sqrt(alpha) * np.eye(q)])
    J = np.vstack([J_data, J_prior])
else:
    res = self.target - model
    J = -A
```

So `chi = S c` is formed *before* entering the forward/Jacobian API, the
residual always is `target - forward_model`, and the analytic data Jacobian is
the negative of the realified/whitened derivative blocks. The pose prior is
exactly `sqrt(alpha)*dx` extra rows in raw pose units, giving the
`B^T B + alpha I` information contribution. Results are cached per `(free,
theta)`.

Optimization coordinates are preconditioned (`PreconditionedModel`, lines
500-533; `whitening_sqrt_root`, lines 493-499): the optimizer variable is `y`
with `theta = theta0 + T y` and `T` chosen so the true-point information metric
is `I` (`T^T K_IS T = I` known-pose, `T^T H_joint T = I` free-pose). The
`least_squares` calls (`fit_least_squares`, lines 534-568) are

```python
sol = least_squares(
    lambda th: model.residual(th, free),
    np.zeros(precond_y_size),            # start at y=0 => true theta
    jac=lambda th: model.jacobian(th, free),
    method="trf", x_scale="jac", max_nfev=1200,
    xtol=1e-10, ftol=1e-10, gtol=1e-10)
theta_hat = model.theta_hat(sol.x)
```

`PreconditionedModel.jacobian(y, free)` returns `J_theta @ T`; recovery is
`theta0 + T @ y`, so `sol.x` is not the estimate returned.

## 4. Theory: K_IS, K_eff, P_known, P_free, and generalized eigenpairs

`linearized_predictions` (lines 368-415), with
`shrinkage_W` (lines 361-366):

```python
def shrinkage_W(B, alpha):
    C = B.T @ B + alpha * np.eye(q)
    return np.eye(B.shape[0]) - B @ np.linalg.solve(C, B.T)

K_IS    = sym(A_stack.T @ A_stack)
W       = shrinkage_W(B_stack, alpha)
K_eff   = sym(A_stack.T @ (W @ A_stack))
P_known = np.linalg.inv(K_IS)
P_free  = np.linalg.inv(K_eff)

# K_IS = L L^T
L = np.linalg.cholesky(sym(K_IS))
X = solve_triangular(L, K_eff, lower=True)              # L^{-1} K_eff
M = sym(solve_triangular(L, X.T, lower=True).T)          # L^{-1} K_eff L^{-T}
eigvals_u, U = np.linalg.eigh(M)
rho_asc = eigvals_u[np.argsort(eigvals_u)]
V = solve_triangular(L.T, U[:, order], lower=False)      # V = L^{-T} U
V = V / np.linalg.norm(V, axis=0)[None, :]               # unit Euclidean norm
retained_dof = float(rho_asc.sum())
```

This is the generalized pair `K_eff v = rho * K_IS v`, with `rho` ascending;
`retained_dof = tr(K_eff K_IS^{-1})`. Symmetry guards, eigenvalue/condition
diagnostics, and inverse covariances are also recorded.

For the free-pose preconditioner and consistency check the same joint
information matrix is used (inside `run_monte_carlo`):

```python
H_joint = block([[K_IS, A_stack.T @ B_stack],
                 [(A_stack.T @ B_stack).T, B_stack.T @ B_stack + alpha*eye(q)]])
```

## 5. Monte Carlo loop

One noise matrix is drawn once in `main()` and shared by both modes
(lines 1170-1180 approximately):

```python
rng = np.random.default_rng(int(cfg["monte_carlo"]["seed"]))  # seed 20260903
noise = rng.normal(size=(n_trials, m_real))                   # (200, 144)
```

where `m_real = blocks["born"]["y_stack"].size` (= 144). Per trial
(`run_monte_carlo`, lines 629-783):

```python
for i in range(n_trials):
    target = blocks["y_stack"] + noise[i]
    base_model = ResidualModel(scene, cfg, blocks, target=target)
    model_k = PreconditionedModel(base_model, c0, T_known)       # known pose
    model_f = PreconditionedModel(base_model,
                                  np.concatenate([c0, zeros(q)]),
                                  T_free)                        # free pose
    sk = fit_least_squares(model_k, free=False, precond_y_size=p)
    sf = fit_least_squares(model_f, free=True,  precond_y_size=p+q)
    ok_known[i], ok_free[i] = sk["success"], sf["success"]
    errs_known[i] = sk["theta_hat"][:p] - c0
    errs_free[i]  = sf["theta_hat"][:p] - c0
```

Note both fits start at the truth (`c0`, `dx=0`); they only retain the map
coefficient block `theta_hat[:p]` as the error.

Empirical covariances and metrics:

```python
def emp_cov(err, mask):
    rows = err[mask]
    if rows.shape[0] < 2: return full((p, p), nan)
    return np.cov(rows, rowvar=False, bias=False)   # sample covariance, /(n-1)

Cov_known = emp_cov(errs_known, ok_known)
Cov_free  = emp_cov(errs_free,  ok_free)
mean_known = errs_known[ok_known].mean(axis=0) if n_ok_known else full(p, nan)
mean_free  = errs_free[ok_free].mean(axis=0)     if n_ok_free  else full(p, nan)
```

Comparison metrics (`covariance_metrics`, lines 569-596), when finite:

```python
rel_fro      = ||Cov_emp - P_pred||_F / ||P_pred||_F
rel_spectral = ||Cov_emp - P_pred||_2 / max(||P_pred||_2, 0.0)
emp_eigvals_desc = sort(eigvalsh(Cov_emp))[::-1]
pred_eigvals_desc = sort(eigvalsh(P_pred))[::-1]
```

For the `n_confounded_directions = 3` smallest-rho directions `v = V[:, i]`,
each row records:

```python
pred_known_var = v @ P_known @ v
pred_free_var  = v @ P_free  @ v
emp_known_var  = v @ Cov_known @ v
emp_free_var   = v @ Cov_free  @ v
predicted_inflation_1_over_rho = 1/rho
predicted_var_ratio_exact = pred_free_var / pred_known_var
empirical_inflation_ratio = emp_free_var / emp_known_var
```

Trace ratios use `tr(P_free)/tr(P_known)` and
`tr(Cov_free)/tr(Cov_known)` with zero-denominator guards. Trial-level success,
status, nfev, cost, message, `c_error_l2`, and (free) `dx_l2` are stored per
row; the raw `errs_*`/`ok_*` arrays are retained in memory and excluded from
the JSON.
