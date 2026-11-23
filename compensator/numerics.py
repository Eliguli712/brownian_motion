"""Reproduce the original payoff experiments, tables and publication figures.

Numerical experiments are diagnostics; the analytic proof is in paper.tex.
All parameters below are illustrative choices, since the recorded problem
specifies arbitrary n, a_i and T. No fitted or measured data are used.
"""
from pathlib import Path
import json
import csv
import platform
import numpy as np
import scipy
from scipy.special import ndtr
from scipy.integrate import quad
from scipy.stats import norm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
FIG = ROOT / "figures"
DATA = ROOT / "data"
FIG.mkdir(exist_ok=True)
DATA.mkdir(exist_ok=True)
plt.rcParams.update({
    "font.family": "STIXGeneral", "mathtext.fontset": "stix", "font.size": 11,
    "axes.titlesize": 13, "axes.labelsize": 11, "legend.fontsize": 9,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": "white", "axes.facecolor": "white",
    "axes.grid": True, "grid.alpha": .17, "pdf.fonttype": 42,
    "savefig.bbox": "tight",
})
BLUE, TEAL, GOLD, GREY = "#174F78", "#138477", "#A47D26", "#536475"
T = 1.0
A = np.array([-.75, 0., .85])
N = len(A)
CUTS = ((1 - A[:, None] - A[None, :]) / 2).ravel()
PAIR = np.maximum(A[None, :] - A[:, None], 0).sum() / (2*N*N)
SEED, PATHS, STEPS = 20261123, 8192, 2048
INVROOT2PI = 1 / np.sqrt(2*np.pi)


def phi(y):
    return INVROOT2PI*np.exp(-.5*np.asarray(y)**2)


def payoff(x):
    x = np.asarray(x)
    return (-np.maximum(x[..., None] + A, 0).sum(axis=-1)/(4*N) + PAIR
            - (x[..., None] < CUTS).sum(axis=-1)/(4*N*N))


def semigroup(t, x):
    x = np.asarray(x)
    if t >= T:
        return payoff(x)
    sigma = np.sqrt(T-t)
    z = x[..., None] + A
    call = z*ndtr(z/sigma) + sigma*phi(z/sigma)
    return -call.sum(axis=-1)/(4*N) + PAIR - ndtr((CUTS-x[..., None])/sigma).sum(axis=-1)/(4*N*N)


def call_primitive(r, a, sigma):
    if sigma == 0:
        z = np.maximum(r + a, 0)
        return z**3/3 - a*z*z/2
    y = (r + a)/sigma
    p, f = ndtr(y), phi(y)
    return sigma**3/3*(y**3*p+(y*y-1)*f) - a*sigma**2/2*((y*y+1)*p+y*f)


def digital_primitive(r, c, sigma):
    if sigma == 0:
        return np.minimum(r, c)**2/2
    y = (c-r)/sigma
    p, f = ndtr(y), phi(y)
    i0 = y*p + f
    i1 = .5*((y*y-1)*p + y*f)
    return -c*sigma*i0 + sigma*sigma*i1


def primitive(t, x):
    x = np.asarray(x)
    sigma = np.sqrt(max(T-t, 0))
    calls = sum(call_primitive(x, a, sigma)-call_primitive(0., a, sigma) for a in A)
    digitals = sum(digital_primitive(x, c, sigma)-digital_primitive(0., c, sigma) for c in CUTS)
    return -calls/(4*N) + PAIR*x*x/2 - digitals/(4*N*N)


def derivatives(t, x):
    x = np.asarray(x)
    v = T-t
    sigma = np.sqrt(v)
    z = (x[..., None]+A)/sigma
    q = (CUTS-x[..., None])/sigma
    sx = -ndtr(z).sum(axis=-1)/(4*N) + phi(q).sum(axis=-1)/(4*N*N*sigma)
    sxx = -phi(z).sum(axis=-1)/(4*N*sigma) + (q*phi(q)).sum(axis=-1)/(4*N*N*v)
    st = phi(z).sum(axis=-1)/(8*N*sigma) - (q*phi(q)).sum(axis=-1)/(8*N*N*v)
    return st, sx, sxx


def save(fig, name):
    fig.savefig(FIG/(name+".pdf"))
    fig.savefig(FIG/(name+".png"), dpi=210)
    plt.close(fig)


def csv_out(name, headings, rows):
    with open(DATA/name, "w", newline="") as out:
        writer = csv.writer(out)
        writer.writerow(headings)
        writer.writerows(rows)


# Independent quadrature checks of the explicit spatial antiderivative.
qerrs, conv_errs = [], []
for t in [0., .3, .85, .995]:
    for x in [-2., -.5, .25, 1.5]:
        exact = quad(lambda r: r*float(semigroup(t, r)), 0, x, epsabs=1e-12, epsrel=1e-12)[0]
        qerrs.append(abs(float(primitive(t, x))-exact))
        sig = np.sqrt(T-t)
        knots = sorted(set([-float(a) for a in A] + list(CUTS)))
        knots = [-np.inf]+knots+[np.inf]
        conv = sum(quad(lambda y: float(payoff(y))*norm.pdf(y, loc=x, scale=sig),
                        left, right, epsabs=2e-12, epsrel=2e-12)[0]
                   for left, right in zip(knots[:-1], knots[1:]))
        conv_errs.append(abs(conv-float(semigroup(t, x))))
assert max(qerrs) < 2e-10, qerrs
assert max(conv_errs) < 2e-10, conv_errs

xs = np.linspace(-2.5, 2.5, 1001)
fig, ax = plt.subplots(figsize=(8.6, 4.6), layout="constrained")
ax.plot(xs, payoff(xs), c=GREY, lw=1.5, label=r"Terminal payoff $g(x)$")
for t, col in zip([0., .6, .95, .995], [BLUE, TEAL, GOLD, "#705F8D"]):
    ax.plot(xs, semigroup(t, xs), lw=1.9, color=col, label=fr"$S(t,x)$, $t={t:g}$")
ax.set(xlabel=r"State $x$", ylabel="Payoff / conditional value", title="The recorded nonsmooth payoff and Gaussian smoothing")
ax.legend(ncol=2)
save(fig, "01_payoff_and_semigroup")

tt = np.linspace(0, .995, 160)
xx = np.linspace(-2.2, 2.2, 230)
sm = np.array([semigroup(t, xx) for t in tt])
bm = sm-.5*np.array([semigroup(t, 0.) for t in tt])[:, None]
fig, axs = plt.subplots(1, 2, figsize=(10, 4.2), layout="constrained")
for ax, z, title in zip(axs, [sm, bm], [r"Heat solution $S(t,x)$", r"Compensator density $b(t,x)$"]):
    im = ax.pcolormesh(xx, tt, z, shading="auto", cmap="viridis")
    ax.set(xlabel=r"State $x$", ylabel=r"Time $t$", title=title)
    fig.colorbar(im, ax=ax, shrink=.87)
save(fig, "02_heat_and_drift")

# Central finite differences provide a numerical PDE test independent of the analytic derivative formulas.
pde_rows = []
for h in [0.02, 0.01, 0.005, 0.0025, 0.00125]:
    residuals = []
    for t in [.1, .4, .75, .9]:
        x = np.linspace(-1.8, 1.8, 181)
        dt = (semigroup(t+h, x)-semigroup(t-h, x))/(2*h)
        dxx = (semigroup(t, x+h)-2*semigroup(t, x)+semigroup(t, x-h))/(h*h)
        residuals.extend(dt+.5*dxx)
    pde_rows.append([h, float(np.max(np.abs(residuals))), float(np.sqrt(np.mean(np.square(residuals))))])
csv_out("pde_residual.csv", ["finite_difference_step", "max_abs_residual", "rms_residual"], pde_rows)

# Nested time grids share exactly the same fine Brownian sample paths.
rng = np.random.default_rng(SEED)
dt = T/STEPS
W = np.zeros((PATHS, STEPS+1))
W[:, 1:] = np.cumsum(rng.normal(0, np.sqrt(dt), (PATHS, STEPS)), axis=1)
times = np.linspace(0, T, STEPS+1)
P = np.empty_like(W)
K = np.zeros_like(W)
Ito = np.zeros_like(W)
density = np.empty_like(W)
for j, t in enumerate(times):
    sw = semigroup(t, W[:, j])
    density[:, j] = sw-.5*float(semigroup(t, 0.))
    P[:, j] = primitive(t, W[:, j])
    if j:
        K[:, j] = K[:, j-1]+.5*dt*(density[:, j-1]+density[:, j])
        Ito[:, j] = Ito[:, j-1]+W[:, j-1]*semigroup(times[j-1], W[:, j-1])*(W[:, j]-W[:, j-1])
M = P-K
meanM, seM = M.mean(axis=0), M.std(axis=0, ddof=1)/np.sqrt(PATHS)
residual = M-Ito
convergence = []
for nsteps in [64, 128, 256, 512, 1024, 2048]:
    stride = STEPS//nsteps
    w = W[:, ::stride]
    d = density[:, ::stride]
    kk = T/nsteps*(.5*d[:, 0]+d[:, 1:-1].sum(axis=1)+.5*d[:, -1])
    ii = np.zeros(PATHS)
    for j in range(nsteps):
        ii += w[:, j]*semigroup(T*j/nsteps, w[:, j])*(w[:, j+1]-w[:, j])
    rr = P[:, -1]-kk-ii
    convergence.append([nsteps, T/nsteps, float(np.sqrt(np.mean(rr*rr))), float(np.mean(rr))])
csv_out("ito_convergence.csv", ["steps", "dt", "rms_terminal_residual", "mean_terminal_residual"], convergence)

pick = 27
fig, axs = plt.subplots(2, 1, figsize=(8.7, 6.6), sharex=True, layout="constrained")
for a, lab, col in [(P, r"$P_t$", BLUE), (K, r"$K_t$", GOLD), (M, r"$M_t=P_t-K_t$", TEAL)]:
    axs[0].plot(times, a[pick], label=lab, color=col, lw=1.6)
axs[0].legend(ncol=3)
axs[0].set(ylabel="Process value", title="One Brownian path: functional, compensator and martingale")
axs[1].plot(times, residual[pick], color=GREY, lw=1.1)
axs[1].axhline(0, color="black", lw=.7)
axs[1].set(xlabel=r"Time $t$", ylabel=r"$P_t-K_t-I_t^{(h)}$", title="Discretization residual against the left-point stochastic sum")
save(fig, "03_path_decomposition")

fig, axs = plt.subplots(1, 2, figsize=(10, 4.2), layout="constrained")
axs[0].fill_between(times, meanM-1.96*seM, meanM+1.96*seM, color=TEAL, alpha=.18, label="Pointwise 95% MC intervals")
axs[0].plot(times, meanM, color=TEAL, label=r"Estimated $\mathbb{E}[M_t]$")
axs[0].axhline(0, color="black", ls="--", lw=1, label="Theorem: zero")
axs[0].set(xlabel="Time", ylabel="Ensemble mean", title="Mean-zero diagnostic")
axs[0].legend()
mid = STEPS//2
bins = np.quantile(W[:, mid], np.linspace(0, 1, 9))
cond_rows = []
for k in range(8):
    mask = (W[:, mid] >= bins[k]) & ((W[:, mid] < bins[k+1]) if k < 7 else (W[:, mid] <= bins[k+1]))
    inc = M[mask, -1]-M[mask, mid]
    cond_rows.append([float(W[mask, mid].mean()), int(mask.sum()), float(inc.mean()), float(inc.std(ddof=1)/np.sqrt(mask.sum()))])
axs[1].errorbar([r[0] for r in cond_rows], [r[2] for r in cond_rows], yerr=[1.96*r[3] for r in cond_rows], fmt="o", color=BLUE, capsize=3)
axs[1].axhline(0, color="black", ls="--", lw=1)
axs[1].set(xlabel=r"Bin mean of $W_{1/2}$", ylabel=r"Mean of $M_1-M_{1/2}$", title="A coarse conditional-mean diagnostic")
save(fig, "04_martingale_diagnostics")
csv_out("conditional_diagnostic.csv", ["bin_mean_state", "count", "mean_increment", "standard_error"], cond_rows)
csv_out("ensemble_timeseries.csv", ["t", "mean_P", "mean_K", "mean_M", "se_M"],
        zip(times[::16], P[:, ::16].mean(0), K[:, ::16].mean(0), meanM[::16], seM[::16]))

fig, axs = plt.subplots(1, 2, figsize=(10, 4.2), layout="constrained")
hs = np.array([r[1] for r in convergence])
errs = np.array([r[2] for r in convergence])
axs[0].loglog(hs, errs, "o-", color=BLUE, label="Nested-path RMS residual")
axs[0].loglog(hs, errs[-1]*np.sqrt(hs/hs[-1]), "--", color=GOLD, label=r"Reference $h^{1/2}$")
axs[0].set(xlabel="Time step h", ylabel="Terminal RMS error", title="Stochastic discretization convergence")
axs[0].legend()
axs[1].loglog([r[0] for r in pde_rows], [r[2] for r in pde_rows], "o-", color=TEAL, label="Centered-difference PDE residual")
axs[1].set(xlabel="Finite-difference step", ylabel="RMS residual", title=r"Checking $S_t+\frac{1}{2}S_{xx}=0$")
axs[1].legend()
save(fig, "05_numerical_convergence")

envelope = np.abs(A).sum()/(4*N)+PAIR+.25
constantC = envelope+.25*np.sqrt(2*T/np.pi)+.25
bound = constantC**2*(T*T+2*T**3)
s00 = float(semigroup(0., 0.))
expectedK = T*s00-.5*quad(lambda t: float(semigroup(t, 0.)), 0, T, epsabs=2e-12, epsrel=2e-12, points=[.99])[0]
knots = [-np.inf]+sorted(set([-float(a) for a in A]+list(CUTS)))+[np.inf]
expectedP = sum(quad(lambda x: float(primitive(T, x))*norm.pdf(x, scale=np.sqrt(T)), l, r,
                     epsabs=2e-12, epsrel=2e-12)[0] for l, r in zip(knots[:-1], knots[1:]))
assert abs(expectedP-expectedK) < 2e-10

# Deterministic Gaussian and time quadrature for E integral H_s^2 ds.
gnodes, gweights = np.polynomial.hermite.hermgauss(100)
lnodes, lweights = np.polynomial.legendre.leggauss(180)
ts = (lnodes+1)*T/2
energy = []
for t in ts:
    w = np.sqrt(2*t)*gnodes
    energy.append(float(np.dot(gweights, (w*semigroup(t, w))**2)/np.sqrt(np.pi)))
energy_integral = float(np.dot(lweights, energy)*T/2)
# Independent adaptive spatial quadrature, split at all payoff breakpoints.
def adaptive_energy(t):
    sd = np.sqrt(t)
    breaks = sorted(set([-float(a)/sd for a in A]+list(CUTS/sd)))
    limits = [-np.inf]+breaks+[np.inf]
    return sum(quad(lambda z: float((sd*z*semigroup(t, sd*z))**2)*norm.pdf(z),
                    l, r, epsabs=2e-10, epsrel=2e-9)[0]
               for l, r in zip(limits[:-1], limits[1:]))
energy_adaptive = float(np.dot(lweights, [adaptive_energy(t) for t in ts])*T/2)
fig, axs = plt.subplots(1, 2, figsize=(10, 4.2), layout="constrained")
axs[0].plot(ts, energy, c=BLUE, label=r"Gauss-Hermite: $\mathbb{E}[H_t^2]$")
axs[0].plot(ts, 2*constantC**2*(ts+3*ts**2), c=GOLD, label="Proved polynomial upper bound")
axs[0].set(xlabel="Time", ylabel="Energy density", title="Square-integrability and its explicit bound")
axs[0].legend()
axs[1].hist(M[:, -1], bins=85, color=TEAL, alpha=.85, density=True)
axs[1].axvline(0, color="black", lw=1, ls="--")
axs[1].set(xlabel=r"Terminal $M_T$", ylabel="Empirical density", title="A martingale can have a skewed terminal law")
save(fig, "06_energy_and_terminal_law")
csv_out("energy_density.csv", ["t", "quadrature_energy", "majorant"], zip(ts, energy, 2*constantC**2*(ts+3*ts**2)))

metrics = {
    "parameters": {"T":T, "n":N, "a":A.tolist(), "seed":SEED, "paths":PATHS, "fine_steps":STEPS},
    "software": {"python":platform.python_version(), "numpy":np.__version__, "scipy":scipy.__version__, "matplotlib":matplotlib.__version__},
    "max_primitive_vs_quad_error":max(qerrs), "max_semigroup_vs_quad_error":max(conv_errs),
    "payoff_envelope_A":float(envelope), "uniform_growth_C":float(constantC), "energy_upper_bound":float(bound),
    "S_0_0":s00, "deterministic_E_P_T":expectedP, "deterministic_E_K_T":expectedK,
    "deterministic_E_M_T":expectedP-expectedK,
    "quadrature_energy":energy_integral, "adaptive_spatial_energy":energy_adaptive,
    "MC_mean_P_T":float(P[:, -1].mean()), "MC_mean_K_T":float(K[:, -1].mean()),
    "MC_mean_M_T":float(meanM[-1]), "MC_se_M_T":float(seM[-1]),
    "MC_variance_M_T":float(M[:, -1].var(ddof=1)),
    "MC_second_moment_M_T":float(np.mean(M[:, -1]**2)),
    "MC_se_second_moment_M_T":float(np.std(M[:, -1]**2, ddof=1)/np.sqrt(PATHS)),
    "MC_mean_I_T":float(Ito[:, -1].mean()), "fine_rms_ito_residual":float(errs[-1]),
    "convergence_slope":float(np.polyfit(np.log(hs), np.log(errs), 1)[0]),
    "PDE_finest_rms":float(pde_rows[-1][2]),
    "convergence":convergence,
}
(DATA/"metrics.json").write_text(json.dumps(metrics, indent=2)+"\n")
np.savez_compressed(DATA/"sample_paths.npz", t=times, W=W[:32], P=P[:32], K=K[:32], M=M[:32], Ito=Ito[:32])
print(json.dumps(metrics, indent=2))
