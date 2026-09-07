"""Two figures for the Levy continuity theorem / tightness discussion.

Figure 1  levy_escape_panels.png
    mu_n = (1 - a) * mu_C + a * Unif[n, 2n]
    row 1  F_n on the compactified axis u = x / (1 + x)
    row 2  Re phi_n(t) on a fixed window, showing the spike at t = 0
    row 3  the same in s = n t, showing the universal boundary layer

Figure 2  pisot_panels.png
    row 1  ||theta^m||  (distance to the nearest integer), Pisot vs non-Pisot
    row 2  upper envelope of |phi_xi(t)| for the Bernoulli convolution
           with ratio xi, xi = 1/3 (1/xi = 3 is Pisot) vs xi = 2/5

Run from the repository root:
    python scripts/levy_examples.py
"""

import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PURPLE, VIOLET, INK = "#57068C", "#8900E1", "#330662"
GRAY, CORAL, TEAL, AMBER = "#81758A", "#C5283D", "#00856A", "#E6A11A"


# --------------------------------------------------------------------------
# Cantor distribution
# --------------------------------------------------------------------------

def cantor_cdf(x, depth=34):
    """Cantor function K(x) = P(C <= x), C = 2 * sum_{k>=1} B_k 3^{-k}.

    Ternary scan: digit 0 -> binary bit 0, digit 2 -> bit 1,
    first digit 1 -> bit 1 and stop.
    """
    x = np.asarray(x, dtype=float)
    out = np.zeros_like(x)
    out[x >= 1] = 1.0
    done = (x <= 0) | (x >= 1)
    frac = np.clip(x, 0.0, 1.0).copy()
    scale = np.full_like(x, 0.5)
    for _ in range(depth):
        frac = frac * 3.0
        d = np.floor(frac)
        frac = frac - d
        hit1 = (d == 1) & (~done)
        out[hit1] += scale[hit1]
        done = done | hit1
        hit2 = (d == 2) & (~done)
        out[hit2] += scale[hit2]
        scale = scale / 2.0
    return out


def re_phi_cantor(t, depth=38):
    """Re phi_C(t) = cos(t/2) * prod_{k>=1} cos(t 3^{-k})."""
    t = np.asarray(t, dtype=float)
    p = np.cos(t / 2.0)
    f = t / 3.0
    for _ in range(depth):
        p = p * np.cos(f)
        f = f / 3.0
    return p


# --------------------------------------------------------------------------
# Escaping uniform slab and the mixture
# --------------------------------------------------------------------------

def unif_cdf(x, n):
    """CDF of Unif[n, 2n]."""
    return np.clip((np.asarray(x, dtype=float) - n) / n, 0.0, 1.0)


def re_phi_unif(t, n):
    """Re of the CF of Unif[n, 2n]:  (sin 2nt - sin nt) / (nt),  = 1 at t = 0."""
    t = np.asarray(t, dtype=float)
    out = np.ones_like(t)
    m = np.abs(t) > 1e-12
    out[m] = (np.sin(2 * n * t[m]) - np.sin(n * t[m])) / (n * t[m])
    return out


def F_mix(x, a, n):
    return (1 - a) * cantor_cdf(x) + a * unif_cdf(x, n)


def re_phi_mix(t, a, n):
    return (1 - a) * re_phi_cantor(t) + a * re_phi_unif(t, n)


def decompactify(u):
    """u in [0, 1] -> x in [0, inf];  u = x / (1 + x)."""
    u = np.asarray(u, dtype=float)
    return np.where(u >= 1, 1e12, u / np.maximum(1 - u, 1e-12))


# --------------------------------------------------------------------------
# Pisot side
# --------------------------------------------------------------------------

def dist_to_int_rational(p, q, M):
    """||(p/q)^m|| for m = 1..M, computed exactly with integers."""
    out = []
    for m in range(1, M + 1):
        num, den = p ** m, q ** m
        r = num % den
        f = r / den
        out.append(max(min(f, 1.0 - f), 1e-9))
    return np.array(out)


def dist_to_int_pisot(conj_abs, M):
    """Distance ||theta^m|| for a quadratic Pisot number.

    If theta^m + theta_conj^m is an integer, theta^m differs from that
    integer by |theta_conj|^m.  For the first few powers this error can exceed
    1/2, so the distance to the *nearest* integer is min(error, 1-error).
    """
    m = np.arange(1, M + 1)
    error = conj_abs ** m
    return np.minimum(error, np.abs(1.0 - error))


def cf_abs_bernoulli(t, xi, depth=55):
    """|phi_xi(t)| = prod_{k>=0} |cos(xi^k t)| for Y_xi = sum eps_k xi^k."""
    t = np.asarray(t, dtype=float)
    p = np.ones_like(t)
    f = t.copy()
    for _ in range(depth):
        p = p * np.abs(np.cos(f))
        f = f * xi
    return np.maximum(p, 1e-12)


def envelope(xi, t0, t1, nbins=340, nsub=260):
    """Upper envelope of |phi_xi|: max over each logarithmic bin.

    Necessary because the resonance peaks near t = 2 pi 3^n have width O(1)
    while sitting at t ~ 3^n, i.e. width ~3^{-n} in log coordinates.
    """
    edges = np.exp(np.linspace(np.log(t0), np.log(t1), nbins + 1))
    lo, hi = edges[:-1], edges[1:]
    frac = np.linspace(0.0, 1.0, nsub)[None, :]
    grid = lo[:, None] + (hi - lo)[:, None] * frac
    vals = cf_abs_bernoulli(grid.ravel(), xi).reshape(grid.shape)
    return np.sqrt(lo * hi), vals.max(axis=1)


# --------------------------------------------------------------------------
# Figure 1
# --------------------------------------------------------------------------

def fig_escape(a=0.35, ns=(1, 6, 40), path="levy_escape_panels.png"):
    fig, axes = plt.subplots(3, 3, figsize=(11.4, 8.4))
    u = np.linspace(0, 1, 4000)
    x = decompactify(u)
    Kx = cantor_cdf(x)
    t = np.linspace(-6, 6, 24000)
    phic = re_phi_cantor(t)
    s = np.linspace(-20, 20, 12000)
    xt = [0, 0.5, 2 / 3, 5 / 6, 1]
    xl = ["0", "1", "2", "5", r"$\infty$"]

    for j, n in enumerate(ns):
        ax = axes[0, j]
        ax.plot(u, (1 - a) * Kx, color=GRAY, lw=1.4, ls="--", label=r"$(1-a)K$")
        ax.plot(u, (1 - a) * Kx + a * unif_cdf(x, n), color=PURPLE, lw=1.8, label=r"$F_n$")
        ax.axhline(1 - a, color=CORAL, lw=0.8, ls=":")
        ax.set_xticks(xt)
        ax.set_xticklabels(xl)
        ax.set_ylim(-0.03, 1.06)
        ax.set_title(f"n = {n}", fontsize=12)
        ax.set_xlabel(r"$x$  (axis $u=x/(1+x)$)", fontsize=9)
        if j == 0:
            ax.set_ylabel(r"$F_n(x)$", fontsize=11)
            ax.legend(fontsize=9, loc="lower right", frameon=False)

        ax = axes[1, j]
        lim = (1 - a) * phic
        ax.plot(t, lim, color=GRAY, lw=1.4, ls="--")
        ax.plot(t, lim + a * re_phi_unif(t, n), color=PURPLE, lw=1.6)
        ax.plot([0], [1.0], "o", color=CORAL, ms=6)
        ax.annotate("", xy=(0, 1.0), xytext=(0, 1 - a),
                    arrowprops=dict(arrowstyle="<->", color=CORAL, lw=1.1))
        ax.text(0.35, 1 - a / 2, f"a = {a}", color=CORAL, fontsize=10, va="center")
        ax.set_xlim(-6, 6)
        ax.set_ylim(-0.55, 1.2)
        ax.set_xlabel(r"$t$", fontsize=9)
        if j == 0:
            ax.set_ylabel(r"$\mathrm{Re}\,\phi_n(t)$", fontsize=11)

        ax = axes[2, j]
        ax.plot(s, re_phi_mix(s / n, a, n), color=PURPLE, lw=1.6)
        ax.plot(s, (1 - a) + a * re_phi_unif(s, 1), color=CORAL, lw=1.0, ls="--")
        ax.axhline(1 - a, color=GRAY, lw=0.8, ls=":")
        ax.set_xlim(-20, 20)
        ax.set_ylim(-0.05, 1.15)
        ax.set_xlabel(r"$s = nt$", fontsize=9)
        if j == 0:
            ax.set_ylabel("rescaled", fontsize=11)

    _tidy(axes)
    fig.suptitle(
        r"Escaping mass $a$ = staircase deficit = jump of $\phi$ at $t=0$   "
        r"($\mu_n=(1-a)\mu_C+a\,\mathrm{Unif}[n,2n]$)", fontsize=12.5)
    fig.tight_layout(rect=[0, 0, 1, 0.955])
    path = Path(path)
    temporary = path.with_name(f".{path.stem}.tmp.png")
    fig.savefig(temporary, dpi=170)
    if not temporary.is_file() or temporary.stat().st_size == 0:
        raise RuntimeError(f"failed to render {path}")
    temporary.replace(path)
    plt.close(fig)
    return path


# --------------------------------------------------------------------------
# Figure 2
# --------------------------------------------------------------------------

def fig_pisot(M=22, path="pisot_panels.png"):
    fig, axes = plt.subplots(2, 1, figsize=(8.6, 7.2))
    m = np.arange(1, M + 1)

    ax = axes[0]
    ax.semilogy(m, dist_to_int_pisot((np.sqrt(5) - 1) / 2, M), "o-",
                color=PURPLE, lw=1.8, ms=4, label=r"$\varphi=(1+\sqrt{5})/2$  (Pisot)")
    ax.semilogy(m, dist_to_int_pisot(np.sqrt(2) - 1, M), "^--",
                color=TEAL, lw=1.8, ms=4, label=r"$1+\sqrt{2}=2.414$  (Pisot)")
    ax.semilogy(m, dist_to_int_rational(3, 2, M), "s-",
                color=CORAL, lw=1.6, ms=4, label=r"$3/2$")
    ax.semilogy(m, dist_to_int_rational(5, 2, M), "x--",
                color=AMBER, lw=1.6, ms=5, label=r"$5/2$")
    ax.set_xlabel(r"$m$", fontsize=10)
    ax.set_ylabel(r"$\|\theta^m\|$", fontsize=11)
    ax.set_ylim(1e-7, 1.0)
    ax.legend(fontsize=9, frameon=False, ncol=2)
    ax.set_title("distance from " + r"$\theta^m$" + " to the nearest integer", fontsize=11)

    ax = axes[1]
    for xi, col, lab, ls in ((1 / 3, PURPLE, r"$\xi=1/3$  envelope", "-"),
                             (0.4, CORAL, r"$\xi=2/5$  envelope", "--")):
        te, ve = envelope(xi, 2.0, 1.4e5)
        ax.loglog(te, ve, ls, color=col, lw=1.0, alpha=0.45, label=lab)
    res_t = np.array([2 * np.pi * 3.0 ** k for k in range(10)])
    res_t = res_t[res_t < 1.4e5]
    ax.loglog(res_t, cf_abs_bernoulli(res_t, 1 / 3), "o-", color=PURPLE, ms=6, lw=1.8,
              label=r"$|\phi_{1/3}(2\pi\cdot 3^{\,n})|$")
    res2 = np.array([2 * np.pi * 2.5 ** k for k in range(14)])
    res2 = res2[res2 < 1.4e5]
    ax.loglog(res2, cf_abs_bernoulli(res2, 0.4), "s--", color=AMBER, ms=5, lw=1.8,
              label=r"$|\phi_{2/5}(2\pi\cdot (5/2)^{\,n})|$")
    ax.set_xlabel(r"$t$", fontsize=10)
    ax.set_ylabel(r"$\sup|\phi_\xi|$ per bin", fontsize=11)
    ax.set_ylim(1e-5, 1.6)
    ax.legend(fontsize=9, frameon=False)
    ax.set_title("finite-range CF envelope and selected resonant scales", fontsize=11)

    _tidy(axes)
    fig.tight_layout()
    path = Path(path)
    temporary = path.with_name(f".{path.stem}.tmp.png")
    fig.savefig(temporary, dpi=170)
    if not temporary.is_file() or temporary.stat().st_size == 0:
        raise RuntimeError(f"failed to render {path}")
    temporary.replace(path)
    plt.close(fig)
    return path


def _tidy(axes):
    for ax in np.ravel(axes):
        ax.tick_params(labelsize=9)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.grid(alpha=0.25, lw=0.6)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "figures",
        help="directory for the two PNG figures (default: repository figures/)",
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    print(fig_escape(path=args.output_dir / "levy_escape_panels.png"))
    print(fig_pisot(path=args.output_dir / "pisot_panels.png"))
