#!/usr/bin/env python3
"""Visual explanation of Brownian density via Fourier inversion and contour shift.

The script produces one six-panel figure.  Its central complex-analysis step is
the contour-shift identity for the entire function

    g_t(z) = exp(-t z**2 / 2).

Because g_t has no singularities, the residue sum is zero.  The argument is
therefore best understood as Cauchy contour deformation (or the zero-residue
case of the residue theorem), not as a non-zero residue calculation.

Run from the repository root:
    python scripts/brownian_contour.py

Output:
    brownian_contour_visualization.png
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch


OUTPUT = (
    Path(__file__).resolve().parents[1]
    / "figures"
    / "brownian_contour_visualization.png"
)


def heat_kernel(x: np.ndarray, t: float) -> np.ndarray:
    """The one-dimensional Brownian transition density p_t(x)."""
    if t <= 0:
        raise ValueError("t must be strictly positive")
    return np.exp(-(x**2) / (2.0 * t)) / np.sqrt(2.0 * np.pi * t)


def truncated_fourier_inverse(
    x: np.ndarray, t: float, cutoff: float, n_omega: int = 6001
) -> np.ndarray:
    """Numerically invert exp(-t omega^2/2) over [-cutoff, cutoff]."""
    omega = np.linspace(-cutoff, cutoff, n_omega)
    integrand = np.exp(-0.5 * t * omega[None, :] ** 2) * np.exp(
        -1j * x[:, None] * omega[None, :]
    )
    return np.trapezoid(integrand, omega, axis=1).real / (2.0 * np.pi)


def vertical_edge_integral(R: float, t: float, h: float, n_y: int = 4001) -> float:
    """Magnitude of integral of g_t along z=R+i y, 0 <= y <= h."""
    y = np.linspace(0.0, h, n_y)
    z = R + 1j * y
    values_dz = np.exp(-0.5 * t * z**2) * 1j
    return float(abs(np.trapezoid(values_dz, y)))


def upper_hitting_probability(
    mu: np.ndarray, a: float = 1.0, b: float = 1.0, sigma: float = 1.0
) -> np.ndarray:
    r"""P_0(tau_b < tau_{-a}) for X_t = mu t + sigma W_t.

    For mu != 0,

        p_+(mu) = [1-exp(-2 mu a/sigma^2)]
                    / [1-exp(-2 mu(a+b)/sigma^2)].

    The continuous extension at mu=0 is a/(a+b).
    """
    mu = np.asarray(mu, dtype=float)
    k = 2.0 * mu / sigma**2
    numerator = -np.expm1(-k * a)
    denominator = -np.expm1(-k * (a + b))
    result = np.empty_like(mu)
    near_zero = np.abs(k) < 1e-8
    result[near_zero] = a / (a + b)
    np.divide(numerator, denominator, out=result, where=~near_zero)
    return result


def add_arrow(
    ax: plt.Axes,
    start: tuple[float, float],
    end: tuple[float, float],
    color: str,
    mutation_scale: float = 15,
) -> None:
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=mutation_scale,
            linewidth=2.0,
            color=color,
            shrinkA=0,
            shrinkB=0,
            zorder=8,
        )
    )


def main() -> None:
    mpl.rcParams.update(
        {
            "figure.facecolor": "#FAF8FC",
            "axes.facecolor": "#FFFFFF",
            "axes.edgecolor": "#B9A7C7",
            "axes.labelcolor": "#330662",
            "axes.titlecolor": "#330662",
            "axes.titlesize": 13,
            "axes.titleweight": "bold",
            "font.family": "DejaVu Sans",
            "font.size": 10.5,
            "grid.color": "#E6DCEC",
            "grid.alpha": 0.75,
            "legend.frameon": True,
            "legend.framealpha": 0.94,
            "legend.edgecolor": "#D9C9E4",
            "mathtext.fontset": "dejavusans",
        }
    )

    navy = "#330662"       # NYU deep violet
    blue = "#57068C"       # NYU Violet
    cyan = "#8900E1"       # bright violet accent
    gold = "#E6A11A"
    coral = "#C5283D"
    green = "#00856A"
    purple = "#7B2CBF"

    fig = plt.figure(figsize=(18.5, 11.7), layout="constrained")
    gs = fig.add_gridspec(2, 3, width_ratios=[1.0, 1.0, 1.12])

    # A. Sample paths and the exact sqrt(t) spatial scale.
    ax = fig.add_subplot(gs[0, 0])
    rng = np.random.default_rng(712)
    horizon, n_steps, n_paths = 2.0, 800, 36
    times = np.linspace(0.0, horizon, n_steps + 1)
    dt = times[1] - times[0]
    increments = np.sqrt(dt) * rng.standard_normal((n_paths, n_steps))
    paths = np.column_stack(
        [np.zeros(n_paths), np.cumsum(increments, axis=1)]
    )
    one_sd = np.sqrt(times)
    ax.fill_between(times, -2 * one_sd, 2 * one_sd, color=blue, alpha=0.08)
    ax.fill_between(times, -one_sd, one_sd, color=blue, alpha=0.13)
    for j in range(n_paths):
        ax.plot(times, paths[j], color=navy, alpha=0.13, linewidth=0.85)
    for j, color in zip([4, 12, 25], [blue, coral, gold]):
        ax.plot(times, paths[j], color=color, linewidth=1.65, alpha=0.95)
    ax.plot(times, one_sd, color=blue, linestyle="--", linewidth=1.1)
    ax.plot(times, -one_sd, color=blue, linestyle="--", linewidth=1.1)
    ax.text(1.95, np.sqrt(1.95) + 0.08, r"$+\sqrt{t}$", color=blue, ha="right")
    ax.text(1.95, -np.sqrt(1.95) - 0.30, r"$-\sqrt{t}$", color=blue, ha="right")
    ax.set(
        title=r"A  Brownian paths: spatial scale $\sqrt{t}$",
        xlabel=r"time $t$",
        ylabel=r"position $W_t$",
        xlim=(0, horizon),
    )
    ax.grid(True)

    # B. Heat kernel at several exact time marginals.
    ax = fig.add_subplot(gs[0, 1])
    x = np.linspace(-4.5, 4.5, 1000)
    times_to_plot = [0.20, 0.50, 1.00, 2.00]
    colors = [gold, coral, blue, purple]
    for time_value, color in zip(times_to_plot, colors):
        density = heat_kernel(x, time_value)
        ax.plot(x, density, color=color, linewidth=2.2, label=fr"$t={time_value:g}$")
    ax.set(
        title=r"B  Exact marginals $W_t\sim N(0,t)$",
        xlabel=r"position $x$",
        ylabel=r"$p_t(x)$",
        xlim=(-4.5, 4.5),
        ylim=(0, 0.94),
    )
    ax.text(
        0.02,
        0.95,
        r"$p_t(x)=\dfrac{1}{\sqrt{2\pi t}}e^{-x^2/(2t)}$",
        transform=ax.transAxes,
        va="top",
        color=navy,
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": "#D9C9E4"},
    )
    ax.legend(loc="upper right")
    ax.grid(True)

    # C. Finite-frequency Fourier inversion converging to the exact kernel.
    ax = fig.add_subplot(gs[0, 2])
    t0 = 1.0
    x_inv = np.linspace(-4.5, 4.5, 501)
    exact = heat_kernel(x_inv, t0)
    for cutoff, color in zip([1.5, 3.0, 6.0], [gold, cyan, purple]):
        approximate = truncated_fourier_inverse(x_inv, t0, cutoff)
        ax.plot(
            x_inv,
            approximate,
            color=color,
            linewidth=1.8,
            label=fr"cutoff $R={cutoff:g}$",
        )
    ax.plot(x_inv, exact, color=navy, linewidth=2.8, label="exact heat kernel")
    ax.set(
        title=r"C  Inverse Fourier integral as $R\to\infty$",
        xlabel=r"position $x$",
        ylabel=r"$\frac{1}{2\pi}\int_{-R}^{R}e^{-i\omega x-t\omega^2/2}\,d\omega$",
        xlim=(-4.5, 4.5),
        ylim=(-0.035, 0.44),
    )
    ax.axhline(0, color="#8090A5", linewidth=0.8)
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(True)

    # D. The rectangle used to shift the line of integration.
    ax = fig.add_subplot(gs[1, 0])
    x0, t0, R0 = 1.5, 1.0, 3.0
    h = x0 / t0
    u = np.linspace(-3.7, 3.7, 500)
    v = np.linspace(-0.35, 1.95, 360)
    U, V = np.meshgrid(u, v)
    log_abs = -0.5 * t0 * (U**2 - V**2)  # log |exp(-t z^2/2)|
    levels = np.linspace(-7.0, 2.0, 19)
    contour = ax.contourf(U, V, np.clip(log_abs, levels[0], levels[-1]), levels=levels, cmap="magma")
    contour_lines = ax.contour(U, V, log_abs, levels=[-4, -2, 0, 1], colors="white", linewidths=0.45, alpha=0.55)
    ax.clabel(contour_lines, inline=True, fontsize=7, fmt=r"$%g$")
    contour_color = "#5BE7FF"
    ax.plot([-R0, R0], [0, 0], color=contour_color, linewidth=2.6)
    ax.plot([R0, R0], [0, h], color=contour_color, linewidth=2.6)
    ax.plot([R0, -R0], [h, h], color=contour_color, linewidth=2.6)
    ax.plot([-R0, -R0], [h, 0], color=contour_color, linewidth=2.6)
    add_arrow(ax, (-0.7, 0), (0.7, 0), contour_color)
    add_arrow(ax, (R0, 0.45), (R0, 1.05), contour_color)
    add_arrow(ax, (0.7, h), (-0.7, h), contour_color)
    add_arrow(ax, (-R0, 1.05), (-R0, 0.45), contour_color)
    ax.text(0, -0.13, r"$C_0:\ \operatorname{Im}z=0$", color="white", ha="center", va="top")
    ax.text(0, h + 0.08, r"$C_h:\ \operatorname{Im}z=h=x/t$", color="white", ha="center", va="bottom")
    ax.text(R0 + 0.10, 0.75, r"$V_R$", color="white", va="center")
    ax.text(-R0 - 0.10, 0.75, r"$V_{-R}$", color="white", va="center", ha="right")
    ax.set(
        title=r"D  Rectangle for $g_t(z)=e^{-tz^2/2}$ (entire)",
        xlabel=r"$\operatorname{Re}z$",
        ylabel=r"$\operatorname{Im}z$",
        xlim=(-3.7, 3.7),
        ylim=(-0.35, 1.95),
    )
    colorbar = fig.colorbar(contour, ax=ax, shrink=0.88, pad=0.02)
    colorbar.set_label(r"$\log |g_t(z)|=-\frac{t}{2}(u^2-v^2)$")

    # E. Direct numerical vertical-edge integral and its rigorous upper bound.
    ax = fig.add_subplot(gs[1, 1])
    R_values = np.linspace(0.55, 5.25, 125)
    edge_values = np.array([vertical_edge_integral(R, t0, h) for R in R_values])
    upper_bound = abs(h) * np.exp(-0.5 * t0 * R_values**2 + 0.5 * t0 * h**2)
    ax.semilogy(R_values, edge_values, color=blue, linewidth=2.2, label=r"$|\int_{V_R}g_t(z)\,dz|$")
    ax.semilogy(
        R_values,
        upper_bound,
        color=coral,
        linestyle="--",
        linewidth=2.0,
        label=r"$|h|e^{-tR^2/2+th^2/2}$",
    )
    ax.set(
        title="E  Why both vertical edges vanish",
        xlabel=r"rectangle half-width $R$",
        ylabel="magnitude (log scale)",
        xlim=(R_values[0], R_values[-1]),
    )
    ax.legend(loc="upper right")
    ax.grid(True, which="both")
    ax.text(
        0.03,
        0.05,
        r"Fixed $h=x/t$; hence the Gaussian $e^{-tR^2/2}$ wins.",
        transform=ax.transAxes,
        color=navy,
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "edgecolor": "#D9C9E4"},
    )

    # F. Corrected boundary-hitting probabilities for the drifted process.
    ax = fig.add_subplot(gs[1, 2])
    mu = np.linspace(-3.0, 3.0, 801)
    p_up = upper_hitting_probability(mu, a=1.0, b=1.0, sigma=1.0)
    ax.plot(mu, p_up, color=green, linewidth=2.7, label=r"$P_0(\tau_{+1}<\tau_{-1})$")
    ax.plot(mu, 1.0 - p_up, color=coral, linewidth=2.7, label=r"$P_0(\tau_{-1}<\tau_{+1})$")
    ax.axvline(0, color="#7C8798", linestyle="--", linewidth=1.2)
    ax.axhline(0.5, color="#7C8798", linestyle=":", linewidth=1.2)
    ax.scatter([0], [0.5], color=navy, s=38, zorder=5)
    ax.set(
        title="F  Drifted Brownian motion: corrected directions",
        xlabel=r"drift $\mu$  (with $\sigma=a=b=1$)",
        ylabel="first-hit probability",
        xlim=(-3, 3),
        ylim=(-0.025, 1.025),
    )
    ax.text(
        0.03,
        0.95,
        r"$p_+(\mu)=\dfrac{1-e^{-2\mu a/\sigma^2}}{1-e^{-2\mu(a+b)/\sigma^2}}$",
        transform=ax.transAxes,
        va="top",
        color=navy,
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": "#D9C9E4"},
    )
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(True)

    fig.suptitle(
        "Brownian density from Fourier inversion and a zero-residue contour shift",
        fontsize=21,
        fontweight="bold",
        color=navy,
    )
    fig.text(
        0.5,
        -0.012,
        "Key point: no poles lie inside the contour, so the residue theorem gives 0; the Gaussian normalization is a separate input.",
        ha="center",
        va="top",
        color="#5B3A6E",
        fontsize=10.5,
    )
    temporary = OUTPUT.with_name(f".{OUTPUT.stem}.tmp.png")
    fig.savefig(temporary, dpi=260, bbox_inches="tight", facecolor=fig.get_facecolor())
    if not temporary.is_file() or temporary.stat().st_size == 0:
        raise RuntimeError(f"failed to render {OUTPUT}")
    temporary.replace(OUTPUT)
    print(f"Saved {OUTPUT}")


if __name__ == "__main__":
    main()
