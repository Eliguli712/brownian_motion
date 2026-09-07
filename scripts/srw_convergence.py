#!/usr/bin/env python3
"""NYU-violet visualization of SRW/ARW convergence to Brownian motion.

The four panels show:
1. diffusion-scale sample paths of a scaled symmetric random walk;
2. exact and empirical first/second moments;
3. locally uniform convergence of characteristic functions on a fixed compact
   frequency window;
4. convergence of discrete gambler's-ruin probabilities to the drifted
   Brownian exit probability.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


OUTPUT = (
    Path(__file__).resolve().parents[1]
    / "figures"
    / "srw_convergence_visualization.png"
)

NYU_VIOLET = "#57068C"
NYU_DEEP = "#330662"
NYU_BRIGHT = "#8900E1"
NYU_LILAC = "#B8A1C8"
GOLD = "#E6A11A"
GREEN = "#00856A"
CORAL = "#C5283D"


def scaled_walk_paths(
    rng: np.random.Generator,
    n_paths: int,
    dt: float,
    horizon: float,
    sigma: float = 1.0,
) -> tuple[np.ndarray, np.ndarray]:
    n_steps = int(round(horizon / dt))
    increments = sigma * np.sqrt(dt) * rng.choice(
        np.array([-1.0, 1.0]), size=(n_paths, n_steps)
    )
    paths = np.column_stack([np.zeros(n_paths), np.cumsum(increments, axis=1)])
    times = np.linspace(0.0, n_steps * dt, n_steps + 1)
    return times, paths


def biased_increment_cf(u: np.ndarray, dt: float, mu: float, sigma: float) -> np.ndarray:
    eps = (mu / sigma) * np.sqrt(dt)
    p = 0.5 * (1.0 + eps)
    q = 1.0 - p
    dx = sigma * np.sqrt(dt)
    return p * np.exp(1j * u * dx) + q * np.exp(-1j * u * dx)


def brownian_cf(u: np.ndarray, t: float, mu: float, sigma: float) -> np.ndarray:
    return np.exp(1j * mu * t * u - 0.5 * sigma**2 * t * u**2)


def brownian_upper_hit(mu: float, sigma: float, lower: float, upper: float) -> float:
    if abs(mu) < 1e-14:
        return lower / (lower + upper)
    num = -np.expm1(-2.0 * mu * lower / sigma**2)
    den = -np.expm1(-2.0 * mu * (lower + upper) / sigma**2)
    return float(num / den)


def discrete_upper_hit(
    dt: float, mu: float, sigma: float, lower: float, upper: float
) -> float:
    dx = sigma * np.sqrt(dt)
    a = int(round(lower / dx))
    b = int(round(upper / dx))
    eps = (mu / sigma) * np.sqrt(dt)
    p = 0.5 * (1.0 + eps)
    q = 1.0 - p
    if abs(p - q) < 1e-14:
        return a / (a + b)
    log_r = np.log(q / p)
    numerator = -np.expm1(a * log_r)
    denominator = -np.expm1((a + b) * log_r)
    return float(numerator / denominator)


def main() -> None:
    mpl.rcParams.update(
        {
            "figure.facecolor": "#FAF8FC",
            "axes.facecolor": "#FFFFFF",
            "axes.edgecolor": NYU_LILAC,
            "axes.labelcolor": NYU_DEEP,
            "axes.titlecolor": NYU_DEEP,
            "axes.titlesize": 13,
            "axes.titleweight": "bold",
            "font.family": "DejaVu Sans",
            "font.size": 10.5,
            "grid.color": "#E6DCEC",
            "grid.alpha": 0.82,
            "legend.framealpha": 0.95,
            "legend.edgecolor": "#D9C9E4",
            "mathtext.fontset": "dejavusans",
        }
    )

    fig, axes = plt.subplots(2, 2, figsize=(15.8, 10.4), layout="constrained")
    rng = np.random.default_rng(712)

    # A. Scaled SRW paths and diffusion envelope.
    ax = axes[0, 0]
    dt, horizon, sigma = 0.01, 2.0, 1.0
    times, paths = scaled_walk_paths(rng, 34, dt, horizon, sigma)
    envelope = sigma * np.sqrt(times)
    ax.fill_between(times, -2 * envelope, 2 * envelope, color=NYU_VIOLET, alpha=0.08)
    ax.fill_between(times, -envelope, envelope, color=NYU_VIOLET, alpha=0.14)
    for path in paths:
        ax.plot(times, path, color=NYU_DEEP, alpha=0.14, linewidth=0.9)
    for j, color in zip([2, 10, 19], [NYU_VIOLET, GOLD, CORAL]):
        ax.plot(times, paths[j], color=color, linewidth=1.8)
    ax.plot(times, envelope, color=NYU_BRIGHT, linestyle="--", linewidth=1.2)
    ax.plot(times, -envelope, color=NYU_BRIGHT, linestyle="--", linewidth=1.2)
    ax.set(
        title=r"A  Scaled SRW: $\Delta X=\pm\sqrt{\Delta t}$",
        xlabel=r"time $t$",
        ylabel=r"$S_t^{(\Delta t)}$",
        xlim=(0, horizon),
    )
    ax.grid(True)

    # B. Exact moments and Monte Carlo estimates for the biased walk.
    ax = axes[0, 1]
    mu, sigma, dt, horizon = 0.65, 1.15, 0.005, 2.0
    n_paths = 30_000
    n_steps = int(round(horizon / dt))
    eps = (mu / sigma) * np.sqrt(dt)
    p = 0.5 * (1.0 + eps)
    bernoulli = rng.random((n_paths, n_steps)) < p
    increments = sigma * np.sqrt(dt) * np.where(bernoulli, 1.0, -1.0)
    sample = np.cumsum(increments, axis=1)
    sample_times = np.arange(1, n_steps + 1) * dt
    stride = 8
    tt = sample_times[stride - 1 :: stride]
    empirical_mean = sample[:, stride - 1 :: stride].mean(axis=0)
    empirical_var = sample[:, stride - 1 :: stride].var(axis=0)
    exact_mean = mu * tt
    exact_var = sigma**2 * tt - mu**2 * tt * dt
    ax.plot(tt, empirical_mean, color=NYU_VIOLET, linewidth=2.2, label="empirical mean")
    ax.plot(tt, exact_mean, color=NYU_DEEP, linestyle="--", linewidth=1.8, label=r"$\mu t$")
    ax.plot(tt, empirical_var, color=GOLD, linewidth=2.2, label="empirical variance")
    ax.plot(
        tt,
        exact_var,
        color=CORAL,
        linestyle="--",
        linewidth=1.8,
        label=r"$\sigma^2t-\mu^2t\Delta t$",
    )
    ax.set(
        title="B  Drift and variance matching",
        xlabel=r"time $t$",
        ylabel="moment value",
        xlim=(0, horizon),
    )
    ax.legend(loc="upper left", ncols=2, fontsize=8.8)
    ax.grid(True)

    # C. Characteristic functions converge uniformly on each fixed compact set.
    ax = axes[1, 0]
    mu, sigma, horizon = 0.7, 1.1, 1.0
    u = np.linspace(-5.0, 5.0, 1301)
    target = brownian_cf(u, horizon, mu, sigma)
    for dt_value, color in zip([0.20, 0.05, 0.0125], [GOLD, NYU_BRIGHT, NYU_VIOLET]):
        n_steps = int(round(horizon / dt_value))
        discrete = biased_increment_cf(u, dt_value, mu, sigma) ** n_steps
        error = np.abs(discrete - target)
        ax.semilogy(u, np.maximum(error, 1e-12), color=color, linewidth=2.0, label=fr"$\Delta t={dt_value:g}$")
    ax.set(
        title=r"C  CF error on compact $u\in[-5,5]$",
        xlabel=r"frequency $u$",
        ylabel=r"$|\phi_{\Delta t}(u)-\phi(u)|$",
        xlim=(-5, 5),
        ylim=(1e-7, 1.2),
    )
    ax.legend(loc="upper right")
    ax.grid(True, which="both")

    # D. Discrete gambler's ruin converges to the Brownian exit formula.
    ax = axes[1, 1]
    mu, sigma, lower, upper = 0.55, 1.0, 1.0, 1.0
    m_values = np.array([5, 8, 12, 18, 25, 40, 60, 90, 130, 180])
    # Choose dx=1/m so boundaries are exact lattice points; then dt=(dx/sigma)^2.
    dt_values = (1.0 / (sigma * m_values)) ** 2
    probs = np.array(
        [discrete_upper_hit(dt_value, mu, sigma, lower, upper) for dt_value in dt_values]
    )
    limit = brownian_upper_hit(mu, sigma, lower, upper)
    ax.semilogx(dt_values, probs, marker="o", color=NYU_VIOLET, linewidth=2.3, label="discrete ARW")
    ax.axhline(limit, color=GREEN, linestyle="--", linewidth=2.0, label="Brownian limit")
    ax.invert_xaxis()
    ax.set(
        title="D  Gambler's ruin: discrete-to-continuous limit",
        xlabel=r"time step $\Delta t$  (decreases to the right)",
        ylabel=r"$P_0(\tau_{+1}<\tau_{-1})$",
    )
    ax.legend(loc="best")
    ax.grid(True, which="both")

    fig.suptitle(
        "From symmetric random walks to Brownian motion",
        fontsize=21,
        fontweight="bold",
        color=NYU_DEEP,
    )
    fig.text(
        0.5,
        -0.012,
        "NYU Courant-style map: moments, characteristic functions, and optional stopping survive the diffusion scaling.",
        ha="center",
        va="top",
        color="#5B3A6E",
    )
    temporary = OUTPUT.with_name(f".{OUTPUT.stem}.tmp.png")
    fig.savefig(temporary, dpi=270, bbox_inches="tight", facecolor=fig.get_facecolor())
    if not temporary.is_file() or temporary.stat().st_size == 0:
        raise RuntimeError(f"failed to render {OUTPUT}")
    temporary.replace(OUTPUT)
    print(f"Saved {OUTPUT}")


if __name__ == "__main__":
    main()
