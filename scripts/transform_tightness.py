#!/usr/bin/env python3
"""NYU-violet visualization of transforms, Cantor CDFs, and tightness.

Panels:
1. Rademacher characteristic and moment generating functions;
2. Cantor staircases compressed toward a point mass at zero;
3. a fixed Cantor component plus mass escaping to infinity;
4. the characteristic-function boundary layer produced by that escape.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


OUTPUT = (
    Path(__file__).resolve().parents[1]
    / "figures"
    / "transform_tightness_visualization.png"
)

NYU_VIOLET = "#57068C"
NYU_DEEP = "#330662"
NYU_BRIGHT = "#8900E1"
NYU_LILAC = "#B8A1C8"
GOLD = "#E6A11A"
GREEN = "#00856A"
CORAL = "#C5283D"
PAPER = "#FAF8FC"


def cantor_cdf(x: np.ndarray, depth: int = 24) -> np.ndarray:
    """Evaluate the standard Cantor distribution function recursively."""
    values = np.asarray(x, dtype=float)
    result = np.zeros_like(values)
    result[values >= 1.0] = 1.0

    active = (values > 0.0) & (values < 1.0)
    y = values.copy()
    scale = np.ones_like(values)

    for _ in range(depth):
        if not np.any(active):
            break
        left = active & (y < 1.0 / 3.0)
        middle = active & (y >= 1.0 / 3.0) & (y <= 2.0 / 3.0)
        right = active & (y > 2.0 / 3.0)

        result[middle] += 0.5 * scale[middle]
        active[middle] = False

        y[left] *= 3.0
        scale[left] *= 0.5

        result[right] += 0.5 * scale[right]
        y[right] = 3.0 * y[right] - 2.0
        scale[right] *= 0.5

    result[active] += scale[active] * (y[active] >= 0.5)
    return np.clip(result, 0.0, 1.0)


def cantor_cf(t: np.ndarray, factors: int = 36) -> np.ndarray:
    """Characteristic function of C=sum 2 B_k/3^k."""
    t = np.asarray(t, dtype=float)
    product = np.ones_like(t)
    for k in range(1, factors + 1):
        product *= np.cos(t / (3.0**k))
    return np.exp(0.5j * t) * product


def uniform_shift_cf(t: np.ndarray, n: int) -> np.ndarray:
    """Characteristic function of Uniform[n,2n], stable at t=0."""
    t = np.asarray(t, dtype=float)
    out = np.ones_like(t, dtype=complex)
    nonzero = np.abs(t) > 1e-14
    out[nonzero] = (
        np.exp(2j * n * t[nonzero]) - np.exp(1j * n * t[nonzero])
    ) / (1j * n * t[nonzero])
    return out


def uniform_shift_cdf(x: np.ndarray, n: int) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    return np.clip((x - n) / n, 0.0, 1.0)


def main() -> None:
    mpl.rcParams.update(
        {
            "figure.facecolor": PAPER,
            "axes.facecolor": "#FFFFFF",
            "axes.edgecolor": NYU_LILAC,
            "axes.labelcolor": NYU_DEEP,
            "axes.titlecolor": NYU_DEEP,
            "axes.titlesize": 12.5,
            "axes.titleweight": "bold",
            "font.family": "DejaVu Sans",
            "font.size": 10.3,
            "grid.color": "#E6DCEC",
            "grid.alpha": 0.82,
            "legend.framealpha": 0.96,
            "legend.edgecolor": "#D9C9E4",
            "mathtext.fontset": "dejavusans",
        }
    )

    fig, axes = plt.subplots(2, 2, figsize=(15.8, 10.4), layout="constrained")

    # A. Rademacher CF and MGF.
    ax = axes[0, 0]
    z = np.linspace(-3.2, 3.2, 1001)
    ax.plot(z, np.cos(z), color=NYU_VIOLET, linewidth=2.5, label=r"CF: $\cos u$")
    ax.plot(z, np.cosh(z), color=GOLD, linewidth=2.5, label=r"MGF: $\cosh k$")
    ax.axhline(0.0, color=NYU_LILAC, linewidth=0.9)
    ax.set(
        title="A  Rademacher transforms: bounded oscillation vs growth",
        xlabel="real transform argument",
        ylabel="transform value",
        xlim=(-3.2, 3.2),
        ylim=(-1.35, 13.0),
    )
    ax.legend(loc="upper center", ncols=2)
    ax.grid(True)

    # B. Compressed Cantor staircases converge weakly to delta_0.
    ax = axes[0, 1]
    x = np.linspace(-0.035, 0.54, 2400)
    for n, color in zip([4, 8, 16], [GOLD, NYU_BRIGHT, NYU_VIOLET]):
        fn = cantor_cdf(n * x - 1.0)
        ax.plot(x, fn, color=color, linewidth=2.0, label=fr"$n={n}$")
    atom_cdf = (x >= 0.0).astype(float)
    ax.plot(x, atom_cdf, color=GREEN, linestyle="--", linewidth=2.2, label=r"$F_{\delta_0}$")
    ax.scatter([0.0], [1.0], color=GREEN, s=34, zorder=5)
    ax.scatter([0.0], [0.0], facecolors="white", edgecolors=NYU_VIOLET, s=36, zorder=5)
    ax.set(
        title="B  A devil's staircase compressed to the origin",
        xlabel=r"state variable $x$",
        ylabel=r"$F_n(x)$",
        xlim=(-0.035, 0.54),
        ylim=(-0.04, 1.04),
    )
    ax.legend(loc="lower right", ncols=2)
    ax.grid(True)

    # C. A Cantor part stays while a uniform part escapes.
    ax = axes[1, 0]
    a = 0.30
    x = np.linspace(-0.25, 15.0, 3200)
    defective = (1.0 - a) * cantor_cdf(x)
    for n, color in zip([2, 4, 7], [GOLD, NYU_BRIGHT, NYU_VIOLET]):
        fn = defective + a * uniform_shift_cdf(x, n)
        ax.plot(x, fn, color=color, linewidth=2.0, label=fr"$n={n}$")
    ax.plot(
        x,
        defective,
        color=CORAL,
        linestyle="--",
        linewidth=2.2,
        label=r"pointwise limit $(1-a)K$",
    )
    ax.axhline(1.0 - a, color=CORAL, linestyle=":", linewidth=1.2)
    ax.annotate(
        r"missing mass $a$",
        xy=(13.2, 1.0 - a),
        xytext=(11.0, 0.47),
        arrowprops={"arrowstyle": "->", "color": CORAL},
        color=CORAL,
        fontsize=10.5,
    )
    ax.set(
        title=r"C  Partial escape: $(1-a)\mu_C+a\,\mathrm{Unif}[n,2n]$",
        xlabel=r"state variable $x$",
        ylabel=r"$F_n(x)$",
        xlim=(-0.25, 15.0),
        ylim=(-0.04, 1.04),
    )
    ax.legend(loc="lower right")
    ax.grid(True)

    # D. The escape creates a shrinking boundary layer at frequency zero.
    ax = axes[1, 1]
    t = np.linspace(-0.65, 0.65, 2601)
    phi_c = cantor_cf(t)
    pointwise_nonzero = (1.0 - a) * phi_c.real
    for n, color in zip([5, 20, 80], [GOLD, NYU_BRIGHT, NYU_VIOLET]):
        phi_n = (1.0 - a) * phi_c + a * uniform_shift_cf(t, n)
        ax.plot(t, phi_n.real, color=color, linewidth=2.0, label=fr"$n={n}$")
    ax.plot(
        t,
        pointwise_nonzero,
        color=CORAL,
        linestyle="--",
        linewidth=2.2,
        label=r"limit for $t\ne0$",
    )
    ax.scatter([0.0], [1.0], color=GREEN, s=42, zorder=6, label=r"$\phi(0)=1$")
    ax.annotate(
        r"jump $a$",
        xy=(0.0, 1.0),
        xytext=(0.12, 0.86),
        arrowprops={"arrowstyle": "->", "color": GREEN},
        color=GREEN,
        fontsize=10.5,
    )
    ax.set(
        title="D  CF boundary layer: discontinuity records escaped mass",
        xlabel=r"frequency $t$",
        ylabel=r"$\operatorname{Re}\phi_n(t)$",
        xlim=(-0.65, 0.65),
        ylim=(0.60, 1.035),
    )
    ax.legend(loc="lower center", ncols=3, fontsize=8.7)
    ax.grid(True)

    fig.suptitle(
        "Transforms, devil's staircases, and tightness",
        fontsize=21,
        fontweight="bold",
        color=NYU_DEEP,
    )
    fig.text(
        0.5,
        -0.012,
        "NYU Courant-style map: bounded Fourier kernels preserve laws; discontinuity at frequency zero detects lost mass.",
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
