#!/usr/bin/env python3
"""Generate every figure used by the Brownian-motion notes.

Run from any directory with::

    python /path/to/brownian_motion/scripts/generate_all.py

The random simulations use fixed seeds, so repeated runs are reproducible.
"""

from __future__ import annotations

from pathlib import Path
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[1]
FIGURE_DIR = PROJECT_ROOT / "figures"
sys.path.insert(0, str(Path(__file__).resolve().parent))

import brownian_contour  # noqa: E402
import levy_examples  # noqa: E402
import srw_convergence  # noqa: E402
import transform_tightness  # noqa: E402


def main() -> None:
    """Regenerate the five canonical PNG figures."""
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    srw_convergence.OUTPUT = FIGURE_DIR / "srw_convergence_visualization.png"
    brownian_contour.OUTPUT = FIGURE_DIR / "brownian_contour_visualization.png"
    transform_tightness.OUTPUT = FIGURE_DIR / "transform_tightness_visualization.png"

    srw_convergence.main()
    brownian_contour.main()
    transform_tightness.main()
    levy_examples.fig_escape(path=FIGURE_DIR / "levy_escape_panels.png")
    levy_examples.fig_pisot(path=FIGURE_DIR / "pisot_panels.png")

    expected = (
        "srw_convergence_visualization.png",
        "brownian_contour_visualization.png",
        "transform_tightness_visualization.png",
        "levy_escape_panels.png",
        "pisot_panels.png",
    )
    for name in expected:
        path = FIGURE_DIR / name
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"figure generation failed: {path}")
        print(f"verified {path.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()

