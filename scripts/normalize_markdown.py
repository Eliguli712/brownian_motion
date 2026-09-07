#!/usr/bin/env python3
"""Normalize Pandoc output for GitHub rendering.

This keeps the LaTeX source canonical while ensuring that generated image paths
resolve from notes/ and that the Markdown begins with a native heading.
"""

from __future__ import annotations

from pathlib import Path
import re


PROJECT_ROOT = Path(__file__).resolve().parents[1]
MARKDOWN = PROJECT_ROOT / "notes" / "brownian_motion_notes.md"


def main() -> None:
    text = MARKDOWN.read_text(encoding="utf-8")
    text = re.sub(r'<img src="([^"/]+\.png)"', r'<img src="../figures/\1"', text)
    title = "# From Random Walks to Brownian Motion\n\n"
    if not text.startswith(title):
        text = title + text
    text = text.replace("<div class=\"titlepage\">", "<div class=\"titlepage\" markdown=\"1\">")
    text = text.replace("<div class=\"keybox\">", "<div class=\"keybox\" markdown=\"1\">")
    MARKDOWN.write_text(text, encoding="utf-8")
    print(f"normalized {MARKDOWN.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()

