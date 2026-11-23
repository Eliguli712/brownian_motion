"""Draw a vector bubble map of the actual proof dependencies and Lean coverage."""
from pathlib import Path
from xml.sax.saxutils import escape
import shutil
import subprocess

ROOT = Path(__file__).resolve().parent
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)
W, H = 1200, 1390
GREEN, AMBER, INK = "#19766F", "#A17E34", "#173756"
items = [
    '<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="1390" viewBox="0 0 1200 1390" role="img" aria-labelledby="title desc">',
    '<title id="title">Proof map: a Brownian functional and its predictable compensator</title>',
    '<desc id="desc">The central spine goes from the original nonsmooth payoff through Gaussian smoothing and the spatial primitive to the Ito decomposition and an L2 martingale. Side branches give kernel derivatives, growth and energy bounds, and the terminal extension. Solid teal bubbles identify compiled component statements under their stated hypotheses. Dashed amber bubbles identify analytic steps whose full Lean formalization remains open.</desc>',
    '<defs><marker id="arrow" markerWidth="10" markerHeight="10" refX="8" refY="4" orient="auto" markerUnits="strokeWidth"><path d="M0,0 L8,4 L0,8 Z" fill="#64798B"/></marker></defs>',
    '<rect width="1200" height="1390" fill="white"/>',
]


def txt(x, y, text, size=20, weight="normal", color=INK, anchor="middle", math=False):
    font = "Latin Modern Math, DejaVu Serif" if math else "DejaVu Sans, sans-serif"
    items.append(f'<text x="{x}" y="{y}" text-anchor="{anchor}" font-family="{font}" font-size="{size}" font-weight="{weight}" fill="{color}">{escape(text)}</text>')


def edge(d, pending=False):
    dash = ' stroke-dasharray="7 6"' if pending else ""
    items.append(f'<path d="{d}" fill="none" stroke="#64798B" stroke-width="2.3"{dash} marker-end="url(#arrow)"/>')


def bubble(cx, cy, rx, ry, title, lines, checked=True, status=None):
    color = GREEN if checked else AMBER
    fill = "#EAF5F4" if checked else "#FFF7E8"
    dash = '' if checked else ' stroke-dasharray="10 7"'
    items.append(f'<ellipse cx="{cx}" cy="{cy}" rx="{rx}" ry="{ry}" fill="{fill}" stroke="{color}" stroke-width="2.4"{dash}/>')
    txt(cx, cy-ry+27, status or ("LEAN-CHECKED COMPONENT" if checked else "ANALYTIC STEP / LEAN OPEN"), 12, "bold", color)
    first = cy-15 if len(lines) < 3 else cy-27
    txt(cx, first, title, 20, "bold")
    for i, (line, ismath) in enumerate(lines):
        txt(cx, first+28*(i+1), line, 21 if ismath else 14, color=INK if ismath else "#536475", math=ismath)


txt(600, 38, "A Brownian functional and its compensator", 30, "bold")
txt(600, 72, "The analytic proof, its dependencies, and the precise Lean verification boundary", 17, color="#536475")
items.append('<ellipse cx="224" cy="105" rx="18" ry="9" fill="#EAF5F4" stroke="#19766F" stroke-width="2"/>')
txt(250, 111, "Compiled component statement", 15, anchor="start")
items.append('<ellipse cx="684" cy="105" rx="18" ry="9" fill="#FFF7E8" stroke="#A17E34" stroke-width="2" stroke-dasharray="4 3"/>')
txt(710, 111, "Analytic link awaiting formalization", 15, anchor="start")

# Draw all connectors before the bubbles.
edge("M600 272 L600 318", True)
edge("M790 200 L840 200")
edge("M350 400 L395 400", True)
edge("M1000 265 L1000 375", True)
edge("M600 482 L600 565", True)
edge("M1000 555 L1000 650")
edge("M600 735 L600 820", True)
edge("M600 980 L600 1065", True)
edge("M1000 810 C1000 1060,1045 1160,842 1160", True)
edge("M355 940 C390 940,358 1040,375 1128", True)
edge("M410 680 C240 680,190 710,190 850", True)

bubble(600, 200, 190, 72, "Original payoff g", [
    ("finite max / indicator sum", False), ("measurability + growth", False)
])
bubble(1000, 200, 160, 65, "Payoff envelope", [
    ("|g(x)| ≤ A + |x|/4", True)
])
bubble(190, 400, 160, 70, "Heat-kernel PDE", [
    ("∂ᵥp = ½ ∂ₓₓp", True), ("actual kernel derivatives", False)
])
bubble(600, 400, 205, 82, "(a) Gaussian smoothing", [
    ("S(t,x) = E[g(x + σ Z)]", True), ("σ² = T − t", True), ("differentiate the original convolution", False)
], checked=False)
bubble(1000, 465, 160, 90, "Pointwise energy", [
    ("H² ≤ 2C²(W² + W⁴)", True), ("under |S(t,x)| ≤ C(1+|x|)", False)
])
bubble(600, 650, 210, 85, "(b) Primitive + generator", [
    ("F(t,x) = ∫₀ˣ r S(t,r) dr", True), ("Fₜ + ½ Fₓₓ = S(t,x) − ½ S(t,0)", True),
    ("FTC / IBP + algebra; given derivatives", False)
])
bubble(1000, 730, 160, 80, "Time majorant", [
    ("∫₀ᵀ 2C²(s + 3s²) ds", True), ("= C²(T² + 2T³)", True)
])
bubble(600, 900, 205, 80, "(c) Itô decomposition", [
    ("Pₜ = Kₜ + ∫₀ᵗ Hₛ dWₛ", True), ("construct the integral + apply Itô", False)
], checked=False)
bubble(190, 940, 165, 90, "Terminal extension", [
    ("Pₜ → ∫₀ᵂᵀ r g(r) dr", True), ("L² limits + finite variation", False)
], checked=False)
bubble(600, 1160, 240, 95, "(d) Compensated martingale", [
    ("Mₜ = Pₜ − Kₜ", True), ("E[Mₜ | ℱₛ] = Mₛ    (s ≤ t ≤ T)", True),
    ("square-integrable on the closed interval [0,T]", False)
], checked=False, status="ANALYTIC THEOREM / LEAN OPEN")

txt(600, 1302, "Kₜ = ∫₀ᵗ [S(s,Wₛ) − ½ S(s,0)] ds", 27, "normal", math=True)
txt(600, 1340, "Solid bubbles verify their listed component claims under explicit hypotheses.", 16, color="#536475")
txt(600, 1365, "The complete P − K martingale theorem remains an open Lean formalization obligation.", 16, color="#536475")
items.append('</svg>')
svg = FIG / "proof_map.svg"
svg.write_text("\n".join(items)+"\n")
inkscape = shutil.which("inkscape")
if inkscape:
    subprocess.run([inkscape, str(svg), "--export-type=pdf", "--export-text-to-path",
                    "--export-filename="+str(FIG/"proof_map.pdf")], check=True)
    subprocess.run([inkscape, str(svg), "--export-type=png", "--export-width=1500",
                    "--export-filename="+str(FIG/"proof_map.png")], check=True)
print(svg)
