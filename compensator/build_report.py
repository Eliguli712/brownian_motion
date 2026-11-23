"""Regenerate numerical macros and compile the English paper with embedded fonts."""
from pathlib import Path
import json
import shutil
import subprocess

root = Path(__file__).resolve().parent
m = json.loads((root / "data/metrics.json").read_text())

def sci(value):
    a, b = f"{value:.3e}".split("e")
    return "$"+a+r"\times10^{"+str(int(b))+"}$"

vals = {
    "SamplePaths": str(m["parameters"]["paths"]),
    "FineSteps": str(m["parameters"]["fine_steps"]),
    "RandomSeed": str(m["parameters"]["seed"]),
    "SemigroupError": sci(m["max_semigroup_vs_quad_error"]),
    "PrimitiveError": sci(m["max_primitive_vs_quad_error"]),
    "ExpectedP": f"${m['deterministic_E_P_T']:.12f}$",
    "ExpectedK": f"${m['deterministic_E_K_T']:.12f}$",
    "MeanM": sci(m["MC_mean_M_T"]),
    "SEMeanM": sci(m["MC_se_M_T"]),
    "CoarseError": f"${m['convergence'][0][2]:.7f}$",
    "FineError": f"${m['fine_rms_ito_residual']:.7f}$",
    "ConvergenceSlope": f"${m['convergence_slope']:.4f}$",
    "PDEError": sci(m["PDE_finest_rms"]),
    "EnergyAdaptive": f"${m['adaptive_spatial_energy']:.10f}$",
    "EnergyBound": f"${m['energy_upper_bound']:.8f}$",
    "ZMean": f"{m['MC_mean_M_T']/m['MC_se_M_T']:.3f}",
}
(root / "data/results_macros.tex").write_text(
    "\n".join("\\newcommand{\\"+k+"}{"+v+"}" for k, v in vals.items())+"\n"
)
(root / "logs").mkdir(exist_ok=True)
if shutil.which("xelatex"):
    fmt = subprocess.run(["kpsewhich", "xelatex.fmt"], text=True, capture_output=True).stdout.strip()
    cmd = ["xelatex"]
    if not fmt:
        build = root / ".build" / "texfmt"
        build.mkdir(parents=True, exist_ok=True)
        fmtfile = build / "xelatex.fmt"
        if not fmtfile.exists():
            with open(root / "logs" / "format_build.log", "w") as out:
                subprocess.run(["xetex", "-ini", "-etex", "-jobname=xelatex",
                                "-progname=xelatex", "xelatex.ini"],
                               cwd=build, stdout=out, stderr=subprocess.STDOUT, check=True)
        cmd.append("-fmt="+str(fmtfile))
    cmd += ["-interaction=nonstopmode", "-halt-on-error", "paper.tex"]
    for i in [1, 2]:
        with open(root / "logs" / f"latex_pass{i}.log", "w") as out:
            subprocess.run(cmd, cwd=root, stdout=out, stderr=subprocess.STDOUT, check=True)
elif shutil.which("tectonic"):
    with open(root / "logs" / "tectonic.log", "w") as out:
        subprocess.run(["tectonic", "paper.tex", "--keep-logs"], cwd=root,
                       stdout=out, stderr=subprocess.STDOUT, check=True)
else:
    raise RuntimeError("A working XeLaTeX or Tectonic installation is required.")
print(root / "paper.pdf")
