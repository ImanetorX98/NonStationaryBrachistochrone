"""The JSON header of the 16 producers: identical everywhere, and harmless outside the repository.

Checks that
  - every producer carries the same header block (only the source expression differs: __file__ for
    Python, sys.argv[0] for .sage files run by the Sage CLI);
  - a producer copied out of the repository, without verification_results.py, runs as in v1.8.3 and
    exits 0 when NSB_RESULT_JSON is unset, and stops with a clear message when it is set;
  - inside the repository the same producer writes a valid JSON result.
The run uses VaidyaMetric/vaidya_tau_assembly.py (about a second; needs NumPy and SymPy).
Run from NonStationaryMetrics/:  python3 tests/test_producer_header.py   (exit status 1 on failure)
"""
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

NSM = Path(__file__).resolve().parents[1]
ROOT = NSM.parent
sys.path.insert(0, str(NSM))
import run_regression as R  # noqa: E402
from verification_results import read_results  # noqa: E402

HEADER = re.compile(r"# Stable verification fields are independent of the human-readable messages\.\n.*?"
                    r"\n    _results = _NoResults\(\)\n", re.S)
SAGE_PRODUCERS = ["KerrScripts/kerr_quasiperiods_bel.sage", "KerrScripts/kerr_thirdkind_theta_closed.sage",
                  "KerrScripts/kerr_tbranch_genus2_periods.sage.py",
                  "NonStationaryMetrics/KerrSessionScripts/kerr_tau_general_genus2_periods.sage.py"]


def main():
    bad = 0
    def check(label, ok):
        nonlocal bad
        bad += not ok
        print(f"{'ok ' if ok else 'KO '} {label}")

    producers = [rel for rel, *_ in R.CHECKS] + SAGE_PRODUCERS
    blocks = {}
    for rel in producers:
        m = HEADER.findall((ROOT / rel).read_text(encoding="utf-8"))
        check(f"one header block in {rel}", len(m) == 1)
        if m:
            text = m[0].replace("_result_sys.argv[0]", "SOURCE").replace("__file__", "SOURCE")
            blocks[rel] = text.replace("  # Sage CLI keeps the source in argv[0].", "")
    check(f"the {len(blocks)} header blocks are identical", len(set(blocks.values())) == 1)

    rel = "NonStationaryMetrics/VaidyaMetric/vaidya_tau_assembly.py"
    with tempfile.TemporaryDirectory() as folder:
        alone = Path(folder) / "vaidya_tau_assembly.py"
        shutil.copy(ROOT / rel, alone)
        env = {k: v for k, v in os.environ.items() if k != "NSB_RESULT_JSON"}
        p = subprocess.run([sys.executable, alone.name], cwd=folder, capture_output=True, text=True, env=env)
        check("outside the repository, no JSON requested: exit 0 and the usual output",
              p.returncode == 0 and "difference" in p.stdout)
        env["NSB_RESULT_JSON"] = str(Path(folder) / "result.json")
        p = subprocess.run([sys.executable, alone.name], cwd=folder, capture_output=True, text=True, env=env)
        check("outside the repository, JSON requested: non-zero exit with a clear message, no file",
              p.returncode != 0 and "verification_results.py not found" in p.stderr
              and not (Path(folder) / "result.json").exists())
        result = Path(folder) / "inside.json"
        env["NSB_RESULT_JSON"] = str(result)
        p = subprocess.run([sys.executable, "vaidya_tau_assembly.py"], cwd=ROOT / Path(rel).parent,
                           capture_output=True, text=True, env=env)
        try:
            read_results(result, rel, {"residual": R.RESULT_SAMPLES["vaidya_tau_assembly.py"]})
            ok = p.returncode == 0
        except (OSError, ValueError):
            ok = False
        check("inside the repository, JSON requested: valid result", ok)
    print("ALL CASES AS EXPECTED" if bad == 0 else f"{bad} CASE(S) WRONG")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
