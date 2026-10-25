"""Parser tests for run_regression.py (no scientific script is executed).

The regression runner reads one named field from each script's output.  In v1.8.3 the producers'
messages were translated into English and three field names changed ("differenza" -> "difference").
These cases check, with the runner's own patterns, that
  - the current English output is read, with the right value and sample count;
  - the old Italian field name is no longer read (a stale producer fails instead of passing silently);
  - malformed or non-finite tokens are refused or reported as non-finite;
  - prose lines ("if diff~1e-12") never become samples;
  - end to end, through run_one with a simulated script output: non-finite values, a wrong sample
    count, a missing field, a non-zero exit and a timeout give FAIL/NOVALUE/ERROR/TIMEOUT, never PASS.
Run from NonStationaryMetrics/:  python3 tests/test_run_regression_parser.py   (exit status 1 on failure)
"""
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import run_regression as R  # noqa: E402

PAT = {rel.split('/')[-1]: pat for rel, kind, pat, bound, label, nexp in R.CHECKS}

CASES = [
    # (script, output text, expected list of parsed values)
    ("pipeline_completa_deltaphi.py", "  difference  = 3.21e-14   <-- chain closed\n", [3.21e-14]),
    ("pipeline_completa_deltaphi.py", "  differenza   = 3.21e-14   <-- catena chiusa\n", []),
    ("pipeline_completa_deltaphi.py", "  difference  = 3.21e-14x\n", []),
    ("pipeline_completa_deltaphi.py", "  difference  = nan\n", [float('nan')]),
    ("pipeline_completa_deltaphi.py", "  difference  = inf\n", [float('inf')]),
    ("pipeline_completa_deltaphi.py", "  difference  = 1e999\n", [float('inf')]),
    ("kerr_psi_explicit_verified.py",
     "  r=11.0: dEF=1e-3  rhs=1e-3  diff=2.10e-16\n"
     "  r=9.0: dEF=1e-3  rhs=1e-3  diff=0.00e+00\n"
     "=> if diff~1e-12: the reduction is CORRECT (Delta handled properly).\n"
     "difference = 4.40e-13   (identity, coefficients GIVEN by the c_k, NO fit)\n",
     [2.10e-16, 0.0, 4.40e-13]),
    ("kerr_psi_explicit_verified.py", "differenza = 4.40e-13\n", []),
    ("vaidya_tau_assembly.py", "  difference = 1.10e-15   <-- Vaidya tau branch, all algebraic\n", [1.10e-15]),
    ("vaidya_tau_assembly.py", "  differenza = 1.10e-15\n", []),
    ("vaidya_tau_assembly.py", "  difference = 0e+\n", []),
]


def same(a, b):
    if len(a) != len(b):
        return False
    return all((x != x and y != y) or x == y for x, y in zip(a, b))


CHECK = {rel.split('/')[-1]: (rel, kind, pat, bound, label, nexp) for rel, kind, pat, bound, label, nexp in R.CHECKS}

# (script, simulated stdout, return code, raise timeout?, expected status)
RUN_CASES = [
    ("vaidya_tau_assembly.py", "  difference = 1.10e-15\n", 0, False, "PASS"),
    ("vaidya_tau_assembly.py", "  difference = 1.10e-3\n", 0, False, "FAIL"),
    ("vaidya_tau_assembly.py", "  difference = nan\n", 0, False, "FAIL"),
    ("vaidya_tau_assembly.py", "  difference = inf\n", 0, False, "FAIL"),
    ("vaidya_tau_assembly.py", "  difference = 1e999\n", 0, False, "FAIL"),
    ("vaidya_tau_assembly.py", "  differenza = 1.10e-15\n", 0, False, "NOVALUE"),
    ("vaidya_tau_assembly.py", "  difference = 1.10e-15\n  difference = 2.0e-15\n", 0, False, "FAIL"),
    ("vaidya_tau_assembly.py", "  difference = 1.10e-15\n", 1, False, "ERROR"),
    ("vaidya_tau_assembly.py", "", 0, True, "TIMEOUT"),
]


def run_case(script, stdout, code, timeout):
    rel, kind, pat, bound, label, nexp = CHECK[script]

    def fake_run(cmd, **kw):
        if timeout:
            raise subprocess.TimeoutExpired(cmd, kw.get("timeout"))
        return subprocess.CompletedProcess(cmd, code, stdout=stdout, stderr="simulated error")
    real = R.subprocess.run
    R.subprocess.run = fake_run
    try:
        return R.run_one(rel, kind, pat, bound, label, 10, nexp)
    finally:
        R.subprocess.run = real


def main() -> int:
    bad = 0
    for script, text, want in CASES:
        got, _ = R.all_numbers(text, PAT[script])
        ok = same(got, want)
        bad += not ok
        print(f"{'ok ' if ok else 'KO '} {script:32s} {text.strip().splitlines()[0][:50]!r:54s} -> {got}")
    print("-- end to end through run_one (simulated script output) --")
    for script, stdout, code, timeout, want in RUN_CASES:
        status, msg, _ = run_case(script, stdout, code, timeout)
        ok = status == want
        bad += not ok
        shown = stdout.strip().splitlines()[0][:40] if stdout.strip() else "(timeout)"
        print(f"{'ok ' if ok else 'KO '} {script:24s} {shown!r:44s} exit={code} -> {status:8s} {msg}")
    print("ALL CASES AS EXPECTED" if bad == 0 else f"{bad} CASE(S) WRONG")
    return 0 if bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
