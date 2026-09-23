#!/usr/bin/env python3
"""Regression runner for the paper's closed-form residuals.

Runs each script listed in the manuscript's "Result-to-script map" table,
extracts its key numerical figure of merit (a residual or a convergence
slope), and checks it against the bound claimed in the paper. Exit code 0
iff every check passes.

Usage:
    python3 run_regression.py            # run all
    python3 run_regression.py --quick    # skip the slow symbolic/ODE scripts
    python3 run_regression.py -k vaidya  # only scripts whose key matches

The figure of merit is pulled from stdout with a per-script regex carrying one
capture group around the field being checked, and the captured token must be
complete -- see the notes on NUM and CHECKS.  A residual must satisfy
|x| <= tol; a slope must lie in [1.7, 2.3] (second-order canonical perturbation
theory).  Each check also declares how many samples it expects, so a
configuration that stops being printed fails loudly instead of shrinking the
sample set in silence.
"""
import argparse
import re
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent  # repo root (Brachistocrona/)

SCI = r"[-+]?\d+\.?\d*[eE][-+]?\d+"

# Any numeric field: integer, decimal, scientific, or a non-finite word.  NaN and
# Inf are matched DELIBERATELY so that the finiteness test downstream can see
# them; a pattern that only matched digits would drop them silently and the test
# would be unreachable code.
# The trailing lookahead forbids a truncated or malformed token: without it the
# pattern matches a PREFIX, so "diff=0junk" and "diff=0e+" both read as 0.0 and
# "diff=1.5e-15x" silently drops the x.  Excluding a following LETTER catches all
# three -- in "0e+" the exponent group needs digits after the sign, so NUM stops
# at "0" and the "e" is what gives it away.  The exclusion is deliberately no
# wider than that: forbidding a following "+" or "-" as well would reject
# "exact slope 2.12+/-0.03", a legitimate value-plus-uncertainty format that one
# of the real checks emits.
# The atomic group is load-bearing.  With a plain group the engine BACKTRACKS
# into a shorter number to satisfy the lookahead, so "diff=1.5e-15x" quietly
# yields 0.15 -- a wrong value rather than a refusal, which is the worse
# failure.  (?>...) forbids that: the number is matched maximally or not at all.
NUM = (r"(?>[-+]?(?:\d+\.?\d*(?:[eE][-+]?\d+)?|nan|inf(?:inity)?))"
       r"(?![A-Za-z])")

# key: (relative path from repo root, kind, pattern, bound, label, nexp)
#   pattern MUST contain exactly one capture group, around the field being
#   checked.  That requirement is the fix for a real hole: the previous version
#   matched a keyword and then took the first (or last) number it could find
#   anywhere on the line, so a line reading
#       error=1.0 tolerance=1.0e-12
#   reported the TOLERANCE, because "1.0" carries no exponent and the
#   scientific-notation pattern skipped over it.  The same happened for
#   "error=1" and for "error=missing".  All three passed a 1e-12 bound.
#   Naming the field explicitly removes the guesswork: if the field is absent or
#   unparsable the pattern simply does not match, the sample count falls short of
#   nexp, and the check fails instead of silently reporting its neighbour.
#   The patterns differ per script because the output formats do -- in the
#   horizon dilogarithm the residual is the LAST column, elsewhere it follows the
#   keyword -- and that is exactly why a single generic rule was unsafe.
CHECKS = [
    ("KerrMetric/pipeline_completa_deltaphi.py",
     "residual", rf"differenza\s*=\s*({NUM})", 1e-12,
     "end-to-end (M,a,E,J)->delta phi", 1),
    # Two field names in this script's output, and a prose line "se diff~1e-12"
    # that the old loose pattern parsed as a fifth residual -- a sample invented
    # from a comment.  The alternation must put the longer name first, and "~"
    # is not "=", so the prose no longer matches.
    ("KerrScripts/kerr_psi_explicit_verified.py",
     "residual", rf"(?:differenza|diff)\s*=\s*({NUM})", 1e-11,
     "TK tau reduction psi", 5),
    ("KerrMetric/kerr_tbranch_psi_assembly.py",
     "residual", rf"diff\s*=\s*({NUM})", 1e-10,
     "TK t-branch clock + assembly", 2),
    ("NonStationaryMetrics/VaidyaMetric/vaidya_dMF_reduction.py",
     "residual", rf"diff\s*=\s*({NUM})", 1e-12, "Vaidya reduction c_k^m", 4),
    ("NonStationaryMetrics/VaidyaMetric/vaidya_fully_explicit.py",
     "residual", rf"diff\s*=\s*({NUM})", 1e-11,
     "Vaidya clock v + assembled delta phi", 1),
    ("NonStationaryMetrics/VaidyaMetric/vaidya_tau_assembly.py",
     "residual", rf"differenza\s*=\s*({NUM})", 1e-12,
     "Vaidya tau-branch assembly", 1),
    ("NonStationaryMetrics/VaidyaMetric/vaidya_horizon_dilog.py",
     "residual", rf"^\s*\d+\s*:.*\s({NUM})\s*$", 1e-10,
     "horizon dilogarithm D_k", 5),
    ("NonStationaryMetrics/VaidyaMetric/vaidya_ell_dilog_match.py",
     "residual", rf"diff\s*=\s*({NUM})", 1e-12, "separatrix elliptic U_0", 10),
    ("NonStationaryMetrics/VaidyaMetric/vaidya_asymmetry.py",
     "residual", rf"diff\s*=\s*({NUM})", 1e-13,
     "accretion/evaporation split", 3),
    ("NonStationaryMetrics/ThakurtaMetric/adiabatic_first_order_exact.py",
     "slope", rf"robust\s+slope\s*=\s*({NUM})", (1.7, 2.3),
     "TK t exact PT, delta phi vs true flow", 1),
    ("NonStationaryMetrics/VaidyaMetric/vaidya_first_order_offshell.py",
     "slope", rf"slope\s+full\s*=\s*({NUM})", (1.7, 2.3),
     "Vaidya v exact PT off-shell", 1),
    ("NonStationaryMetrics/ThakurtaMetric/fig_phi_validation_corrected.py",
     "slope", rf"exact\s+slope\s+({NUM})", (1.7, 2.3),
     "true-dynamics validation (t,tau)", 3),
]

# scripts that are slow (symbolic / ODE integration) -- skipped under --quick
SLOW = {
    "kerr_tbranch_psi_assembly.py",
    "vaidya_tau_assembly.py",
    "fig_phi_validation_corrected.py",
    "adiabatic_first_order_exact.py",
}


def all_numbers(text, pattern, floats=False):
    """Return the captured field from every line matching `pattern`.

    `pattern` must carry exactly one capture group around the quantity being
    checked; see the note on CHECKS for why guessing which number on the line is
    the right one was unsafe.  Returns (values, n_matching_lines) so the caller
    can also report lines that matched a keyword without carrying the field.
    """
    hits, n_lines = [], 0
    for line in text.splitlines():
        m = re.search(pattern, line, re.IGNORECASE)
        if not m:
            continue
        n_lines += 1
        if m.lastindex:
            hits.append(float(m.group(1)))
    return hits, n_lines


def run_one(rel, kind, pattern, bound, label, timeout, nexp=None):
    path = ROOT / rel
    if not path.exists():
        return "MISSING", f"{rel} not found", None
    t0 = time.time()
    try:
        p = subprocess.run([sys.executable, str(path)],
                           capture_output=True, text=True, timeout=timeout,
                           cwd=path.parent)
    except subprocess.TimeoutExpired:
        return "TIMEOUT", f"exceeded {timeout}s", time.time() - t0
    dt = time.time() - t0
    if p.returncode != 0:
        tail = p.stderr.strip().splitlines()[-1:] or ["(no stderr)"]
        return "ERROR", f"exit {p.returncode}: {tail[-1]}", dt
    out = p.stdout
    vals, n_lines = all_numbers(out, pattern, floats=(kind == "slope"))
    if not vals:
        return "NOVALUE", f"no '{pattern}' line in output", dt
    # A matching line carrying no number is NOT an error: the per-script
    # patterns are deliberately loose, and prose trips them -- "Kerr" contains
    # "err", "DIFFERENZA" contains "diff", "residuo" contains "resid".  Failing
    # on that flagged three healthy scripts.  What guards against a dropped
    # sample is instead that nan/inf are now parsed rather than skipped, plus
    # the expected-count test below.
    prose = n_lines - len(vals)
    if nexp is not None and len(vals) != nexp:
        # A check that silently reports on fewer samples than the script is
        # supposed to produce is not a check.  Counting them is the cheapest
        # way to notice a configuration that stopped being printed.
        return ("FAIL", f"expected {nexp} sample(s), parsed {len(vals)}", dt)
    bad = [v for v in vals if v != v or abs(v) == float("inf")]
    if bad:
        return "FAIL", f"{len(bad)} non-finite of {len(vals)} matches: {bad[:3]}", dt
    if kind == "residual":
        worst = max(abs(v) for v in vals)          # the largest, not the last
        ok = worst <= bound
        return ("PASS" if ok else "FAIL",
                f"worst residual={worst:.2e} of {len(vals)} "
                f"(bound {bound:.0e}{', '+str(prose)+' prose line(s) skipped' if prose else ''})", dt)
    else:  # slope
        lo, hi = bound
        ok = all(lo <= v <= hi for v in vals)      # every configuration, not the last
        worst = min(vals, key=lambda v: min(abs(v - lo), abs(v - hi)) if not (lo <= v <= hi) else 1e9)
        return ("PASS" if ok else "FAIL",
                f"{len(vals)} slope(s), worst={worst:.3f} (want [{lo},{hi}]{', '+str(prose)+' prose line(s) skipped' if prose else ''})", dt)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="skip slow scripts")
    ap.add_argument("-k", metavar="SUBSTR", help="only checks whose label/path match")
    ap.add_argument("--timeout", type=int, default=1800)
    args = ap.parse_args()

    checks = CHECKS
    if args.quick:
        checks = [c for c in checks if Path(c[0]).name not in SLOW]
    if args.k:
        s = args.k.lower()
        checks = [c for c in checks if s in c[0].lower() or s in c[4].lower()]

    print(f"Regression runner -- {len(checks)} checks (root {ROOT})\n" + "=" * 78)
    n_pass = n_fail = 0
    for rel, kind, pattern, bound, label, nexp in checks:
        status, detail, dt = run_one(rel, kind, pattern, bound, label,
                                     args.timeout, nexp)
        tstr = f"{dt:5.1f}s" if dt else "  -  "
        mark = {"PASS": "✓"}.get(status, "✗")
        print(f"[{mark} {status:7}] {tstr}  {label}")
        print(f"            {detail}")
        if status == "PASS":
            n_pass += 1
        else:
            n_fail += 1
    print("=" * 78)
    print(f"{n_pass} passed, {n_fail} failed / {len(checks)} total")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
