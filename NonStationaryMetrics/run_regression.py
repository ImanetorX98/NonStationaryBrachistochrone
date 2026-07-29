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

The figure of merit is pulled from stdout with a per-script regex. For
residual checks we take the LAST scientific-notation number on a matching
line and require |x| <= tol. For slope checks we require the value to lie
in [1.7, 2.3] (second-order canonical perturbation theory).
"""
import argparse
import re
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent  # repo root (Brachistocrona/)

SCI = r"[-+]?\d+\.?\d*[eE][-+]?\d+"

# key: (relative path from repo root, kind, pattern, bound, label)
#   kind "residual": last SCI number on a line matching pattern must be <= bound
#   kind "slope":    last float on a line matching pattern must be within bound (lo,hi)
CHECKS = [
    ("KerrMetric/pipeline_completa_deltaphi.py",
     "residual", r"differenza", 1e-12, "end-to-end (M,a,E,J)->delta phi"),
    ("KerrScripts/kerr_psi_explicit_verified.py",
     "residual", r"(?:diff|resid|err)", 1e-11, "TK tau reduction psi"),
    ("KerrMetric/kerr_tbranch_psi_assembly.py",
     "residual", r"(?:diff|resid|err)", 1e-10, "TK t-branch clock + assembly"),
    ("NonStationaryMetrics/VaidyaMetric/vaidya_dMF_reduction.py",
     "residual", r"(?:diff|resid|err)", 1e-12, "Vaidya reduction c_k^m"),
    ("NonStationaryMetrics/VaidyaMetric/vaidya_fully_explicit.py",
     "residual", r"(?:diff|resid|err)", 1e-11, "Vaidya clock v + assembled delta phi"),
    ("NonStationaryMetrics/VaidyaMetric/vaidya_tau_assembly.py",
     "residual", r"(?:diff|resid|err)", 1e-12, "Vaidya tau-branch assembly"),
    ("NonStationaryMetrics/VaidyaMetric/vaidya_horizon_dilog.py",
     "residual", r"^\s*\d+\s*:", 1e-10, "horizon dilogarithm D_k"),
    ("NonStationaryMetrics/VaidyaMetric/vaidya_ell_dilog_match.py",
     "residual", r"(?:diff|resid|err)", 1e-12, "separatrix elliptic U_0"),
    ("NonStationaryMetrics/VaidyaMetric/vaidya_asymmetry.py",
     "residual", r"(?:diff|resid|err)", 1e-13, "accretion/evaporation split"),
    ("NonStationaryMetrics/ThakurtaMetric/adiabatic_first_order_exact.py",
     "slope", r"robust\s+slope", (1.7, 2.3), "TK t exact PT, delta phi vs true flow"),
    ("NonStationaryMetrics/VaidyaMetric/vaidya_first_order_offshell.py",
     "slope", r"slope\s+full", (1.7, 2.3), "Vaidya v exact PT off-shell"),
    ("NonStationaryMetrics/ThakurtaMetric/fig_phi_validation_corrected.py",
     "slope", r"exact\s+slope", (1.7, 2.3), "true-dynamics validation (t,tau)"),
]

# scripts that are slow (symbolic / ODE integration) -- skipped under --quick
SLOW = {
    "kerr_tbranch_psi_assembly.py",
    "vaidya_tau_assembly.py",
    "fig_phi_validation_corrected.py",
    "adiabatic_first_order_exact.py",
}


def last_number(text, pattern, floats=False):
    """Return a numeric token from the last line matching `pattern`.

    For residual checks (floats=False) take the last sci-notation number on
    the line. For slope checks (floats=True) take the FIRST number appearing
    after the matched keyword, so 'exact slope 2.12+/-0.03' yields 2.12.
    """
    num = SCI if not floats else r"[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?"
    hit = None
    for line in text.splitlines():
        m = re.search(pattern, line, re.IGNORECASE)
        if not m:
            continue
        tail = line[m.end():]
        nums = re.findall(num, tail if floats else line)
        if nums:
            hit = float(nums[0] if floats else nums[-1])
    return hit


def run_one(rel, kind, pattern, bound, label, timeout):
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
    val = last_number(out, pattern, floats=(kind == "slope"))
    if val is None:
        return "NOVALUE", f"no '{pattern}' line in output", dt
    if kind == "residual":
        ok = abs(val) <= bound
        return ("PASS" if ok else "FAIL",
                f"residual={val:.2e} (bound {bound:.0e})", dt)
    else:  # slope
        lo, hi = bound
        ok = lo <= val <= hi
        return ("PASS" if ok else "FAIL",
                f"slope={val:.3f} (want [{lo},{hi}])", dt)


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
    for rel, kind, pattern, bound, label in checks:
        status, detail, dt = run_one(rel, kind, pattern, bound, label, args.timeout)
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
