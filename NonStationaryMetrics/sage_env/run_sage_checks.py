#!/usr/bin/env python3
"""Run the four Sage scripts cited in Paper II and check their printed values.

    python3 run_sage_checks.py [--env sage-nsb | --prefix /path/to/conda/env]

The scripts themselves print numbers and do not all assert a tolerance, so this
runner does: it parses each quantity quoted in the paper and fails (exit 1) if a
run errors, times out, does not print the quantity, or exceeds the tolerance.
The tolerances are the orders of magnitude quoted in Paper II, with a margin for
platform-dependent rounding; they are not the printed values themselves.

`.sage` files are run with the environment's `sage`; the already preparsed
`.sage.py` files with the environment's Python (conda's `sage` has no `-python`).
"""
from __future__ import annotations

import argparse
import math
import re
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
NSM = HERE.parent
ROOT = NSM.parent
# A field is captured whole (up to whitespace, comma or semicolon) and then validated
# as a COMPLETE token: "2.652e-12junk" or "0e+" must fail, not be read as a prefix.
FIELD = r"([^\s,;]+)"
FLOAT = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")


def parse_float(tok: str):
    """Return the finite float of a complete numeric token, else None."""
    if not FLOAT.fullmatch(tok):
        return None
    v = float(tok)
    return v if math.isfinite(v) else None


def parse_list(body: str):
    """Every token of a bracketed list must be a finite number; else None."""
    toks = [t for t in re.split(r"[,\s]+", body.strip()) if t]
    vals = [parse_float(t) for t in toks]
    return None if not toks or any(v is None for v in vals) else vals

# (script, runner, [(label, regex with one group, kind, bound)])
#   kind "max": one finite value, |value| <= bound;
#   kind "pos": exactly `bound` finite values (the genus), all > 0
CHECKS = [
    (NSM / "KerrSessionScripts/kerr_tau_general_genus2_periods.sage.py", "python", [
        ("Riemann symmetry |tau-tau^T| (paper: 7e-12)", r"\|tau-tau\^T\| = " + FIELD, "max", 1e-11),
        ("Im(tau) positive definite", r"eigenvalues of Im\(tau\) = \[([^\]]*)\]", "pos", 2),
    ]),
    (ROOT / "KerrScripts/kerr_tbranch_genus2_periods.sage.py", "python", [
        ("Riemann symmetry |tau-tau^T| (paper: 3e-12)", r"\|tau - tau\^T\| = " + FIELD, "max", 1e-11),
        ("Im(tau) positive definite", r"eigenvalues of Im\(tau\) = \[([^\]]*)\]", "pos", 2),
    ]),
    (ROOT / "KerrScripts/kerr_quasiperiods_bel.sage", "sage", [
        ("quasi-period symmetry of kappa (paper: 1e-12)", r"kappa SYM ERR\s*=\s*" + FIELD, "max", 1e-11),
    ]),
    (ROOT / "KerrScripts/kerr_thirdkind_theta_closed.sage", "sage", [
        ("third-kind theta closed form, best relative residual (paper: 7e-5)",
         r"BEST: delta #\d+, relative residual " + FIELD, "max", 1e-4),
    ]),
]


def env_prefix(args) -> Path:
    if args.prefix:
        return Path(args.prefix)
    base = subprocess.run(["conda", "info", "--base"], capture_output=True, text=True,
                          check=True).stdout.strip()
    return Path(base) / "envs" / args.env


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--env", default="sage-nsb", help="conda environment name")
    ap.add_argument("--prefix", help="conda environment prefix (overrides --env)")
    ap.add_argument("--timeout", type=int, default=3600)
    args = ap.parse_args()
    pre = env_prefix(args)
    runners = {"python": str(pre / "bin" / "python"), "sage": str(pre / "bin" / "sage")}
    for r in runners.values():
        if not Path(r).exists():
            print(f"missing {r}: build the environment with setup_sage_env.sh")
            return 1

    n_fail = n_pass = 0
    for k, (script, kind, items) in enumerate(CHECKS, 1):
        rel = script.relative_to(ROOT)
        print(f"[{k}/{len(CHECKS)}] running {rel} ...", flush=True)
        t0 = time.time()
        try:
            p = subprocess.run([runners[kind], script.name], cwd=script.parent,
                               capture_output=True, text=True, timeout=args.timeout)
            print(f"      finished in {time.time() - t0:.0f} s", flush=True)
        except subprocess.TimeoutExpired:
            print(f"[FAIL] {rel}: timeout {args.timeout}s"); n_fail += 1; continue
        if p.returncode != 0:
            tail = (p.stderr.strip().splitlines() or ["(no stderr)"])[-1]
            print(f"[FAIL] {rel}: exit {p.returncode}: {tail}"); n_fail += 1; continue
        for label, rx, how, bound in items:
            m = re.search(rx, p.stdout)
            if not m:
                print(f"[FAIL] {rel}: {label}: value not printed"); n_fail += 1; continue
            if how == "max":
                v = parse_float(m.group(1))
                ok = v is not None and abs(v) <= bound
                shown = f"{v:.3e}" if v is not None else f"invalid token {m.group(1)!r}"
                print(f"[{'PASS' if ok else 'FAIL'}] {rel}: {label}: {shown} (bound {bound:.0e})")
            else:
                vals = parse_list(m.group(1))
                ok = vals is not None and len(vals) == bound and all(x > 0 for x in vals)
                shown = vals if vals is not None else f"invalid list {m.group(1)!r}"
                print(f"[{'PASS' if ok else 'FAIL'}] {rel}: {label}: {shown} (expect {bound} finite > 0)")
            n_fail += not ok
            n_pass += ok
    n_total = sum(len(c[2]) for c in CHECKS)
    print(f"{n_pass}/{n_total} checks passed")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
