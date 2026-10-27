#!/usr/bin/env python3
"""Run the four Sage scripts cited in Paper II and check their JSON values.

    python3 run_sage_checks.py [--env sage-nsb | --prefix /path/to/conda/env]

The runner reads versioned, finite, identified JSON measurements and applies
its own tolerances. --format legacy retains stdout parsing for comparisons.
A missing, incomplete or invalid JSON result fails without a stdout fallback.
The tolerances are the orders of magnitude quoted in Paper II, with a margin for
platform-dependent rounding; they are not the printed values themselves.

`.sage` files are run with the environment's `sage`; the already preparsed
`.sage.py` files with the environment's Python (conda's `sage` has no `-python`).
"""
from __future__ import annotations

import argparse
import math
import os
import re
import subprocess
import sys
import time
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
NSM = HERE.parent
ROOT = NSM.parent
sys.path.insert(0, str(NSM))
from verification_results import read_results, write_json

RESULT_SPECS = {
    "kerr_tau_general_genus2_periods.sage.py": {
        "symmetry": ("norm",), "positive_eigenvalues": ("eigenvalue_0", "eigenvalue_1")},
    "kerr_tbranch_genus2_periods.sage.py": {
        "symmetry": ("norm",), "positive_eigenvalues": ("eigenvalue_0", "eigenvalue_1")},
    "kerr_quasiperiods_bel.sage": {"symmetry": ("kappa_norm",)},
    "kerr_thirdkind_theta_closed.sage": {"relative_residual": ("best_characteristic",)},
}
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
        ("Riemann symmetry |tau-tau^T| (paper: 7e-12)", r"\|tau-tau\^T\| = " + FIELD, "max", 1e-11, "symmetry"),
        ("Im(tau) positive definite", r"eigenvalues of Im\(tau\) = \[([^\]]*)\]", "pos", 2, "positive_eigenvalues"),
    ]),
    (ROOT / "KerrScripts/kerr_tbranch_genus2_periods.sage.py", "python", [
        ("Riemann symmetry |tau-tau^T| (paper: 3e-12)", r"\|tau - tau\^T\| = " + FIELD, "max", 1e-11, "symmetry"),
        ("Im(tau) positive definite", r"eigenvalues of Im\(tau\) = \[([^\]]*)\]", "pos", 2, "positive_eigenvalues"),
    ]),
    (ROOT / "KerrScripts/kerr_quasiperiods_bel.sage", "sage", [
        ("quasi-period symmetry of kappa (paper: 1e-12)", r"kappa SYM ERR\s*=\s*" + FIELD, "max", 1e-11, "symmetry"),
    ]),
    (ROOT / "KerrScripts/kerr_thirdkind_theta_closed.sage", "sage", [
        ("third-kind theta closed form, best relative residual (paper: 7e-5)",
         r"BEST: delta #\d+, relative residual " + FIELD, "max", 1e-4, "relative_residual"),
    ]),
]

# Each check names the JSON metric it reads; the names must cover RESULT_SPECS exactly,
# so that reordering either table cannot silently swap checks and bounds.
for _script, _kind, _items in CHECKS:
    if sorted(item[4] for item in _items) != sorted(RESULT_SPECS[_script.name]):
        raise SystemExit(f"run_sage_checks: metric names of {_script.name} do not match RESULT_SPECS")


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
    ap.add_argument("--format", choices=("json", "legacy"), default="json")
    ap.add_argument("--report-json", type=Path, help="write a structured run summary")
    args = ap.parse_args()
    pre = env_prefix(args)
    runners = {"python": str(pre / "bin" / "python"), "sage": str(pre / "bin" / "sage")}
    for r in runners.values():
        if not Path(r).exists():
            print(f"missing {r}: build the environment with setup_sage_env.sh")
            return 1

    n_fail = n_pass = 0
    records = []
    for k, (script, kind, items) in enumerate(CHECKS, 1):
        rel = script.relative_to(ROOT)
        print(f"[{k}/{len(CHECKS)}] running {rel} ...", flush=True)
        t0 = time.time()
        record = {"producer": rel.as_posix(), "checks": [], "status": "ERROR"}
        records.append(record)
        try:
            with tempfile.TemporaryDirectory(prefix="nsb-sage-result-") as folder:
                result_path = Path(folder) / "result.json"
                env = dict(os.environ)
                env.pop("NSB_RESULT_JSON", None)
                fields = None
                if args.format == "json":
                    env["NSB_RESULT_JSON"] = str(result_path)
                p = subprocess.run([runners[kind], script.name], cwd=script.parent,
                                   capture_output=True, text=True, timeout=args.timeout, env=env)
                if p.returncode == 0 and args.format == "json":
                    fields = read_results(result_path, rel.as_posix(), RESULT_SPECS[script.name])
            record["seconds"] = time.time() - t0
            print(f"      finished in {time.time() - t0:.0f} s", flush=True)
        except subprocess.TimeoutExpired:
            record.update(status="TIMEOUT", detail=f"timeout {args.timeout}s")
            print(f"[FAIL] {rel}: timeout {args.timeout}s"); n_fail += 1; continue
        except (OSError, ValueError, KeyError, OverflowError) as exc:
            record.update(status="INVALID", detail=str(exc))
            print(f"[FAIL] {rel}: invalid result: {exc}"); n_fail += 1; continue
        if p.returncode != 0:
            tail = (p.stderr.strip().splitlines() or ["(no stderr)"])[-1]
            record["detail"] = f"exit {p.returncode}: {tail}"
            print(f"[FAIL] {rel}: exit {p.returncode}: {tail}"); n_fail += 1; continue
        record["status"] = "PASS"
        for label, rx, how, bound, metric in items:
            if fields is not None:
                vals = [fields[metric][sid] for sid in RESULT_SPECS[script.name][metric]]
                v = vals[0]
            else:
                m = re.search(rx, p.stdout)
                if not m:
                    record["status"] = "FAIL"
                    record["checks"].append({"label": label, "passed": False, "detail": "value not printed"})
                    print(f"[FAIL] {rel}: {label}: value not printed"); n_fail += 1; continue
                v = parse_float(m.group(1)) if how == "max" else None
                vals = parse_list(m.group(1)) if how == "pos" else None
            if how == "max":
                ok = v is not None and abs(v) <= bound
                shown = f"{v:.3e}" if v is not None else f"invalid token {m.group(1)!r}"
                print(f"[{'PASS' if ok else 'FAIL'}] {rel}: {label}: {shown} (bound {bound:.0e})")
            else:
                ok = vals is not None and len(vals) == bound and all(x > 0 for x in vals)
                shown = vals if vals is not None else f"invalid list {m.group(1)!r}"
                print(f"[{'PASS' if ok else 'FAIL'}] {rel}: {label}: {shown} (expect {bound} finite > 0)")
            n_fail += not ok
            n_pass += ok
            record["checks"].append({"label": label, "kind": how, "bound": bound,
                                      "values": [v] if how == "max" else vals, "passed": ok})
            if not ok:
                record["status"] = "FAIL"
    n_total = sum(len(c[2]) for c in CHECKS)
    print(f"{n_pass}/{n_total} checks passed")
    if args.report_json:
        write_json(args.report_json, {"schema_version": 1, "runner": "sage-checks",
                                     "format": args.format, "complete": True,
                                     "passed": n_fail == 0, "checks": records})
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
