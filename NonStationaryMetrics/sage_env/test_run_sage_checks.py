#!/usr/bin/env python3
"""Adversarial test of run_sage_checks.py: malformed output must never PASS.

    python3 test_run_sage_checks.py

No Sage is needed: the four scripts are replaced by simulated outputs (fixtures),
so this tests the parser and the failure handling, not the mathematics.  The
valid fixture reproduces the values of the real run of 21 October 2026.
"""
from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_sage_checks as R  # noqa: E402

VALID = {
    "kerr_tau_general_genus2_periods.sage.py":
        "[ 7.4s] |tau-tau^T| = 7.29e-12 (simmetrica); autovalori Im(tau) = "
        "[0.18961075409624165, 1.7573154131112518] (>0)\n",
    "kerr_tbranch_genus2_periods.sage.py":
        "|tau - tau^T| = 2.27e-12\neigenvalues of Im(tau) = [0.235, 1.565]  (basis-dependent)\n",
    "kerr_quasiperiods_bel.sage": "kappa SYM ERR      = 2.651609927327309e-12\n",
    "kerr_thirdkind_theta_closed.sage":
        "MIGLIORE: delta #1, residuo rel 7.25e-05, c0=23.4368-1.9615j\n",
}


def run_with(outputs: dict, rc: dict | None = None, timeout: set | None = None) -> int:
    rc, timeout = rc or {}, timeout or set()

    def fake_run(cmd, cwd, capture_output, text, timeout_=None, **kw):
        name = cmd[-1]
        if name in timeout:
            raise subprocess.TimeoutExpired(cmd, 1)
        return SimpleNamespace(returncode=rc.get(name, 0), stdout=outputs.get(name, ""),
                               stderr="simulated error" if rc.get(name) else "")

    def fake_run_kw(cmd, **kw):
        return fake_run(cmd, kw.get("cwd"), True, True)

    with tempfile.TemporaryDirectory() as d:
        (Path(d) / "bin").mkdir()
        for exe in ("python", "sage"):
            (Path(d) / "bin" / exe).write_text("")
        old_run, old_argv = R.subprocess.run, sys.argv
        R.subprocess.run, sys.argv = fake_run_kw, ["x", "--prefix", d]
        try:
            return R.main()
        finally:
            R.subprocess.run, sys.argv = old_run, old_argv


def with_(name, text):
    out = dict(VALID)
    out[name] = text
    return out


TAU, TB, KAP, TH = list(VALID)
CASES = [
    ("valid output", VALID, {}, set(), 0),
    ("eigenvalues [nan 1.7573]", with_(TB, "|tau - tau^T| = 2.27e-12\neigenvalues of Im(tau) = [nan 1.7573]\n"), {}, set(), 1),
    ("one eigenvalue [1.7573]", with_(TB, "|tau - tau^T| = 2.27e-12\neigenvalues of Im(tau) = [1.7573]\n"), {}, set(), 1),
    ("infinite eigenvalue [0.1896 1e999]", with_(TB, "|tau - tau^T| = 2.27e-12\neigenvalues of Im(tau) = [0.1896 1e999]\n"), {}, set(), 1),
    ("negative eigenvalue", with_(TB, "|tau - tau^T| = 2.27e-12\neigenvalues of Im(tau) = [0.1896, -1.7573]\n"), {}, set(), 1),
    ("three eigenvalues", with_(TB, "|tau - tau^T| = 2.27e-12\neigenvalues of Im(tau) = [0.1, 0.2, 0.3]\n"), {}, set(), 1),
    ("empty eigenvalue list", with_(TB, "|tau - tau^T| = 2.27e-12\neigenvalues of Im(tau) = []\n"), {}, set(), 1),
    ("kappa 2.652e-12junk", with_(KAP, "kappa SYM ERR      = 2.652e-12junk\n"), {}, set(), 1),
    ("kappa 0e+", with_(KAP, "kappa SYM ERR      = 0e+\n"), {}, set(), 1),
    ("kappa nan", with_(KAP, "kappa SYM ERR      = nan\n"), {}, set(), 1),
    ("kappa above bound", with_(KAP, "kappa SYM ERR      = 3e-10\n"), {}, set(), 1),
    ("theta residual 7.25e-05x", with_(TH, "MIGLIORE: delta #1, residuo rel 7.25e-05x, c0=1\n"), {}, set(), 1),
    ("symmetry inf", with_(TAU, "|tau-tau^T| = inf (simmetrica); autovalori Im(tau) = [0.1, 1.7]\n"), {}, set(), 1),
    ("field missing", with_(KAP, "nothing here\n"), {}, set(), 1),
    ("non-zero exit", VALID, {KAP: 1}, set(), 1),
    ("timeout", VALID, {}, {TH}, 1),
]


def main() -> int:
    bad = 0
    for label, outputs, rc, to, want in CASES:
        got = run_with(outputs, rc, to)
        ok = got == want
        bad += not ok
        print(f"{'ok ' if ok else 'KO '} {label:38s} exit {got} (expected {want})", flush=True)
    print("ALL CASES AS EXPECTED" if bad == 0 else f"{bad} CASE(S) WRONG")
    return 0 if bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
