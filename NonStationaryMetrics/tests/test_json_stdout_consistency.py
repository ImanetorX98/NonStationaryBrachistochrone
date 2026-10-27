"""The JSON results and the printed log must report the same numbers (slow: runs the producers).

The runners decide on the JSON values; a reader looks at the printed log. This test runs each of the
12 Python regression producers once with NSB_RESULT_JSON set, reads every sample from the JSON, reads
the field that the legacy stdout parser would read, and requires agreement within the rounding of the
printed token (half a unit in its last digit); printed tokens must be complete finite numbers, and a
printed list must match the JSON samples one to one, in any order.  With --sage-prefix it does the same for the four Sage
producers of sage_env/run_sage_checks.py.

Run from NonStationaryMetrics/:
    python3 tests/test_json_stdout_consistency.py                 # 12 Python producers (a few minutes)
    python3 tests/test_json_stdout_consistency.py --sage-prefix ~/miniconda3/envs/sage-nsb
Exit status 1 on any disagreement, missing value, failed run or invalid JSON.
"""
import argparse
import itertools
import os
import re
import subprocess
import sys
import tempfile
from decimal import Decimal
from pathlib import Path

NSM = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(NSM))
sys.path.insert(0, str(NSM / "sage_env"))
import run_regression as R  # noqa: E402
from verification_results import read_results  # noqa: E402


def half_unit(token):
    """Half a unit in the last printed digit of a numeric token."""
    return float(Decimal(1).scaleb(Decimal(token).as_tuple().exponent)) / 2


def agree(value, token):
    printed = float(token)
    return abs(value - printed) <= half_unit(token) * (1 + 1e-9) + 1e-12 * abs(printed)


def value_agrees(value, token):
    """One printed token: it must be a complete finite number that matches the JSON value."""
    import run_sage_checks as S
    return S.parse_float(token) is not None and agree(value, token)


def list_agrees(values, body):
    """A printed list against JSON samples: every element a complete finite number, same count,
    and a one-to-one match (each printed element used once), in any order."""
    import run_sage_checks as S
    if S.parse_list(body) is None:
        return False
    tokens = [t for t in re.split(r"[,\s]+", body.strip()) if t]
    if len(tokens) != len(values):
        return False
    return any(all(agree(v, t) for v, t in zip(values, perm)) for perm in itertools.permutations(tokens))


def run(cmd, cwd, result):
    env = dict(os.environ, NSB_RESULT_JSON=str(result), MPLBACKEND="Agg")
    return subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=3600, env=env)


def python_producers():
    bad = 0
    for rel, kind, pattern, bound, label, nexp in R.CHECKS:
        path = R.ROOT / rel
        ids = R.RESULT_SAMPLES[path.name]
        with tempfile.TemporaryDirectory() as folder:
            result = Path(folder) / "result.json"
            p = run([sys.executable, path.name], path.parent, result)
            if p.returncode != 0:
                print(f"KO  {rel}: exit {p.returncode}"); bad += 1; continue
            values = read_results(result, rel, {kind: ids})[kind]
        tokens = []
        for line in p.stdout.splitlines():
            m = re.search(pattern, line, re.IGNORECASE)
            if m and m.lastindex:
                tokens.append(m.group(1))
        if len(tokens) != len(ids):
            print(f"KO  {rel}: {len(ids)} JSON samples but {len(tokens)} printed fields"); bad += 1; continue
        wrong = [(sid, values[sid], tok) for sid, tok in zip(ids, tokens) if not value_agrees(values[sid], tok)]
        if wrong:
            print(f"KO  {rel}: JSON and printed values differ: {wrong}"); bad += 1
        else:
            print(f"ok  {rel}: {len(ids)} sample(s) agree with the printed log")
    return bad


def sage_producers(prefix):
    import run_sage_checks as S
    runners = {"python": str(Path(prefix) / "bin" / "python"), "sage": str(Path(prefix) / "bin" / "sage")}
    bad = 0
    for script, kind, items in S.CHECKS:
        rel = script.relative_to(S.ROOT).as_posix()
        with tempfile.TemporaryDirectory() as folder:
            result = Path(folder) / "result.json"
            p = run([runners[kind], script.name], script.parent, result)
            if p.returncode != 0:
                print(f"KO  {rel}: exit {p.returncode}"); bad += 1; continue
            fields = read_results(result, rel, S.RESULT_SPECS[script.name])
        for label, rx, how, bound, metric in items:
            m = re.search(rx, p.stdout)
            if not m:
                print(f"KO  {rel}: {label}: field not printed"); bad += 1; continue
            vals = [fields[metric][sid] for sid in S.RESULT_SPECS[script.name][metric]]
            if how == "pos":   # the printed list and the JSON samples need not share an order
                ok = list_agrees(vals, m.group(1))
            else:
                ok = value_agrees(vals[0], m.group(1))
            print(f"{'ok ' if ok else 'KO '} {rel}: {label}")
            bad += not ok
    return bad


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--sage-prefix", help="conda prefix of the Sage environment (optional)")
    args = ap.parse_args()
    bad = python_producers()
    if args.sage_prefix:
        bad += sage_producers(os.path.expanduser(args.sage_prefix))
    print("ALL VALUES AGREE" if bad == 0 else f"{bad} DISAGREEMENT(S)")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
