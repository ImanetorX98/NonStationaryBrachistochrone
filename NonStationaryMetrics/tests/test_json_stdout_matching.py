"""Fast cases for the matching rules of test_json_stdout_consistency.py (no producer is run).

A printed list must contain only complete finite numbers and match the JSON samples one to one, in
any order; a single printed value must be a complete finite number within its printed rounding.
Run from NonStationaryMetrics/:  python3 tests/test_json_stdout_matching.py   (exit status 1 on failure)
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_json_stdout_consistency import list_agrees, value_agrees  # noqa: E402

CASES = [
    # (description, function, JSON values, printed text, expected)
    ("list in the same order", list_agrees, [0.1, 1.7], "0.1, 1.7", True),
    ("list in reversed order", list_agrees, [0.1, 1.7], "1.7, 0.1", True),
    ("equal multiplicities", list_agrees, [0.1, 0.1], "0.1, 0.1", True),
    ("multiplicities differ (one printed value reused)", list_agrees, [0.1, 0.1], "0.1, 1.7", False),
    ("junk suffix on an element", list_agrees, [0.1, 1.7], "0.1junk, 1.7", False),
    ("extra nan element", list_agrees, [0.1, 1.7], "nan, 0.1, 1.7", False),
    ("fewer printed elements", list_agrees, [0.1, 1.7], "0.1", False),
    ("value outside the printed rounding", list_agrees, [0.1, 1.76], "0.1, 1.7", False),
    ("empty list", list_agrees, [0.1], "", False),
    ("single value within rounding", value_agrees, [2.2702e-12], "2.270e-12", True),
    ("single value outside rounding", value_agrees, [2.2712e-12], "2.270e-12", False),
    ("single nan token", value_agrees, [0.0], "nan", False),
    ("single token with suffix", value_agrees, [1.0], "1.0x", False),
]


def main():
    bad = 0
    for label, fn, values, text, want in CASES:
        got = fn(values, text) if fn is list_agrees else fn(values[0], text)
        ok = got == want
        bad += not ok
        print(f"{'ok ' if ok else 'KO '} {label:50s} -> {got}")
    print("ALL CASES AS EXPECTED" if bad == 0 else f"{bad} CASE(S) WRONG")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
