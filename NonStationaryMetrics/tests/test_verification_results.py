"""Adversarial JSON transport and runner tests; no scientific runtime required.

Run: python3 tests/test_verification_results.py
These tests establish failure handling, not the scientific identities.
"""
import copy
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

NSM = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(NSM))
import run_regression as R
from verification_results import ResultWriter, read_results, write_json

spec = importlib.util.spec_from_file_location("sage_runner", NSM / "sage_env/run_sage_checks.py")
S = importlib.util.module_from_spec(spec)
spec.loader.exec_module(S)


def document(producer, expected, value=0.0):
    return {"schema_version": 1, "producer": producer, "complete": True,
            "metrics": [{"id": metric, "samples": [{"id": sid, "value": value} for sid in samples]}
                        for metric, samples in expected.items()]}


class TransportTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name) / "result.json"
        self.producer = "a.py"
        self.expected = {"residual": ("a", "b")}
        self.doc = document(self.producer, self.expected)

    def read(self, data):
        self.path.write_text(json.dumps(data))
        return read_results(self.path, self.producer, self.expected)

    def test_valid_reordered_samples(self):
        self.doc["metrics"][0]["samples"].reverse()
        self.assertEqual(self.read(self.doc), {"residual": {"a": 0.0, "b": 0.0}})

    def test_invalid_documents(self):
        variants = []
        for key, value in [("schema_version", 2), ("schema_version", True),
                           ("schema_version", 1.0), ("complete", False),
                           ("complete", 1), ("producer", "other.py"), ("metrics", {})]:
            d = copy.deepcopy(self.doc); d[key] = value; variants.append(d)
        d = copy.deepcopy(self.doc); del d["complete"]; variants.append(d)
        d = copy.deepcopy(self.doc); d["passed"] = True; variants.append(d)
        for value in [True, None, "0.0", float("nan"), float("inf"), 1e999]:
            d = copy.deepcopy(self.doc); d["metrics"][0]["samples"][0]["value"] = value; variants.append(d)
        for key, value in [("samples", []), ("samples", {}), ("id", "other_metric")]:
            d = copy.deepcopy(self.doc); d["metrics"][0][key] = value; variants.append(d)
        d = copy.deepcopy(self.doc); d["metrics"] *= 2; variants.append(d)
        d = copy.deepcopy(self.doc); d["metrics"][0]["samples"] = [d["metrics"][0]["samples"][0]] * 2; variants.append(d)
        d = copy.deepcopy(self.doc); d["metrics"][0]["samples"][0]["id"] = "substituted"; variants.append(d)
        d = copy.deepcopy(self.doc); d["metrics"][0]["samples"].append({"id": "extra", "value": 0}); variants.append(d)
        d = copy.deepcopy(self.doc); d["metrics"][0]["samples"][0]["bound"] = 1; variants.append(d)
        for index, d in enumerate(variants):
            with self.subTest(index=index), self.assertRaises(ValueError):
                self.read(d)

    def test_truncation_duplicate_keys_and_exponent_overflow(self):
        for text in ['{"schema_version":', json.dumps(self.doc).replace('"complete": true', '"complete": true, "complete": true'),
                     json.dumps(self.doc).replace('"value": 0.0', '"value": 1e999', 1)]:
            with self.subTest(text=text), self.assertRaises(ValueError):
                self.path.write_text(text)
                read_results(self.path, self.producer, self.expected)

    def test_absent_file_and_no_completion_before_finish(self):
        with patch.dict(os.environ, NSB_RESULT_JSON=str(self.path)):
            w = ResultWriter(NSM / "run_regression.py")
            w.add("residual", "a", 0.0)
            self.assertFalse(self.path.exists())
            w.finish()
        self.assertEqual(read_results(self.path, "NonStationaryMetrics/run_regression.py", {"residual": ("a",)}),
                         {"residual": {"a": 0.0}})

    def test_writer_rejects_nonfinite_and_duplicate_samples(self):
        w = ResultWriter(NSM / "run_regression.py")
        for value in [True, float("nan"), float("inf")]:
            with self.assertRaises(ValueError):
                w.add("r", "x", value)
        w.add("r", "x", 0)
        with self.assertRaises(ValueError):
            w.add("r", "x", 1)

    def test_failed_serialisation_preserves_existing_file(self):
        write_json(self.path, {"old": True})
        with self.assertRaises(ValueError):
            write_json(self.path, {"bad": float("nan")})
        self.assertEqual(json.loads(self.path.read_text()), {"old": True})


class RunnerTests(unittest.TestCase):
    def python_run(self, value=0.0, mutate=None, code=0, mode="write", kind="residual"):
        check = next(c for c in R.CHECKS if c[1] == kind)
        rel, _, pattern, bound, label, nexp = check
        def fake_run(cmd, **kw):
            if mode == "timeout":
                raise subprocess.TimeoutExpired(cmd, 1)
            if mode == "spawn_error":
                raise OSError("cannot execute")
            d = document(rel, {kind: R.RESULT_SAMPLES[Path(rel).name]}, value)
            if mutate:
                mutate(d)
            if mode == "write":
                Path(kw["env"]["NSB_RESULT_JSON"]).write_text(json.dumps(d))
            # A misleading human success line must never override JSON/exit status.
            return subprocess.CompletedProcess(cmd, code, stdout="ALL CHECKS PASSED difference = 0", stderr="failed")
        with patch.object(R.subprocess, "run", fake_run):
            return R.run_one(rel, kind, pattern, bound, label, 10, nexp)[0]

    def test_python_uses_json_not_stdout(self):
        self.assertEqual(self.python_run(), "PASS")
        self.assertEqual(self.python_run(value=1), "FAIL")
        self.assertEqual(self.python_run(value=2.0, kind="slope"), "PASS")
        self.assertEqual(self.python_run(value=1.0, kind="slope"), "FAIL")

    def test_python_failures_do_not_fall_back(self):
        for mode, code, want in [("missing", 0, "INVALID"), ("write", 1, "ERROR"),
                                  ("timeout", 0, "TIMEOUT"), ("spawn_error", 0, "ERROR")]:
            with self.subTest(mode=mode, code=code):
                self.assertEqual(self.python_run(code=code, mode=mode), want)
        self.assertEqual(self.python_run(mutate=lambda d: d.update(complete=False)), "INVALID")
        self.assertEqual(self.python_run(value=float("nan")), "INVALID")

    def sage_run(self, mutate=None, mode="write", code=0):
        with tempfile.TemporaryDirectory() as folder:
            pre = Path(folder); (pre / "bin").mkdir()
            for name in ("sage", "python"):
                (pre / "bin" / name).touch()
            report = pre / "summary.json"
            def fake_run(cmd, **kw):
                name = cmd[-1]
                if mode == "timeout":
                    raise subprocess.TimeoutExpired(cmd, 1)
                script = next(c[0] for c in S.CHECKS if c[0].name == name)
                d = document(script.relative_to(S.ROOT).as_posix(), S.RESULT_SPECS[name], 1e-12)
                for m in d["metrics"]:
                    if m["id"] == "positive_eigenvalues":
                        for item in m["samples"]:
                            item["value"] = 0.5
                if mutate:
                    mutate(d)
                if mode == "write":
                    Path(kw["env"]["NSB_RESULT_JSON"]).write_text(json.dumps(d))
                return subprocess.CompletedProcess(cmd, code, stdout="ALL CHECKS PASSED", stderr="failed")
            with patch.object(S.subprocess, "run", fake_run), patch.object(sys, "argv", ["runner", "--prefix", folder, "--report-json", str(report)]):
                result = S.main()
            return result, json.loads(report.read_text())

    def test_sage_json_without_matching_prose(self):
        code, report = self.sage_run()
        self.assertEqual(code, 0)
        self.assertTrue(report["passed"])

    def test_sage_failures(self):
        for mode, code in [("missing", 0), ("timeout", 0), ("write", 1)]:
            with self.subTest(mode=mode, code=code):
                rc, report = self.sage_run(mode=mode, code=code)
                self.assertEqual(rc, 1)
                self.assertFalse(report["passed"])
        for mutate in [lambda d: d.update(schema_version=999),
                       lambda d: d["metrics"][0]["samples"][0].update(value=1),
                       lambda d: d["metrics"][0]["samples"][0].update(value=float("nan")),
                       lambda d: d["metrics"][0]["samples"][0].update(id="substituted"),
                       lambda d: d["metrics"].extend(copy.deepcopy(d["metrics"]))]:
            with self.subTest(mutate=mutate):
                rc, report = self.sage_run(mutate=mutate)
                self.assertEqual(rc, 1)
                self.assertFalse(report["passed"])

    def test_negative_eigenvalue_fails(self):
        def mutate(d):
            for m in d["metrics"]:
                if m["id"] == "positive_eigenvalues":
                    m["samples"][0]["value"] = -0.1
        code, report = self.sage_run(mutate=mutate)
        self.assertEqual(code, 1)
        self.assertFalse(report["passed"])


if __name__ == "__main__":
    unittest.main()
