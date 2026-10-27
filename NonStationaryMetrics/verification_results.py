"""Versioned, strict JSON transport for numerical verification results.

Producers optionally write to NSB_RESULT_JSON, after successful completion.
Consumers supply their own expected metric/sample IDs and acceptance bounds.
This module uses only Python's standard library (also available inside Sage).
"""
from __future__ import annotations

import json
import math
import os
import tempfile
from pathlib import Path

SCHEMA_VERSION = 1
ROOT = Path(__file__).resolve().parent.parent


def write_json(path, document):
    """Atomically replace a JSON document; never serialise NaN or infinity."""
    path = Path(path)
    text = json.dumps(document, indent=2, allow_nan=False) + "\n"
    with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=path.parent,
                                     prefix=path.name + ".", delete=False) as f:
        temp = Path(f.name)
        try:
            f.write(text)
            f.flush()
            os.fsync(f.fileno())
        except BaseException:
            temp.unlink(missing_ok=True)
            raise
    try:
        os.replace(temp, path)
    finally:
        temp.unlink(missing_ok=True)


class ResultWriter:
    """Collect numerical fields, preserving the human-readable stdout.

    Call finish() explicitly at the end of a successful program, never through
    atexit: an exception or early exit must not leave a completed document.
    Values are finite binary64 numbers, not arbitrary-precision certificates.
    """

    def __init__(self, script):
        self.producer = Path(script).resolve().relative_to(ROOT).as_posix()
        self.metrics = {}

    def add(self, metric, sample, value):
        if not isinstance(metric, str) or not metric or not isinstance(sample, str) or not sample:
            raise ValueError("metric and sample IDs must be non-empty strings")
        if isinstance(value, bool):
            raise ValueError("boolean result is not a numerical measurement")
        value = float(value)
        if not math.isfinite(value):
            raise ValueError("non-finite verification result")
        samples = self.metrics.setdefault(metric, {})
        if sample in samples:
            raise ValueError("duplicate sample: " + metric + "/" + sample)
        samples[sample] = value

    def finish(self):
        path = os.environ.get("NSB_RESULT_JSON")
        if path:
            write_json(path, {
                "schema_version": SCHEMA_VERSION,
                "producer": self.producer,
                "complete": True,
                "metrics": [
                    {"id": key, "samples": [{"id": k, "value": v} for k, v in values.items()]}
                    for key, values in self.metrics.items()
                ],
            })


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("duplicate JSON key: " + key)
        result[key] = value
    return result


def _keys(obj, expected):
    if not isinstance(obj, dict) or set(obj) != set(expected):
        raise ValueError("unexpected or missing JSON fields; expected " + str(expected))


def read_results(path, producer, expected):
    """Read only the exact producer, schema, metrics and sample sets requested.

    expected maps metric IDs to sample IDs; order is irrelevant. No fallback
    to stdout is permitted when JSON is absent, malformed or incomplete.
    """
    def reject_constant(token):
        raise ValueError("non-standard JSON number: " + token)

    data = json.loads(Path(path).read_text(encoding="utf-8"),
                      object_pairs_hook=_unique_object, parse_constant=reject_constant)
    _keys(data, ("schema_version", "producer", "complete", "metrics"))
    if type(data["schema_version"]) is not int or data["schema_version"] != SCHEMA_VERSION:
        raise ValueError("unsupported result schema version")
    if data["producer"] != producer or data["complete"] is not True:
        raise ValueError("wrong producer or incomplete execution")
    if not isinstance(data["metrics"], list):
        raise ValueError("metrics must be a list")
    result = {}
    for metric in data["metrics"]:
        _keys(metric, ("id", "samples"))
        key = metric["id"]
        if not isinstance(key, str) or key not in expected or key in result:
            raise ValueError("unknown or duplicate metric ID")
        if not isinstance(metric["samples"], list):
            raise ValueError("samples must be a list")
        values = {}
        for sample in metric["samples"]:
            _keys(sample, ("id", "value"))
            sid, value = sample["id"], sample["value"]
            if not isinstance(sid, str) or sid not in expected[key] or sid in values:
                raise ValueError("unknown or duplicate sample ID")
            if type(value) not in (int, float) or not math.isfinite(value):
                raise ValueError("sample value must be a finite JSON number")
            values[sid] = value
        if set(values) != set(expected[key]):
            raise ValueError("missing expected samples for " + key)
        result[key] = values
    if set(result) != set(expected):
        raise ValueError("missing expected metrics")
    return result
