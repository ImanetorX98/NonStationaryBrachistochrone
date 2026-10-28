# Structured numerical verification results

The Python regression runner and the Sage runner read JSON results by default.
Human-readable stdout remains available, but changing its language, spacing or
layout does not change which numerical fields the runners validate.

## Running the checks

From `NonStationaryMetrics/`, with the interpreters documented in the README:

```bash
python3 run_regression.py --report-json /tmp/nsb-regression.json
python3 sage_env/run_sage_checks.py --env sage-nsb --report-json /tmp/nsb-sage.json
python3 tests/test_verification_results.py
python3 tests/test_producer_header.py           # producer header; standalone use
python3 tests/test_json_stdout_matching.py       # fast: the matching rules of that comparison
python3 tests/test_json_stdout_consistency.py    # slow: JSON values against the printed logs
```

Each runner launches a producer with `NSB_RESULT_JSON` pointing to a fresh,
private temporary file. It requires exit status zero and a valid completed
document, then applies its own acceptance bounds and exact expected sample IDs.
The temporary producer results are removed after use. `--report-json` optionally
retains the measurements, bounds and verdicts in a summary; its parent directory
must already exist. Producer documents and runner summaries are distinct formats.

To retain a producer result directly, for example:

```bash
NSB_RESULT_JSON=/tmp/nsb-assembly.json python3 ../KerrScripts/deltaphi_full_pipeline.py
```

The producers still print their usual output when this variable is absent, also
when a script is copied out of the repository: without `verification_results.py`
it runs as before and exports nothing. If `NSB_RESULT_JSON` is set but the module
cannot be found, the producer stops with a clear message instead of losing the
result. `tests/test_producer_header.py` checks both cases and that the header block
is identical in all producers.
The two `.sage` producers use the source filename supplied by `sage file.sage`;
the Sage runner invokes that supported CLI entry point, rather than `load()`.

The runners decide on the JSON values; the printed log is what a reader sees.
`tests/test_json_stdout_consistency.py` runs the producers and requires every JSON
sample to match the field the legacy parser reads, within the rounding of the
printed digits (with `--sage-prefix` also for the Sage producers). Printed tokens
must be complete finite numbers, and a printed list (the eigenvalues of Im tau)
must match the JSON samples one to one, in any order; the fast cases of
`tests/test_json_stdout_matching.py` cover reordering, equal and unequal
multiplicities, suffixed tokens and extra `nan` elements. In the Sage
runner each check names the JSON metric it reads, and the names must cover the
expected metrics exactly, so reordering either table cannot swap checks and bounds.

## Producer format

```json
{
  "schema_version": 1,
  "producer": "KerrScripts/deltaphi_full_pipeline.py",
  "complete": true,
  "metrics": [
    {
      "id": "residual",
      "samples": [{"id": "assembly", "value": 5.3290705182007514e-14}]
    }
  ]
}
```

The example illustrates the format; the value is not a stored reference used
to decide whether a run passes. Actual results are computed on each execution.

`verification_results.py` implements the standard-library writer and strict
reader; `schemas/verification-result-v1.schema.json` documents the structural
schema. The reader additionally rejects duplicate object keys, duplicate metric
or sample IDs, unknown/missing fields, wrong producer names, unsupported versions,
non-finite values, booleans and strings used as numbers, and missing/unexpected
metrics or samples. Sample order is immaterial. Expected IDs and bounds are
consumer-owned in the runners, not supplied by the producer.

Values are exported as finite binary64 numbers before display rounding. This
preserves more precision than the printed log, but does not turn a numerical
residual into an exact symbolic result or an interval certificate. The writer
atomically replaces the destination only when the producer explicitly calls
`finish()` at the end; an exception or early exit does not complete the run.
Exit status is checked even when a result file exists.

## Extending or changing a check

1. Record the already computed numerical value with
   `ResultWriter.add(metric_id, sample_id, value)` at the corresponding check.
2. Call `finish()` after all required work succeeds.
3. Independently update the expected sample IDs and numerical bounds in the
   consuming runner. Keep IDs stable across prose edits; changes to the
   experiment's parameters or sample set require explicit review.
4. Add failure cases and compare the scientific outputs with the previous
   release. Update provenance through its generator when hashed sources change.

No extra third-party dependency is required. `--format legacy` in both runners
retains the previous stdout parser for migration comparisons; it is not an
automatic fallback. The earlier parser tests explicitly exercise that mode.

This migration covers the 12 Python regression programs and four Sage programs
(six Sage checks). Other consumers, including the richer adiabatic provenance
generator and Wolfram checks, retain their existing interfaces. Their migration
requires a separate contract for their additional fields.
