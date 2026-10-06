#!/usr/bin/env python3
"""Compile the private Via B project and audit every public theorem's axioms.

Dependencies must already be installed. This script never invokes lake update
or uploads source files. Use --lean-bin for an isolated Lean installation.
The axiom audit does not check whether geometric hypotheses have been proved.
"""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--lean-bin", default=os.environ.get("VIA_B_LEAN_BIN"))
args = parser.parse_args()
env = os.environ.copy()
if args.lean_bin:
    env["PATH"] = str(Path(args.lean_bin).resolve()) + os.pathsep + env.get("PATH", "")


def run(*command):
    result = subprocess.run(command, cwd=ROOT, env=env, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if result.returncode:
        raise SystemExit(result.stdout)
    return result.stdout


version = run("lean", "--version").strip()
if not re.search(r"version 4\.24\.0(?:,|\s|\))", version):
    raise SystemExit(f"Expected Lean 4.24.0, got: {version}")

sources = sorted((ROOT / "ViaB").glob("*.lean"))
theorems = set()
for source in sources:
    content = source.read_text()
    if re.search(r"\b(?:sorry|admit|axiom)\b", content):
        raise SystemExit(f"Unproved declaration or placeholder in {source.name}")
    theorems.update("ViaB." + name for name in re.findall(r"^theorem\s+(\w+)", content, re.M))

audit_source = (ROOT / "AxiomAudit.lean").read_text()
audited = set(re.findall(r"^#print axioms (ViaB\.\w+)$", audit_source, re.M))
if audited != theorems:
    raise SystemExit(f"Axiom audit coverage differs: {sorted(audited ^ theorems)}")

build = run("lake", "build")
axioms = run("lake", "env", "lean", "AxiomAudit.lean")
entries = re.findall(r"'([^']+)' depends on axioms: \[([^\]]*)\]", axioms)
if set(name for name, _ in entries) != theorems or len(entries) != len(theorems):
    raise SystemExit("Incomplete or duplicate transitive axiom output")
allowed = {"propext", "Classical.choice", "Quot.sound"}
for name, dependencies in entries:
    found = {item.strip() for item in dependencies.split(",") if item.strip()}
    if not found <= allowed:
        raise SystemExit(f"Unexpected axiom dependencies for {name}: {sorted(found - allowed)}")

logs = ROOT / "verification"
logs.mkdir(exist_ok=True)
(logs / "build.log").write_text(build)
(logs / "axioms.log").write_text(axioms)
project_files = sources + [ROOT / name for name in
    ("ViaB.lean", "AxiomAudit.lean", "lakefile.toml", "lean-toolchain", "lake-manifest.json", "verify.py")]
report = {
    "checked_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    "lean": version,
    "public_theorems_checked": len(theorems),
    "allowed_standard_axioms": sorted(allowed),
    "full_geometric_maxwell_theorem_formalized": False,
    "unbounded_radius_escape_for_declared_future_ode_formalized": True,
    "bounded_exterior_lingering_for_declared_future_ode_formalized": True,
    "escape_stability_under_uniform_local_field_bound_formalized": True,
    "uniform_field_bound_on_geometric_finite_time_tube_formalized": False,
    "uniform_field_bound_on_regular_finite_time_reference_formalized": True,
    "escape_stability_from_initial_continuity_for_declared_flows_formalized": True,
    "geometric_solution_continuation_formalized": False,
    "sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
               for path in project_files},
}
(logs / "summary.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
print(build.strip())
print(f"PASS: {len(theorems)} public theorems; standard axioms only; no proof placeholders.")
print("Full geometric Maxwell exclusion: NOT YET FORMALIZED. See README.md.")
