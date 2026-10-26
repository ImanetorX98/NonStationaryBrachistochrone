#!/usr/bin/env python3
"""Run every deterministic experimental proof-of-concept generator."""

from __future__ import annotations

from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent
SCRIPTS = (
    "active_particle_poc.py",
    "draining_vortex_poc.py",
    "ring_bec_poc.py",
)


def main() -> None:
    for script in SCRIPTS:
        path = HERE / script
        print(f"\n=== {script} ===", flush=True)
        subprocess.run([sys.executable, str(path)], cwd=HERE, check=True)
    print("\n=== validate_pocs.py ===", flush=True)
    subprocess.run([sys.executable, str(HERE / "validate_pocs.py")], cwd=HERE, check=True)
    print(f"\nAll proof-of-concept figures are in {HERE / 'output'}")


if __name__ == "__main__":
    main()
