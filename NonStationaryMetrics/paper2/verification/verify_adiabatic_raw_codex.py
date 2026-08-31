#!/usr/bin/env python3
"""Independent refit of the archived adiabatic convergence pairs.

This script does not import or execute the generators that produced the physical
trajectories.  It treats their NPZ output as experimental data and tests the
asymptotic claim by several fits that are different from the provenance script.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np


RAW = Path(__file__).resolve().parents[2] / "ThakurtaMetric" / "adiabatic_convergence_raw.npz"


def slope(eps: np.ndarray, residual: np.ndarray) -> float:
    x = np.log(np.asarray(eps, dtype=float))
    y = np.log(np.asarray(residual, dtype=float))
    x0 = x - x.mean()
    return float(np.dot(x0, y - y.mean()) / np.dot(x0, x0))


def main() -> int:
    data = np.load(RAW, allow_pickle=True)
    indices = sorted(int(k.split("_")[-1]) for k in data.files if k.startswith("eps_"))
    checks: list[tuple[str, bool]] = []

    print(f"dataset: {RAW}")
    for i in indices:
        eps = np.asarray(data[f"eps_{i}"], dtype=float)
        lead = np.asarray(data[f"res_leading_{i}"], dtype=float)
        exact = np.asarray(data[f"res_exact_{i}"], dtype=float)
        label = str(data[f"label_{i}"])

        full_lead = slope(eps, lead)
        full_exact = slope(eps, exact)
        shrinking = [slope(eps[:n], exact[:n]) for n in range(len(eps), 2, -1)]
        pairwise = np.log(exact[1:] / exact[:-1]) / np.log(eps[1:] / eps[:-1])

        # If R=c2 eps^2+c3 eps^3+..., then R/eps^2 approaches c2.
        scaled = exact / eps**2
        coeff = np.linalg.lstsq(
            np.column_stack((eps**2, eps**3)), exact, rcond=None
        )[0]
        model = coeff[0] * eps**2 + coeff[1] * eps**3
        model_rel = float(np.linalg.norm(exact - model) / np.linalg.norm(exact))

        print(f"\n{label}")
        print(f"  leading slope                 {full_lead:.6f}")
        print(f"  complete slope                {full_exact:.6f}")
        print("  shrinking-window slopes       " + " -> ".join(f"{v:.6f}" for v in shrinking))
        print("  adjacent Richardson exponents " + ", ".join(f"{v:.6f}" for v in pairwise))
        print("  R/eps^2                       " + ", ".join(f"{v:.8g}" for v in scaled))
        print(f"  c2 eps^2+c3 eps^3 rel. error  {model_rel:.3e}")

        checks.extend(
            [
                (f"{i}: leading residual is first order", 0.90 < full_lead < 1.15),
                (f"{i}: complete residual has exponent near two", 1.90 < full_exact < 2.25),
                (
                    f"{i}: narrowing the window moves the exponent toward two",
                    abs(shrinking[-1] - 2.0) < abs(shrinking[0] - 2.0),
                ),
                (f"{i}: quadratic+cubic asymptotic model fits", model_rel < 2.0e-3),
            ]
        )

    for label, ok in checks:
        print(f"{'PASS' if ok else 'FAIL'}  {label}")
    passed = sum(ok for _, ok in checks)
    print(f"\nsummary: {passed}/{len(checks)} checks passed")
    return 0 if passed == len(checks) else 1


if __name__ == "__main__":
    raise SystemExit(main())
