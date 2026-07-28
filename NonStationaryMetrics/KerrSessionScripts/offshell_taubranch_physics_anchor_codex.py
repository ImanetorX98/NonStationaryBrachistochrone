# -*- coding: utf-8 -*-
"""
Physics anchor for the generic tau-branch off-shell special-function closure.

This script verifies the complete chain on the same physical frozen orbit:

  special-function assembly
      == direct dilation wrap
      == dilation sub-piece extracted from canonical perturbation theory,

and then checks that the full first-order correction containing that sub-piece
has log-log residual slope two against the original-variable non-autonomous
Hamiltonian flow.

Sign convention
---------------
offshell_taubranch_closed_form_codex.py uses the positive algebraic radical
p_plus=+sqrt(S)/(Delta*D).  The ingoing physical shell has
P_r=-sqrt(S)/(Delta*D), while the on-shell kernel is +A/sqrt(S).  Therefore
the physical ingoing dilation wrap is the negative of the positive-radical
assembly returned by validate_case().
"""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import sympy as sp
from scipy.integrate import cumulative_trapezoid as ct
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from tests.test_adiabatic_noreg import _TK_TAU, _tk_slope
from KerrSessionScripts.offshell_taubranch_closed_form_codex import (
    E,
    J,
    M,
    a,
    validate_case,
)


def frozen_tau_orbit():
    """Return lambda,r,P_r and the frozen canonical derivatives."""

    tau = _TK_TAU
    e0, j0, r0 = tau.Ehat, tau.J0, tau.r0
    solution = solve_ivp(
        lambda parameter, state: [
            tau.Hp(state[0], state[1], e0, j0),
            -tau.Hr(state[0], state[1], e0, j0),
            tau.HJ(state[0], state[1], e0, j0),
        ],
        [0.0, 300.0],
        [r0, tau.ingoing_pr(r0, e0, j0), 0.0],
        rtol=1.0e-12,
        atol=1.0e-14,
        max_step=0.005,
        dense_output=True,
        events=tau._event(),
    )
    # The nested source is extracted with two cumulative integrations.  A dense
    # sampling keeps this independent trapezoidal link below the 1e-8 level.
    parameter = np.linspace(0.0, solution.t[-1], 60000)
    radius, momentum, angle = solution.sol(parameter)
    h_pr = tau.Hp(radius, momentum, e0, j0)
    g_pr = tau.Gpr(radius, momentum, e0, j0)
    return parameter, radius, momentum, angle, h_pr, g_pr


def physical_dilation_subpiece(radius_target: float) -> float:
    """Extract -int (G_Pr/H_Pr) int(P_r H_Pr dlambda) dr."""

    parameter, radius, momentum, _, h_pr, g_pr = frozen_tau_orbit()
    source_dilation = ct(momentum * h_pr, parameter, initial=0.0)
    coefficient = ct(-(g_pr / h_pr) * source_dilation, radius, initial=0.0)
    return float(
        interp1d(
            radius,
            coefficient,
            fill_value="extrapolate",
            bounds_error=False,
        )(radius_target)
    )


def residual_data() -> tuple[np.ndarray, np.ndarray, float]:
    """Residuals of the full corrected tau first order vs original flow."""

    tau = _TK_TAU
    phi0, correction = tau.delta_phi(corrected=True)
    radius_eval = np.linspace(8.0, 11.0, 1500)
    eps_values = np.array([1.0e-3, 2.0e-3, 4.0e-3, 8.0e-3])
    residuals = []
    for epsilon in eps_values:
        radius, phi = tau.flow_original(epsilon)
        phi_true = interp1d(
            radius,
            phi,
            fill_value="extrapolate",
            bounds_error=False,
        )(radius_eval)
        residuals.append(
            np.nanmax(
                np.abs(
                    phi_true
                    - (phi0(radius_eval) + epsilon * correction(radius_eval))
                )
            )
        )
    residuals = np.asarray(residuals)
    slope = float(
        np.polyfit(np.log(eps_values), np.log(residuals), 1)[0]
    )
    return eps_values, residuals, slope


def representation_error(epsilon: float = 0.02) -> float:
    """Original canonical variables vs normalized variables with dilation."""

    tau = _TK_TAU
    radius_original, phi_original = tau.flow_original(epsilon)
    radius_normalized, phi_normalized = tau.flow_normalized(
        epsilon, damping=True
    )
    radius_eval = np.linspace(8.0, 11.0, 1500)
    original = interp1d(
        radius_original,
        phi_original,
        fill_value="extrapolate",
        bounds_error=False,
    )(radius_eval)
    normalized = interp1d(
        radius_normalized,
        phi_normalized,
        fill_value="extrapolate",
        bounds_error=False,
    )(radius_eval)
    return float(np.nanmax(np.abs(original - normalized)))


def main() -> None:
    substitutions = {
        M: sp.Integer(1),
        a: sp.Rational(9, 10),
        E: sp.Rational(7, 5),
        J: sp.Rational(5, 2),
    }
    special = validate_case(substitutions, 12.0, 6.0)
    physical_closed = -special["assembled"]
    physical_machine = physical_dilation_subpiece(6.0)

    print("tau dilation wrap on the same ingoing frozen orbit, r0=12 -> r=6:")
    print(f"  special-function closed = {physical_closed:+.12e}")
    print(f"  canonical-PT sub-piece  = {physical_machine:+.12e}")
    print(f"  link error              = {abs(physical_closed-physical_machine):.3e}")

    corrected_slope = _tk_slope(_TK_TAU, corrected=True)
    old_slope = _tk_slope(_TK_TAU, corrected=False)
    eps_values, residuals, direct_slope = residual_data()
    repr_error = representation_error()

    print("\nfull tau first order against the original-variable physical flow:")
    for epsilon, residual in zip(eps_values, residuals):
        print(f"  eps={epsilon:.3e} residual={residual:.12e}")
    print(f"  direct log-log slope    = {direct_slope:.6f}")
    print(f"  oracle corrected slope  = {corrected_slope:.6f}")
    print(f"  old no-dilation slope   = {old_slope:.6f}")
    print(f"  original/normalized err = {repr_error:.3e}")

    if abs(physical_closed - physical_machine) > 2.0e-8:
        raise SystemExit("tau closed-form/physics link failed")
    if direct_slope < 1.8 or corrected_slope < 1.8:
        raise SystemExit("tau physical slope is not second order")
    if old_slope > 1.3:
        raise SystemExit("tau anti-regression control did not expose the old bug")
    if repr_error > 1.0e-10:
        raise SystemExit("tau canonical representations disagree")

    print("\nTau chain passed: closed form = PT sub-piece; full correction = O(eps^2).")


if __name__ == "__main__":
    main()
