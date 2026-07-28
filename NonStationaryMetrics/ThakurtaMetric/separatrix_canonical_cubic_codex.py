# -*- coding: utf-8 -*-
"""
Canonical cubic-pole term for the physical t/eta separatrix.

The generic, launch-anchored dilation wrap is

    dPhi_dil,L/dr = -K_off(r) I_L(r),
    I_L(r) = integral_{r0}^{r} p_0(u) du.

At a double root r_d, K_off ~ kappa_3/(r-r_d)^3 while
I_L(r_d)=I_d is finite and nonzero.  The free-arrival/separatrix-following
canonical problem is instead terminally anchored.  Enforcing the terminal
shell/transversality condition replaces I_L by

    I_T(r) = I_L(r)-I_d = integral_{r_d}^{r} p_0(u) du

and therefore adds the canonical term

    dPhi_can/dr = I_d K_off(r).

Its cubic coefficient is exactly opposite to the launch-anchored dilation
coefficient.  Since p_0=O(r-r_d), I_T=O((r-r_d)^2), so the combined
differential has only a simple pole and its primitive is logarithmic.

This is a boundary-condition statement:

* it applies to the terminally anchored canonical separatrix BVP;
* it must not be added to the fixed-launch IVP, where I_L(r0)=0 is the
  correct integration constant and the O(sqrt(epsilon)) inner layer remains.

Run:
    python3 ThakurtaMetric/separatrix_canonical_cubic_codex.py
"""

from __future__ import annotations

from dataclasses import dataclass

import mpmath as mp
import sympy as sp


r, M, a, E, J, rd, qd = sp.symbols(
    "r M a E J r_d q_d", real=True, nonzero=True
)

D = sp.expand((E**2 - 1) * r + 2 * M)
DELTA = sp.expand(r**2 - 2 * M * r + a**2)
Q2 = sp.expand(
    2 * E**2 * J**2 * M * r
    - E**2 * J**2 * r**2
    - 4 * E**2 * J * M * a * r
    + 2 * E**2 * M * a**2 * r
    + E**2 * a**2 * r**2
    + E**2 * r**4
    + 4 * J**2 * M**2
    - 4 * J**2 * M * r
    + J**2 * r**2
    - 8 * J * M**2 * a
    + 4 * J * M * a * r
    + 4 * M**2 * a**2
)


def symbolic_laurent_data() -> dict[str, sp.Expr]:
    """Return the exact leading coefficients at Q2=(r-r_d)^2 q(r)."""

    d_at_root = sp.expand(D.subs(r, rd))
    delta_at_root = sp.expand(DELTA.subs(r, rd))
    t_at_root = sp.expand(rd * d_at_root)

    # Positive algebraic sheet for y=(r-r_d)Y, Y^2=r D q.
    kappa_three = sp.factor(
        E**2
        * J
        * rd**4
        * d_at_root
        / (qd * sp.sqrt(t_at_root * qd))
    )

    # Physical ingoing momentum p_0=-(r-r_d)Y/(Delta D).
    momentum_linear = sp.factor(
        -sp.sqrt(t_at_root * qd) / (delta_at_root * d_at_root)
    )

    # I_T = integral_{r_d}^r p_0 du
    #     = (momentum_linear/2) (r-r_d)^2 + ...
    terminal_action_quadratic = sp.factor(momentum_linear / 2)

    # -K_off I_T = C_log/(r-r_d)+O(1).
    logarithmic_differential = sp.factor(
        -kappa_three * terminal_action_quadratic
    )
    expected_logarithmic = sp.factor(
        E**2 * J * rd**4 / (2 * qd * delta_at_root)
    )
    if sp.simplify(logarithmic_differential - expected_logarithmic) != 0:
        raise AssertionError("The symbolic simple-pole coefficient did not reduce")

    return {
        "D_d": d_at_root,
        "Delta_d": delta_at_root,
        "q_d": sp.diff(Q2, r, 2).subs(r, rd) / 2,
        "kappa_3": kappa_three,
        "p_1": momentum_linear,
        "I_T_x2": terminal_action_quadratic,
        "C_log": logarithmic_differential,
    }


@dataclass(frozen=True)
class BranchData:
    label: str
    energy: mp.mpf
    angular_momentum: mp.mpf
    double_root: mp.mpf
    d_jc_d_e: mp.mpf
    d_rd_d_e: mp.mpf
    q_double: mp.mpf
    kappa_three: mp.mpf
    momentum_linear: mp.mpf
    action_at_root: mp.mpf
    launch_cubic: mp.mpf
    canonical_cubic: mp.mpf
    logarithmic_differential: mp.mpf


def _q2_numeric(
    radius: mp.mpf,
    angular_momentum: mp.mpf,
    energy: mp.mpf,
    mass: mp.mpf,
    spin: mp.mpf,
) -> mp.mpf:
    return (
        2 * energy**2 * angular_momentum**2 * mass * radius
        - energy**2 * angular_momentum**2 * radius**2
        - 4 * energy**2 * angular_momentum * mass * spin * radius
        + 2 * energy**2 * mass * spin**2 * radius
        + energy**2 * spin**2 * radius**2
        + energy**2 * radius**4
        + 4 * angular_momentum**2 * mass**2
        - 4 * angular_momentum**2 * mass * radius
        + angular_momentum**2 * radius**2
        - 8 * angular_momentum * mass**2 * spin
        + 4 * angular_momentum * mass * spin * radius
        + 4 * mass**2 * spin**2
    )


def _double_root_family_derivatives(
    double_root: mp.mpf,
    angular_momentum: mp.mpf,
    energy: mp.mpf,
    mass: mp.mpf,
    spin: mp.mpf,
) -> tuple[mp.mpf, mp.mpf]:
    """Differentiate Q2=Q2_r=0 along the double-root family."""

    q_e = mp.diff(
        lambda value: _q2_numeric(
            double_root, angular_momentum, value, mass, spin
        ),
        energy,
    )
    q_j = mp.diff(
        lambda value: _q2_numeric(
            double_root, value, energy, mass, spin
        ),
        angular_momentum,
    )
    d_jc_d_e = -q_e / q_j

    q_re = mp.diff(
        lambda value: mp.diff(
            lambda radius: _q2_numeric(
                radius, angular_momentum, value, mass, spin
            ),
            double_root,
        ),
        energy,
    )
    q_rj = mp.diff(
        lambda value: mp.diff(
            lambda radius: _q2_numeric(
                radius, value, energy, mass, spin
            ),
            double_root,
        ),
        angular_momentum,
    )
    q_rr = mp.diff(
        lambda radius: _q2_numeric(
            radius, angular_momentum, energy, mass, spin
        ),
        double_root,
        2,
    )
    d_rd_d_e = -(q_re + q_rj * d_jc_d_e) / q_rr
    return d_jc_d_e, d_rd_d_e


def branch_data(
    label: str,
    root_guess: str,
    angular_guess: str,
    launch_radius: str = "20",
) -> tuple[BranchData, dict[str, object]]:
    """Compute a separatrix branch and return functions used in the limit test."""

    mp.mp.dps = 80
    mass = mp.mpf(1)
    spin = mp.mpf("0.9")
    energy = mp.mpf("1.2")
    r0 = mp.mpf(launch_radius)

    double_root, angular_momentum = mp.findroot(
        lambda radius, angular: (
            _q2_numeric(radius, angular, energy, mass, spin),
            mp.diff(
                lambda value: _q2_numeric(
                    value, angular, energy, mass, spin
                ),
                radius,
            ),
        ),
        (mp.mpf(root_guess), mp.mpf(angular_guess)),
    )

    q_rr = mp.diff(
        lambda radius: _q2_numeric(
            radius, angular_momentum, energy, mass, spin
        ),
        double_root,
        2,
    )
    q_rrr = mp.diff(
        lambda radius: _q2_numeric(
            radius, angular_momentum, energy, mass, spin
        ),
        double_root,
        3,
    )
    q_rrrr = mp.diff(
        lambda radius: _q2_numeric(
            radius, angular_momentum, energy, mass, spin
        ),
        double_root,
        4,
    )
    q_double = q_rr / 2

    def d_function(radius: mp.mpf) -> mp.mpf:
        return (energy**2 - 1) * radius + 2 * mass

    def delta_function(radius: mp.mpf) -> mp.mpf:
        return radius**2 - 2 * mass * radius + spin**2

    # Q2=(r-r_d)^2 q(r); Taylor form avoids catastrophic cancellation.
    def q_factor(radius: mp.mpf) -> mp.mpf:
        x = radius - double_root
        return q_rr / 2 + (q_rrr / 6) * x + (q_rrrr / 24) * x**2

    def elliptic_radical(radius: mp.mpf) -> mp.mpf:
        return mp.sqrt(radius * d_function(radius) * q_factor(radius))

    def momentum(radius: mp.mpf) -> mp.mpf:
        x = radius - double_root
        return (
            -x
            * elliptic_radical(radius)
            / (delta_function(radius) * d_function(radius))
        )

    def kernel(radius: mp.mpf) -> mp.mpf:
        x = radius - double_root
        return (
            energy**2
            * angular_momentum
            * radius**4
            * d_function(radius)
            / (
                x**3
                * q_factor(radius)
                * elliptic_radical(radius)
            )
        )

    action_at_root = mp.quad(
        momentum,
        [r0, 10, 6, 4, 2, double_root]
        if double_root < 2
        else [r0, 10, 6, 4, double_root],
    )
    t_at_root = double_root * d_function(double_root)
    kappa_three = (
        energy**2
        * angular_momentum
        * double_root**4
        * d_function(double_root)
        / (
            q_double
            * mp.sqrt(t_at_root * q_double)
        )
    )
    momentum_linear = (
        -mp.sqrt(t_at_root * q_double)
        / (
            delta_function(double_root)
            * d_function(double_root)
        )
    )
    logarithmic_differential = -kappa_three * momentum_linear / 2
    d_jc_d_e, d_rd_d_e = _double_root_family_derivatives(
        double_root, angular_momentum, energy, mass, spin
    )

    data = BranchData(
        label=label,
        energy=energy,
        angular_momentum=angular_momentum,
        double_root=double_root,
        d_jc_d_e=d_jc_d_e,
        d_rd_d_e=d_rd_d_e,
        q_double=q_double,
        kappa_three=kappa_three,
        momentum_linear=momentum_linear,
        action_at_root=action_at_root,
        launch_cubic=-kappa_three * action_at_root,
        canonical_cubic=kappa_three * action_at_root,
        logarithmic_differential=logarithmic_differential,
    )
    functions: dict[str, object] = {
        "momentum": momentum,
        "kernel": kernel,
        "r0": r0,
    }
    return data, functions


def verify_limits(
    data: BranchData, functions: dict[str, object]
) -> dict[str, mp.mpf]:
    """Numerically check the cubic cancellation and remaining simple pole."""

    momentum = functions["momentum"]
    kernel = functions["kernel"]
    r0 = functions["r0"]
    assert callable(momentum)
    assert callable(kernel)
    assert isinstance(r0, mp.mpf)

    last: dict[str, mp.mpf] = {}
    for exponent in (2, 3, 4, 5):
        x = mp.mpf(10) ** (-exponent)
        radius = data.double_root + x
        terminal_action = mp.quad(
            momentum, [data.double_root, radius]
        )
        launch_action = data.action_at_root + terminal_action
        kernel_value = kernel(radius)

        launch_integrand = -kernel_value * launch_action
        canonical_integrand = kernel_value * data.action_at_root
        combined_integrand = launch_integrand + canonical_integrand

        last = {
            "x": x,
            "x3_kernel": x**3 * kernel_value,
            "terminal_action_over_x2": terminal_action / x**2,
            "x3_launch": x**3 * launch_integrand,
            "x3_canonical": x**3 * canonical_integrand,
            "x3_combined": x**3 * combined_integrand,
            "x_combined": x * combined_integrand,
        }

    tolerances = {
        "kernel": abs(last["x3_kernel"] - data.kappa_three),
        "action": abs(
            last["terminal_action_over_x2"]
            - data.momentum_linear / 2
        ),
        "launch": abs(last["x3_launch"] - data.launch_cubic),
        "canonical": abs(
            last["x3_canonical"] - data.canonical_cubic
        ),
        "cubic_cancel": abs(last["x3_combined"]),
        "simple_pole": abs(
            last["x_combined"] - data.logarithmic_differential
        ),
    }
    last.update({f"error_{key}": value for key, value in tolerances.items()})
    return last


def _print_branch(data: BranchData, limits: dict[str, mp.mpf]) -> None:
    print(f"\n{data.label}:")
    for name, value in (
        ("J_c", data.angular_momentum),
        ("r_d", data.double_root),
        ("dJ_c/dE", data.d_jc_d_e),
        ("dr_d/dE", data.d_rd_d_e),
        ("q_d", data.q_double),
        ("kappa_3", data.kappa_three),
        ("p_1", data.momentum_linear),
        ("I_d", data.action_at_root),
        ("launch cubic coefficient", data.launch_cubic),
        ("canonical cubic coefficient", data.canonical_cubic),
        ("remaining simple-pole coefficient", data.logarithmic_differential),
    ):
        print(f"  {name:34s} = {mp.nstr(value, 28)}")

    print("  limit checks at x=1e-5:")
    for name in (
        "x3_kernel",
        "terminal_action_over_x2",
        "x3_launch",
        "x3_canonical",
        "x3_combined",
        "x_combined",
    ):
        print(f"    {name:31s} = {mp.nstr(limits[name], 20)}")
    print("  absolute errors:")
    for name in (
        "error_kernel",
        "error_action",
        "error_launch",
        "error_canonical",
        "error_cubic_cancel",
        "error_simple_pole",
    ):
        print(f"    {name:31s} = {mp.nstr(limits[name], 8)}")


def main() -> None:
    symbolic = symbolic_laurent_data()
    print("symbolic separatrix Laurent coefficients:")
    print("  q_d    =", symbolic["q_d"])
    print("  kappa3 =", symbolic["kappa_3"])
    print("  p_1    =", symbolic["p_1"])
    print("  C_log  =", symbolic["C_log"])
    print("\ncanonical re-anchoring:")
    print("  dPhi_can/dr = +I_d K_off")
    print("  dPhi_total/dr = -K_off (I_L-I_d)")
    print("  cubic coefficient: -kappa3 I_d + kappa3 I_d = 0")

    branches = (
        ("prograde ergospheric physical branch", "1.5123", "2.93635"),
        ("retrograde external physical branch", "3.5139", "-8.05352"),
    )
    all_ok = True
    for label, root_guess, angular_guess in branches:
        data, functions = branch_data(label, root_guess, angular_guess)
        limits = verify_limits(data, functions)
        _print_branch(data, limits)
        scale = max(
            mp.mpf(1),
            abs(data.kappa_three),
            abs(data.action_at_root),
            abs(data.logarithmic_differential),
        )
        all_ok = all_ok and (
            limits["error_cubic_cancel"] < mp.mpf("1e-8") * scale
            and limits["error_simple_pole"] < mp.mpf("1e-3") * scale
        )

    print("\nPASS" if all_ok else "\nFAILED")
    if not all_ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
