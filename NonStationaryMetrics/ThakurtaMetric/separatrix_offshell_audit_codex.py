# -*- coding: utf-8 -*-
"""
Off-shell separatrix audit after the generic t/eta and tau closures.

The script separates two mathematically different degenerations:

1. tau ergospheric algebraic degeneration J=a/E:
   one factor r-2M comes from T and one from Q3.  The kernel numerator carries
   the same factor, so the off-shell kernel has only a simple pole.  The curve
   reduces to genus one and the wrap closes in elliptic Abelian integrals and
   elliptic dilogarithms.

2. physical t/eta moving-double-root separatrix:
   Q2=(r-r_d)**2*q2 with a nonzero kernel numerator.  Hence
   K_off~kappa_3/(r-r_d)**3 while the accumulated radial action tends to a
   finite nonzero I_d.  The outer first-order primitive has a double pole.
   A smooth translation r_d -> r_d+epsilon*delta_r_d changes the logarithmic
   frozen shape only by a simple pole, and therefore cannot by itself cancel
   this double pole.  A uniform inner problem on |r-r_d|=O(sqrt(epsilon)) is
   required unless the canonical boundary condition produces an additional
   genuine cubic-pole term.  The latter is now constructed for the
   terminally anchored separatrix BVP in
   separatrix_canonical_cubic_codex.py: re-anchoring the radial action at
   r_d adds +I_d*K_off and cancels the cubic coefficient exactly.  It is not
   available in the distinct fixed-launch IVP.

This is a classification/obstruction result; it deliberately does not call
the tau r=2M degeneration a physical external separatrix.
"""

from __future__ import annotations

import mpmath as mp
import sympy as sp


r, M, a, E, J = sp.symbols("r M a E J", nonzero=True)
D = sp.expand((E**2 - 1) * r + 2 * M)
Delta = sp.expand(r**2 - 2 * M * r + a**2)


def tau_symbolic_data() -> dict[str, sp.Expr]:
    """Exact genus-two to genus-one degeneration at J=a/E."""

    q3 = sp.expand(r * Delta - J**2 * D)
    s_tau = sp.expand(r * (r - 2 * M) * D * q3)
    kernel_numerator = sp.expand(r**2 * D * (E * J * r - 2 * M * a))
    j_degenerate = a / E

    q3_separatrix = sp.factor(q3.subs(J, j_degenerate))
    expected_q3 = sp.factor(
        (r - 2 * M) * (E**2 * r**2 + a**2) / E**2
    )
    if sp.simplify(q3_separatrix - expected_q3) != 0:
        raise AssertionError("Tau J=a/E factorization failed")

    s_separatrix = sp.factor(s_tau.subs(J, j_degenerate))
    expected_s = sp.factor(
        r
        * (r - 2 * M) ** 2
        * D
        * (E**2 * r**2 + a**2)
        / E**2
    )
    if sp.simplify(s_separatrix - expected_s) != 0:
        raise AssertionError("Tau separatrix curve factorization failed")

    reduced_kernel_a = sp.factor(
        sp.cancel(
            kernel_numerator.subs(J, j_degenerate) / q3_separatrix
        )
    )
    expected_kernel_a = sp.factor(
        E**2 * a * r**2 * D / (E**2 * r**2 + a**2)
    )
    if sp.simplify(reduced_kernel_a - expected_kernel_a) != 0:
        raise AssertionError("Tau separatrix kernel cancellation failed")

    # y=(r-2M)Y, Y**2=Q4.  Express the radial-action differential as
    # [polynomial + Delta-pole]/Y.
    q4 = sp.factor(r * D * (E**2 * r**2 + a**2) / E**2)
    action_numerator = sp.expand(
        r * (r - 2 * M) * (E**2 * r**2 + a**2) / E**2
    )
    quotient, remainder = sp.div(
        sp.Poly(action_numerator, r), sp.Poly(Delta, r)
    )
    if sp.expand(
        action_numerator
        - Delta * quotient.as_expr()
        - remainder.as_expr()
    ) != 0:
        raise AssertionError("Tau elliptic action division failed")

    return {
        "Q3_separatrix": q3_separatrix,
        "S_separatrix": s_separatrix,
        "Q4": q4,
        "kernel_A_reduced": reduced_kernel_a,
        "action_polynomial": sp.factor(quotient.as_expr()),
        "action_remainder": sp.factor(remainder.as_expr()),
    }


def tau_numeric_classification() -> list[tuple[float, float, str]]:
    """Classify genuine Q3-double-root degenerations for the validation case."""

    substitutions = {
        M: sp.Integer(1),
        a: sp.Rational(9, 10),
        E: sp.Rational(7, 5),
    }
    q3 = sp.expand((r * Delta - J**2 * D).subs(substitutions))
    discriminant = sp.resultant(q3, sp.diff(q3, r), r)
    j_candidates = sorted(
        float(sp.re(root))
        for root in sp.nroots(discriminant, maxsteps=200)
        if abs(float(sp.im(root))) < 1.0e-10 and float(sp.re(root)) > 0
    )
    horizon = 1.0 + (1.0 - 0.9**2) ** 0.5
    classifications = []
    for j_value in j_candidates:
        roots = [
            complex(root)
            for root in sp.nroots(q3.subs(J, j_value), maxsteps=200)
        ]
        closest_pair = min(
            (
                (first, second)
                for index, first in enumerate(roots)
                for second in roots[index + 1 :]
            ),
            key=lambda pair: abs(pair[0] - pair[1]),
        )
        r_double = float((closest_pair[0].real + closest_pair[1].real) / 2)
        if r_double > horizon:
            label = "external"
        elif r_double > 0:
            label = "inside seed null surface"
        else:
            label = "negative-radius algebraic degeneration"
        classifications.append((j_value, r_double, label))
    return classifications


def t_physical_separatrix_data() -> dict[str, mp.mpf]:
    """High-precision Laurent data for E=1.2,a=0.9,M=1."""

    mp.mp.dps = 70
    mass = mp.mpf(1)
    spin = mp.mpf("0.9")
    energy = mp.mpf("1.2")

    def d_function(radius):
        return (energy**2 - 1) * radius + 2 * mass

    def delta_function(radius):
        return radius**2 - 2 * mass * radius + spin**2

    def q2(radius, angular):
        return (
            2 * energy**2 * angular**2 * mass * radius
            - energy**2 * angular**2 * radius**2
            - 4 * energy**2 * angular * mass * spin * radius
            + 2 * energy**2 * mass * spin**2 * radius
            + energy**2 * spin**2 * radius**2
            + energy**2 * radius**4
            + 4 * angular**2 * mass**2
            - 4 * angular**2 * mass * radius
            + angular**2 * radius**2
            - 8 * angular * mass**2 * spin
            + 4 * angular * mass * spin * radius
            + 4 * mass**2 * spin**2
        )

    r_double, j_critical = mp.findroot(
        lambda radius, angular: (
            q2(radius, angular),
            mp.diff(lambda value: q2(value, angular), radius),
        ),
        (mp.mpf("1.5123"), mp.mpf("2.93635")),
    )
    q_double = mp.diff(
        lambda value: q2(value, j_critical), r_double, 2
    ) / 2
    t_value = r_double * d_function(r_double)
    kappa_three = (
        energy**2
        * j_critical
        * r_double**4
        * d_function(r_double)
        / (q_double * mp.sqrt(t_value * q_double))
    )

    def spectral_curve(radius):
        return radius * d_function(radius) * q2(radius, j_critical)

    def physical_momentum(radius):
        value = spectral_curve(radius)
        if value < 0 and abs(value) < mp.mpf("1e-45"):
            value = mp.mpf(0)
        return -mp.sqrt(value) / (
            delta_function(radius) * d_function(radius)
        )

    r0 = mp.mpf(20)
    action_at_double_root = mp.quad(
        lambda value: physical_momentum(value),
        [r0, 10, 5, 2, r_double],
    )

    # Derivatives of the double-root family Q2=Q2_r=0.
    rs, js, es = sp.symbols("rs js es")
    q_symbolic = (
        2 * es**2 * js**2 * rs
        - es**2 * js**2 * rs**2
        - 4 * es**2 * js * sp.Rational(9, 10) * rs
        + 2 * es**2 * sp.Rational(81, 100) * rs
        + es**2 * sp.Rational(81, 100) * rs**2
        + es**2 * rs**4
        + 4 * js**2
        - 4 * js**2 * rs
        + js**2 * rs**2
        - 8 * js * sp.Rational(9, 10)
        + 4 * js * sp.Rational(9, 10) * rs
        + 4 * sp.Rational(81, 100)
    )
    derivative_functions = [
        sp.lambdify((rs, js, es), expression, "mpmath")
        for expression in (
            sp.diff(q_symbolic, es),
            sp.diff(q_symbolic, js),
            sp.diff(sp.diff(q_symbolic, rs), es),
            sp.diff(sp.diff(q_symbolic, rs), js),
            sp.diff(q_symbolic, rs, 2),
        )
    ]
    q_e, q_j, q_re, q_rj, q_rr = [
        function(r_double, j_critical, energy)
        for function in derivative_functions
    ]
    dj_de = -q_e / q_j
    dr_de = -(q_re + q_rj * dj_de) / q_rr

    return {
        "r_double": r_double,
        "j_critical": j_critical,
        "dJc_dE": dj_de,
        "drd_dE": dr_de,
        "kappa_3": kappa_three,
        "I_d": action_at_double_root,
        "integrand_cubic_coefficient": kappa_three
        * action_at_double_root,
        "primitive_double_coefficient": kappa_three
        * action_at_double_root
        / 2,
    }


def main() -> None:
    tau = tau_symbolic_data()
    print("tau J=a/E algebraic degeneration:")
    print("  Q3 =", tau["Q3_separatrix"])
    print("  S  =", tau["S_separatrix"])
    print("  Q4 =", tau["Q4"])
    print("  reduced kernel A =", tau["kernel_A_reduced"])
    print("  elliptic action polynomial =", tau["action_polynomial"])
    print("  elliptic action Delta remainder =", tau["action_remainder"])
    print("  pole order: K_off ~ (r-2M)^-1, wrap ~ log(r-2M)")

    print("\ntau genuine Q3-double-root classification at M=1,a=0.9,E=1.4:")
    for j_value, r_double, label in tau_numeric_classification():
        print(f"  J={j_value:.12f}, r_d={r_double:+.12f}: {label}")
    print("  J=a/E=0.642857... is the separate r=2M factor-collision threshold.")

    t_data = t_physical_separatrix_data()
    print("\nphysical t/eta double-root data at M=1,a=0.9,E=1.2:")
    for key in (
        "j_critical",
        "r_double",
        "dJc_dE",
        "drd_dE",
        "kappa_3",
        "I_d",
        "integrand_cubic_coefficient",
        "primitive_double_coefficient",
    ):
        print(f"  {key} = {mp.nstr(t_data[key], 30)}")

    print("\ntracking audit:")
    print("  frozen shape       ~ c*log|x|, x=r-r_d")
    print("  smooth root shift  ~ -c*delta_r_d/x       (simple pole)")
    print("  dilation wrap      ~ C/x^2                 (double pole)")
    print("  therefore a coordinate/root shift alone cannot cancel the dilation pole.")
    print("  terminal canonical anchoring supplies +I_d*K_off and cancels it;")
    print("  in the fixed-launch IVP the inner scale remains |x|=O(sqrt(epsilon)).")
    print("  see separatrix_canonical_cubic_codex.py for the two-BVP distinction.")


if __name__ == "__main__":
    main()
