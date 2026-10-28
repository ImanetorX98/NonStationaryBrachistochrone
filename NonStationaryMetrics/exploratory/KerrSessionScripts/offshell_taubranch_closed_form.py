# -*- coding: utf-8 -*-
"""
Generic special-function closed form of the Thakurta--Kerr tau-branch
off-shell dilation wrap.

This file starts from the physical proper-time Hamiltonian

    H_tau = ptilde*phi0 + sqrt(Delta*v/Pbar)
            *sqrt((Delta/r**2)*P_r**2 + ptilde**2/Pbar) - f/E,

not from the older surrogate with final term -1.  On its frozen shell,

    y**2 = S = r(r-2M)D Q3,
    D = (E**2-1)r+2M,
    Q3 = r*Delta-J**2*D,
    P_r = -y/(Delta*D),

and the off-shell kernel is exactly

    (G_Pr/H_Pr)|shell = A/y,
    A = r**2*D*(E*J*r-2*M*a)/Q3.

The complete dilation wrap

    Phi_tau,dil = -int[r0,r] (A/y)(x) I(x) dx,
    I(x) = int[r0,x] P_r(z) dz,

is assembled into genus-two Abelian and iterated-Abelian letters with all
coefficients symbolic in (M,a,E,J).  The Hermite coefficients are specified by
an exact 3x3 quotient-ring system and may be expanded by Cramer's rule.

Run this file for parameter-general symbolic checks and two independent exact-
parameter numerical validations against the defining nested quadrature.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import math

import numpy as np
import sympy as sp
from scipy.integrate import quad


r = sp.symbols("r")
M, a, E, J = sp.symbols("M a E J", nonzero=True)

D = sp.expand((E**2 - 1) * r + 2 * M)
Delta = sp.expand(r**2 - 2 * M * r + a**2)
Q3 = sp.expand(r * Delta - J**2 * D)
T = sp.expand(r * (r - 2 * M) * D)
S = sp.expand(T * Q3)

# Rational numerator before division by Q3:
# (G_Pr/H_Pr)|shell = KERNEL_NUMERATOR/(Q3*sqrt(S)).
KERNEL_NUMERATOR = sp.expand(r**2 * D * (E * J * r - 2 * M * a))
KERNEL_QUOTIENT_POLY, KERNEL_REMAINDER_POLY = sp.div(
    sp.Poly(KERNEL_NUMERATOR, r), sp.Poly(Q3, r)
)
KERNEL_QUOTIENT = sp.expand(KERNEL_QUOTIENT_POLY.as_expr())
KERNEL_REMAINDER = sp.expand(KERNEL_REMAINDER_POLY.as_expr())

# S/(Delta*D) = sum a_k r^k + C_Delta/Delta.
ACTION_COEFFICIENTS = (
    -2 * M * J**2,
    -(E**2 - 1) * J**2,
    -2 * M,
    sp.Integer(1),
)
C_DELTA = sp.expand(J**2 * a**2 * D)


def _coeff_vector(polynomial: sp.Expr, degree: int) -> sp.Matrix:
    polynomial = sp.expand(polynomial)
    return sp.Matrix([polynomial.coeff(r, power) for power in range(degree + 1)])


def hermite_linear_system() -> tuple[sp.Matrix, sp.Matrix]:
    """Return H,b for P*T*Q3' == -2*N modulo Q3.

    For generic parameters P=P0+P1*r+P2*r**2 is the unique solution of
    H*(P0,P1,P2)^T=b.  Thus P_i=det(H_i)/det(H) are explicit rational
    functions of all four physical parameters.
    """

    q3_prime = sp.diff(Q3, r)
    columns = []
    for power in range(3):
        reduced = sp.rem(sp.expand(r**power * T * q3_prime), Q3, r)
        columns.append(_coeff_vector(reduced, 2))
    return sp.Matrix.hstack(*columns), -2 * _coeff_vector(KERNEL_REMAINDER, 2)


@dataclass(frozen=True)
class HermiteData:
    substitutions: dict[sp.Symbol, sp.Expr]
    p_polynomial: sp.Expr
    remainder_polynomial: sp.Expr
    g_polynomial: sp.Expr
    g_coefficients: tuple[sp.Expr, ...]


def hermite_data_exact(substitutions: dict[sp.Symbol, sp.Expr]) -> HermiteData:
    """Compute the Hermite data with exact arithmetic after substitution."""

    q3 = sp.expand(Q3.subs(substitutions))
    t_poly = sp.expand(T.subs(substitutions))
    n_poly = sp.expand(KERNEL_REMAINDER.subs(substitutions))
    q3_prime = sp.diff(q3, r)
    t_prime = sp.diff(t_poly, r)

    inverse = sp.invert(
        sp.Poly(sp.expand(q3_prime * t_poly), r),
        sp.Poly(q3, r),
    )
    p_polynomial = sp.expand(
        sp.rem(sp.expand(-2 * n_poly * inverse.as_expr()), q3, r)
    )
    numerator = sp.expand(
        2 * n_poly
        - 2 * sp.diff(p_polynomial, r) * q3 * t_poly
        - p_polynomial * q3 * t_prime
        + p_polynomial * q3_prime * t_poly
    )
    quotient, remainder = sp.div(
        sp.Poly(numerator, r), sp.Poly(2 * q3, r)
    )
    if remainder.as_expr() != 0:
        raise AssertionError("Non-exact tau Hermite division")

    rem_polynomial = sp.expand(quotient.as_expr())
    g_polynomial = sp.expand(
        KERNEL_QUOTIENT.subs(substitutions) + rem_polynomial
    )
    g_coefficients = tuple(
        sp.Poly(g_polynomial, r).coeff_monomial(r**power)
        for power in range(5)
    )

    identity_numerator = sp.expand(
        2 * n_poly
        - 2 * sp.diff(p_polynomial, r) * q3 * t_poly
        - p_polynomial * q3 * t_prime
        + p_polynomial * q3_prime * t_poly
        - 2 * rem_polynomial * q3
    )
    if identity_numerator != 0:
        raise AssertionError("Tau Hermite differential identity failed")

    return HermiteData(
        substitutions=dict(substitutions),
        p_polynomial=p_polynomial,
        remainder_polynomial=rem_polynomial,
        g_polynomial=g_polynomial,
        g_coefficients=g_coefficients,
    )


def horizon_data() -> tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr]:
    """Simple-pole data of the radial action for M**2 != a**2."""

    gap = sp.sqrt(M**2 - a**2)
    r_plus = M + gap
    r_minus = M - gap
    rho_plus = sp.simplify(C_DELTA.subs(r, r_plus) / (r_plus - r_minus))
    rho_minus = sp.simplify(C_DELTA.subs(r, r_minus) / (r_minus - r_plus))
    return r_plus, r_minus, rho_plus, rho_minus


def elementary_primitive_data(
    p_polynomial: sp.Expr,
    substitutions: dict[sp.Symbol, sp.Expr],
) -> tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr, sp.Expr]:
    """Return q,r_+,r_-,c_+,c_- for

       P*r*(r-2M)/Delta = q + c_+/(r-r_+) + c_-/(r-r_-).
    """

    delta = sp.expand(Delta.subs(substitutions))
    mass = M.subs(substitutions)
    numerator = sp.expand(p_polynomial * r * (r - 2 * mass))
    quotient, remainder = sp.div(sp.Poly(numerator, r), sp.Poly(delta, r))
    if sp.degree(remainder.as_expr(), r) > 1:
        raise AssertionError("Unexpected tau elementary remainder degree")

    gap = sp.sqrt((M**2 - a**2).subs(substitutions))
    r_plus = sp.simplify(mass + gap)
    r_minus = sp.simplify(mass - gap)
    c_plus = sp.simplify(
        r_plus
        * (r_plus - 2 * mass)
        * p_polynomial.subs(r, r_plus)
        / (r_plus - r_minus)
    )
    c_minus = sp.simplify(
        r_minus
        * (r_minus - 2 * mass)
        * p_polynomial.subs(r, r_minus)
        / (r_minus - r_plus)
    )
    return quotient.as_expr(), r_plus, r_minus, c_plus, c_minus


def symbolic_checks() -> None:
    """Parameter-general identities independent of any numerical orbit."""

    action_polynomial = sum(
        ACTION_COEFFICIENTS[power] * r**power for power in range(4)
    )
    if sp.expand(S - Delta * D * action_polynomial - D * C_DELTA) != 0:
        raise AssertionError("Generic tau action identity failed")

    quotient, remainder = sp.div(
        sp.Poly(KERNEL_NUMERATOR, r), sp.Poly(Q3, r)
    )
    if sp.expand(quotient.as_expr() - KERNEL_QUOTIENT) != 0:
        raise AssertionError("Generic tau kernel quotient failed")
    if sp.expand(remainder.as_expr() - KERNEL_REMAINDER) != 0:
        raise AssertionError("Generic tau kernel remainder failed")

    h_matrix, b_vector = hermite_linear_system()
    p_symbols = sp.Matrix(sp.symbols("P0:3"))
    formal_p = sum(p_symbols[power] * r**power for power in range(3))
    congruence = sp.rem(
        sp.expand(formal_p * T * sp.diff(Q3, r) + 2 * KERNEL_REMAINDER),
        Q3,
        r,
    )
    if sp.expand(
        _coeff_vector(congruence, 2) - (h_matrix * p_symbols - b_vector)
    ) != sp.zeros(3, 1):
        raise AssertionError("Tau Hermite matrix does not encode the congruence")


def _float(expression: sp.Expr) -> float:
    value = complex(sp.N(expression, 17))
    if abs(value.imag) > 1.0e-12:
        raise ValueError(f"Expected a real value, got {value}")
    return float(value.real)


def _quad(function, start: float, stop: float) -> float:
    return quad(
        function,
        start,
        stop,
        epsabs=2.0e-10,
        epsrel=2.0e-10,
        limit=240,
    )[0]


def validate_case(
    substitutions: dict[sp.Symbol, sp.Expr],
    r0: float,
    rf: float,
) -> dict[str, float]:
    """Validate the full special-function assembly against direct quadrature."""

    data = hermite_data_exact(substitutions)
    q3_expr = sp.expand(Q3.subs(substitutions))
    s_expr = sp.expand(S.subs(substitutions))
    d_expr = sp.expand(D.subs(substitutions))
    delta_expr = sp.expand(Delta.subs(substitutions))
    kernel_numerator_expr = sp.expand(KERNEL_NUMERATOR.subs(substitutions))
    c_delta_expr = sp.expand(C_DELTA.subs(substitutions))

    q3_n = sp.lambdify(r, q3_expr, "numpy")
    s_n = sp.lambdify(r, s_expr, "numpy")
    d_n = sp.lambdify(r, d_expr, "numpy")
    delta_n = sp.lambdify(r, delta_expr, "numpy")
    kernel_numerator_n = sp.lambdify(r, kernel_numerator_expr, "numpy")
    c_delta_n = sp.lambdify(r, c_delta_expr, "numpy")
    p_n = sp.lambdify(r, data.p_polynomial, "numpy")
    rem_n = sp.lambdify(r, data.remainder_polynomial, "numpy")

    action_coefficients = [
        _float(value.subs(substitutions)) for value in ACTION_COEFFICIENTS
    ]
    g_coefficients = [_float(value) for value in data.g_coefficients]

    def y(x: float) -> float:
        value = float(s_n(x))
        if value <= 0:
            raise ValueError(f"S_tau({x})={value} is not positive")
        return math.sqrt(value)

    @lru_cache(maxsize=32768)
    def u_cached(x: float, power: int) -> float:
        return _quad(lambda t: t**power / y(t), r0, x)

    def u(x: float, power: int) -> float:
        return u_cached(float(x), power)

    mass = _float(M.subs(substitutions))
    spin = _float(a.subs(substitutions))
    gap = math.sqrt(mass**2 - spin**2)
    roots = (mass + gap, mass - gap)
    residues = (
        float(c_delta_n(roots[0])) / (roots[0] - roots[1]),
        float(c_delta_n(roots[1])) / (roots[1] - roots[0]),
    )

    def pi_q(x: float, root: float) -> float:
        return _quad(lambda t: 1.0 / ((t - root) * y(t)), r0, x)

    def pi_total(x: float) -> float:
        return sum(
            residues[index] * pi_q(x, roots[index]) for index in range(2)
        )

    def action_basis(x: float) -> float:
        return (
            sum(action_coefficients[k] * u(x, k) for k in range(4))
            + pi_total(x)
        )

    def action_direct(x: float) -> float:
        return _quad(lambda t: y(t) / (delta_n(t) * d_n(t)), r0, x)

    def kernel(x: float) -> float:
        return float(kernel_numerator_n(x)) / (float(q3_n(x)) * y(x))

    direct = -_quad(lambda t: kernel(t) * action_direct(t), r0, rf)

    def w(j: int, k: int, x: float) -> float:
        return _quad(
            lambda t: (t**j * u(t, k) - t**k * u(t, j)) / y(t),
            r0,
            x,
        )

    block_a = 0.0
    for j in range(5):
        for k in range(4):
            block_a += (
                -0.5
                * g_coefficients[j]
                * action_coefficients[k]
                * (u(rf, j) * u(rf, k) + w(j, k, rf))
            )

    phi_g = sum(g_coefficients[j] * u(rf, j) for j in range(5))

    def d_letter(j: int, root: float, x: float) -> float:
        return _quad(
            lambda t: u(t, j) / ((t - root) * y(t)),
            r0,
            x,
        )

    block_b = -phi_g * pi_total(rf)
    for index, root in enumerate(roots):
        block_b += residues[index] * sum(
            g_coefficients[j] * d_letter(j, root, rf) for j in range(5)
        )

    quotient, r_plus, r_minus, c_plus, c_minus = elementary_primitive_data(
        data.p_polynomial, substitutions
    )
    quotient_primitive = sp.integrate(quotient, r)
    quotient_primitive_n = sp.lambdify(r, quotient_primitive, "numpy")
    rp = _float(r_plus)
    rm = _float(r_minus)
    cp = _float(c_plus)
    cm = _float(c_minus)
    elementary_closed = (
        float(quotient_primitive_n(rf) - quotient_primitive_n(r0))
        + cp * math.log((rf - rp) / (r0 - rp))
        + cm * math.log((rf - rm) / (r0 - rm))
    )
    elementary_direct = _quad(
        lambda t: float(p_n(t))
        * t
        * (t - 2 * mass)
        / float(delta_n(t)),
        r0,
        rf,
    )

    b_at_rf = float(p_n(rf)) * y(rf) / float(q3_n(rf))
    block_c = -b_at_rf * action_basis(rf) + elementary_closed
    assembled = block_a + block_b + block_c

    n_n = sp.lambdify(r, KERNEL_REMAINDER.subs(substitutions), "numpy")
    b_at_r0 = float(p_n(r0)) * y(r0) / float(q3_n(r0))
    kernel_pole_direct = _quad(
        lambda t: float(n_n(t)) / (float(q3_n(t)) * y(t)), r0, rf
    )
    kernel_pole_closed = (
        b_at_rf
        - b_at_r0
        + _quad(lambda t: float(rem_n(t)) / y(t), r0, rf)
    )

    return {
        "action_error": abs(action_direct(rf) - action_basis(rf)),
        "hermite_error": abs(kernel_pole_direct - kernel_pole_closed),
        "elementary_error": abs(elementary_direct - elementary_closed),
        "direct": direct,
        "assembled": assembled,
        "assembly_error": abs(direct - assembled),
    }


def main() -> None:
    symbolic_checks()
    h_matrix, _ = hermite_linear_system()
    print("Generic tau symbolic identities: OK")
    print("Hermite system shape:", h_matrix.shape)
    print("kernel polynomial quotient =", sp.factor(KERNEL_QUOTIENT))
    print("action coefficients =", tuple(sp.factor(x) for x in ACTION_COEFFICIENTS))
    print("third-kind poles: Delta=0 only")

    cases = (
        (
            {
                M: sp.Integer(1),
                a: sp.Rational(9, 10),
                E: sp.Rational(7, 5),
                J: sp.Rational(5, 2),
            },
            12.0,
            6.0,
        ),
        (
            {
                M: sp.Rational(6, 5),
                a: sp.Rational(4, 5),
                E: sp.Rational(3, 2),
                J: sp.Integer(3),
            },
            14.0,
            8.0,
        ),
    )
    for index, (substitutions, r0, rf) in enumerate(cases, start=1):
        result = validate_case(substitutions, r0, rf)
        print(f"\ncase {index}: {substitutions}, r0={r0}, r={rf}")
        print(f"  radial action error  = {result['action_error']:.3e}")
        print(f"  Hermite error        = {result['hermite_error']:.3e}")
        print(f"  elementary error     = {result['elementary_error']:.3e}")
        print(f"  off-shell direct     = {result['direct']:+.12e}")
        print(f"  off-shell closed     = {result['assembled']:+.12e}")
        print(f"  assembly error       = {result['assembly_error']:.3e}")
        if max(
            result["action_error"],
            result["hermite_error"],
            result["elementary_error"],
            result["assembly_error"],
        ) > 2.0e-8:
            raise SystemExit("tau validation tolerance exceeded")

    print("\nAll tau closed-form identities and numerical validations passed.")


if __name__ == "__main__":
    main()
