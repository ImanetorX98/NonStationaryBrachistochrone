# -*- coding: utf-8 -*-
"""
Generic closed form for the Thakurta--Kerr t/eta-branch off-shell dilation wrap.

The result is a canonical genus-two iterated-Abelian closed form.  It is not an
elementary primitive.  All coefficients are exact functions of (M, a, E, J);
none of those four parameters is frozen.

For a consistently chosen sheet y**2 = S and base point r0, define

    U_k(r)       = int[r0,r] x**k dx/y(x),
    Pi_q(r)      = int[r0,r] dx/((x-q)y(x)),
    W_jk(r)      = int[r0,r] (U_k dU_j - U_j dU_k),
    D_jq(r)      = int[r0,r] U_j(x) dx/((x-q)y(x)).

The radial action is

    I = sum(k=0..3) a_k U_k + sum(q=r_+,r_-) rho_q Pi_q,

and the off-shell contribution is

    Phi_off =
        -B I + E_P
        -1/2 sum(j=0..4,k=0..3) g_j a_k (U_j U_k + W_jk)
        -Phi_G Pi
        +sum(q=r_+,r_-) rho_q sum(j=0..4) g_j D_jq,

where B=P*y/Q2, Phi_G=sum g_j U_j, Pi=sum rho_q Pi_q and
E_P=int P(r)r/Delta(r) dr.  E_P is explicitly polynomial plus logarithms.

The Hermite polynomial P is the unique degree <= 3 solution of

    P*T*Q2' = -2*N  (mod Q2),          T=r*((E**2-1)r+2M),

and can equivalently be obtained from the exact 4x4 linear system returned by
hermite_linear_system().  This avoids a very large, fragile expanded rational
expression while remaining an explicit algebraic formula (Cramer's rule gives
each coefficient).

Run this file to perform:
  * parameter-general symbolic identities;
  * exact Hermite reductions at two independent rational parameter sets;
  * numerical comparison of the closed assembly with the defining nested
    integral.
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
T = sp.expand(r * D)
S = sp.expand(T * Q2)

# A/y is the on-shell outer kernel.  Its polynomial part is J*D.
A = sp.expand(E**2 * J * r**4 * D)
N = sp.expand(sp.rem(A, Q2, r))

# S/(Delta*D) = sum a_k r^k + C/Delta.  The factor D in the
# polynomial-division remainder cancels exactly, so there is no pole at D=0.
ACTION_COEFFICIENTS = (
    2 * M * (4 * E**2 * M**2 - J**2 + 2 * (1 - E**2) * J * a),
    4 * E**2 * M**2 - (E**2 - 1) * J**2,
    2 * E**2 * M,
    E**2,
)
C_DELTA = sp.expand(
    (
        E**2 * J**2 * a**2
        - 8 * E**2 * J * M**2 * a
        + 16 * E**2 * M**4
        - 4 * E**2 * M**2 * a**2
        - J**2 * a**2
        + 4 * M**2 * a**2
    )
    * r
    + 4 * E**2 * J * M * a**3
    - 8 * E**2 * M**3 * a**2
    + 2 * J**2 * M * a**2
    - 4 * J * M * a**3
)


def _coeff_vector(poly: sp.Expr, degree: int) -> sp.Matrix:
    expanded = sp.expand(poly)
    return sp.Matrix([expanded.coeff(r, k) for k in range(degree + 1)])


def hermite_linear_system() -> tuple[sp.Matrix, sp.Matrix]:
    """Return H,b such that H*(P0,P1,P2,P3)^T=b.

    H is the matrix of multiplication by T*Q2' in the quotient ring
    K[r]/(Q2), K=Q(M,a,E,J).  Thus this is an exact, finite specification
    of every coefficient of P for generic nonsingular parameters.
    """

    q2_prime = sp.diff(Q2, r)
    columns = []
    for power in range(4):
        reduced = sp.rem(sp.expand(r**power * T * q2_prime), Q2, r)
        columns.append(_coeff_vector(reduced, 3))
    h_matrix = sp.Matrix.hstack(*columns)
    b_vector = -2 * _coeff_vector(N, 3)
    return h_matrix, b_vector


@dataclass(frozen=True)
class HermiteData:
    substitutions: dict[sp.Symbol, sp.Expr]
    p_polynomial: sp.Expr
    remainder_polynomial: sp.Expr
    g_polynomial: sp.Expr
    g_coefficients: tuple[sp.Expr, ...]


def hermite_data_exact(substitutions: dict[sp.Symbol, sp.Expr]) -> HermiteData:
    """Compute P, Rem and g_j exactly after an exact parameter substitution."""

    q2 = sp.expand(Q2.subs(substitutions))
    t_poly = sp.expand(T.subs(substitutions))
    n_poly = sp.expand(N.subs(substitutions))
    q2_prime = sp.diff(q2, r)
    t_prime = sp.diff(t_poly, r)

    inverse = sp.invert(
        sp.Poly(sp.expand(q2_prime * t_poly), r),
        sp.Poly(q2, r),
    )
    p_poly = sp.expand(
        sp.rem(sp.expand(-2 * n_poly * inverse.as_expr()), q2, r)
    )

    numerator = sp.expand(
        2 * n_poly
        - 2 * sp.diff(p_poly, r) * q2 * t_poly
        - p_poly * q2 * t_prime
        + p_poly * q2_prime * t_poly
    )
    quotient, division_remainder = sp.div(
        sp.Poly(numerator, r), sp.Poly(2 * q2, r)
    )
    if division_remainder.as_expr() != 0:
        raise AssertionError("Hermite reduction is not an exact polynomial division")

    rem_poly = sp.expand(quotient.as_expr())
    g_poly = sp.expand((J * D).subs(substitutions) + rem_poly)
    coefficients = tuple(
        sp.Poly(g_poly, r).coeff_monomial(r**power) for power in range(5)
    )

    hermite_numerator = sp.expand(
        2 * n_poly
        - 2 * sp.diff(p_poly, r) * q2 * t_poly
        - p_poly * q2 * t_prime
        + p_poly * q2_prime * t_poly
        - 2 * rem_poly * q2
    )
    if hermite_numerator != 0:
        raise AssertionError("Hermite differential identity failed")

    return HermiteData(
        substitutions=dict(substitutions),
        p_polynomial=p_poly,
        remainder_polynomial=rem_poly,
        g_polynomial=g_poly,
        g_coefficients=coefficients,
    )


def horizon_data() -> tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr]:
    """Return r_+, r_-, rho_+, rho_- for the simple-horizon case."""

    root_gap = sp.sqrt(M**2 - a**2)
    r_plus = M + root_gap
    r_minus = M - root_gap
    rho_plus = sp.simplify(C_DELTA.subs(r, r_plus) / (r_plus - r_minus))
    rho_minus = sp.simplify(C_DELTA.subs(r, r_minus) / (r_minus - r_plus))
    return r_plus, r_minus, rho_plus, rho_minus


def elementary_primitive_data(
    p_polynomial: sp.Expr, substitutions: dict[sp.Symbol, sp.Expr]
) -> tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr, sp.Expr]:
    """Return q(r), r_+, r_-, c_+, c_- for

        P(r)r/Delta(r) = q(r) + c_+/(r-r_+) + c_-/(r-r_-).
    """

    delta = sp.expand(Delta.subs(substitutions))
    quotient, remainder = sp.div(
        sp.Poly(sp.expand(r * p_polynomial), r), sp.Poly(delta, r)
    )
    if sp.degree(remainder.as_expr(), r) > 1:
        raise AssertionError("Unexpected elementary-remainder degree")

    root_gap = sp.sqrt((M**2 - a**2).subs(substitutions))
    mass = M.subs(substitutions)
    r_plus = sp.simplify(mass + root_gap)
    r_minus = sp.simplify(mass - root_gap)
    c_plus = sp.simplify(
        (r_plus * p_polynomial.subs(r, r_plus)) / (r_plus - r_minus)
    )
    c_minus = sp.simplify(
        (r_minus * p_polynomial.subs(r, r_minus)) / (r_minus - r_plus)
    )
    return quotient.as_expr(), r_plus, r_minus, c_plus, c_minus


def symbolic_checks() -> None:
    """Cheap parameter-general identities; no parameter is specialized."""

    action_polynomial = sum(
        ACTION_COEFFICIENTS[k] * r**k for k in range(4)
    )
    action_identity = sp.expand(
        S - Delta * D * action_polynomial - D * C_DELTA
    )
    if action_identity != 0:
        raise AssertionError("Generic radial-action reduction failed")

    kernel_quotient, kernel_remainder = sp.div(
        sp.Poly(A, r), sp.Poly(Q2, r)
    )
    if sp.expand(kernel_quotient.as_expr() - J * D) != 0:
        raise AssertionError("Generic kernel polynomial quotient failed")
    if sp.expand(kernel_remainder.as_expr() - N) != 0:
        raise AssertionError("Generic kernel remainder failed")

    h_matrix, b_vector = hermite_linear_system()
    p_symbols = sp.Matrix(sp.symbols("P0:4"))
    formal_p = sum(p_symbols[k] * r**k for k in range(4))
    congruence = sp.rem(
        sp.expand(formal_p * T * sp.diff(Q2, r) + 2 * N), Q2, r
    )
    congruence_vector = _coeff_vector(congruence, 3)
    if sp.expand(congruence_vector - (h_matrix * p_symbols - b_vector)) != sp.zeros(4, 1):
        raise AssertionError("Hermite 4x4 system does not encode the congruence")


def _float(expr: sp.Expr) -> float:
    value = complex(sp.N(expr, 17))
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
        limit=220,
    )[0]


def validate_case(
    substitutions: dict[sp.Symbol, sp.Expr],
    r0: float,
    rf: float,
) -> dict[str, float]:
    """Compare every block with independent real-axis quadratures."""

    data = hermite_data_exact(substitutions)
    q2_expr = sp.expand(Q2.subs(substitutions))
    s_expr = sp.expand(S.subs(substitutions))
    d_expr = sp.expand(D.subs(substitutions))
    delta_expr = sp.expand(Delta.subs(substitutions))
    a_expr = sp.expand(A.subs(substitutions))
    c_expr = sp.expand(C_DELTA.subs(substitutions))

    q2_n = sp.lambdify(r, q2_expr, "numpy")
    s_n = sp.lambdify(r, s_expr, "numpy")
    d_n = sp.lambdify(r, d_expr, "numpy")
    delta_n = sp.lambdify(r, delta_expr, "numpy")
    a_n = sp.lambdify(r, a_expr, "numpy")
    c_n = sp.lambdify(r, c_expr, "numpy")
    p_n = sp.lambdify(r, data.p_polynomial, "numpy")
    g_values = [_float(value) for value in data.g_coefficients]
    action_values = [
        _float(value.subs(substitutions)) for value in ACTION_COEFFICIENTS
    ]

    def y(x: float) -> float:
        value = float(s_n(x))
        if value <= 0:
            raise ValueError(f"S({x})={value} is not positive on the test arc")
        return math.sqrt(value)

    @lru_cache(maxsize=32768)
    def u_cached(x: float, power: int) -> float:
        return _quad(lambda t: t**power / y(t), r0, x)

    def u(x: float, power: int) -> float:
        return u_cached(float(x), power)

    mass = _float(M.subs(substitutions))
    spin = _float(a.subs(substitutions))
    root_gap = math.sqrt(mass**2 - spin**2)
    roots = (mass + root_gap, mass - root_gap)

    c_at_roots = [float(c_n(root)) for root in roots]
    residues = (
        c_at_roots[0] / (roots[0] - roots[1]),
        c_at_roots[1] / (roots[1] - roots[0]),
    )

    def pi_q(x: float, root: float) -> float:
        return _quad(lambda t: 1.0 / ((t - root) * y(t)), r0, x)

    def pi_total(x: float) -> float:
        return sum(
            residues[index] * pi_q(x, roots[index]) for index in range(2)
        )

    def action_basis(x: float) -> float:
        return sum(action_values[k] * u(x, k) for k in range(4)) + pi_total(x)

    def action_direct(x: float) -> float:
        return _quad(lambda t: y(t) / (delta_n(t) * d_n(t)), r0, x)

    def kernel(x: float) -> float:
        return float(a_n(x)) / (float(q2_n(x)) * y(x))

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
                * g_values[j]
                * action_values[k]
                * (u(rf, j) * u(rf, k) + w(j, k, rf))
            )

    phi_g = sum(g_values[j] * u(rf, j) for j in range(5))

    def d_letter(j: int, root: float, x: float) -> float:
        return _quad(
            lambda t: u(t, j) / ((t - root) * y(t)),
            r0,
            x,
        )

    block_b = -phi_g * pi_total(rf)
    for index, root in enumerate(roots):
        block_b += residues[index] * sum(
            g_values[j] * d_letter(j, root, rf) for j in range(5)
        )

    quotient, r_plus, r_minus, c_plus, c_minus = elementary_primitive_data(
        data.p_polynomial, substitutions
    )
    quotient_primitive = sp.integrate(quotient, r)
    qp_n = sp.lambdify(r, quotient_primitive, "numpy")
    rp = _float(r_plus)
    rm = _float(r_minus)
    cp = _float(c_plus)
    cm = _float(c_minus)
    elementary_closed = (
        float(qp_n(rf) - qp_n(r0))
        + cp * math.log((rf - rp) / (r0 - rp))
        + cm * math.log((rf - rm) / (r0 - rm))
    )
    elementary_direct = _quad(
        lambda t: float(p_n(t)) * t / float(delta_n(t)), r0, rf
    )

    b_at_rf = float(p_n(rf)) * y(rf) / float(q2_n(rf))
    block_c = -b_at_rf * action_basis(rf) + elementary_closed
    assembled = block_a + block_b + block_c

    # Independent check of the Hermite primitive itself.
    n_n = sp.lambdify(r, N.subs(substitutions), "numpy")
    rem_n = sp.lambdify(r, data.remainder_polynomial, "numpy")
    b_at_r0 = float(p_n(r0)) * y(r0) / float(q2_n(r0))
    kernel_pole_direct = _quad(
        lambda t: float(n_n(t)) / (float(q2_n(t)) * y(t)), r0, rf
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
    print("Generic symbolic identities: OK")
    print("Hermite system shape:", h_matrix.shape)
    print("No D=0 pole: R_action = D*C_Delta exactly")

    cases = (
        (
            {M: sp.Integer(1), a: sp.Rational(9, 10), E: sp.Rational(7, 5), J: sp.Integer(6)},
            11.0,
            6.5,
        ),
        (
            {M: sp.Rational(6, 5), a: sp.Rational(4, 5), E: sp.Rational(3, 2), J: sp.Integer(7)},
            13.0,
            7.0,
        ),
    )
    for index, (substitutions, r0, rf) in enumerate(cases, start=1):
        result = validate_case(substitutions, r0, rf)
        print(f"\ncase {index}: {substitutions}, r0={r0}, r={rf}")
        print(
            "  radial action error  = "
            f"{result['action_error']:.3e}"
        )
        print(
            "  Hermite error        = "
            f"{result['hermite_error']:.3e}"
        )
        print(
            "  elementary error     = "
            f"{result['elementary_error']:.3e}"
        )
        print(
            "  off-shell direct     = "
            f"{result['direct']:+.12e}"
        )
        print(
            "  off-shell closed     = "
            f"{result['assembled']:+.12e}"
        )
        print(
            "  assembly error       = "
            f"{result['assembly_error']:.3e}"
        )

        if max(
            result["action_error"],
            result["hermite_error"],
            result["elementary_error"],
            result["assembly_error"],
        ) > 2.0e-8:
            raise SystemExit("validation tolerance exceeded")

    print("\nAll closed-form identities and numerical validations passed.")


if __name__ == "__main__":
    main()
