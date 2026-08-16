"""End-to-end numerical audit of eq:vaidya-full (section 4.4 closed form).

LHS : delta_phi/mdot = int_{r0}^{r} (dF/dm) * v dr        [eq:vaidya-deltaphi]
RHS : the fully assembled block expression                [eq:vaidya-full]

with F = phi' = K/sqrt(S), K = J D_E, and the reduction coefficients A^m, c_k
obtained by solving the 11x11 system of Lemma I.4 (nothing fitted).
"""
import mpmath as mp
import sympy as sp

mp.mp.dps = 40

# ---------------------------------------------------------------- parameters
M, EE, JJ = sp.Integer(1), sp.Rational(7, 5), sp.Rational(5, 2)
r = sp.Symbol('r')
m, E, J = sp.symbols('m E J')

D_E = (E**2 - 1) * r + 2 * m
S_s = sp.expand(r * (r - 2 * m) * D_E * (r**2 * (r - 2 * m) - J**2 * D_E))
K_s = J * D_E
N_s = sp.expand(S_s * sp.diff(K_s, m) - K_s * sp.diff(S_s, m) / 2)

sub = {m: M, E: EE, J: JJ}
S_e = sp.expand(S_s.subs(sub))
N_e = sp.expand(N_s.subs(sub))
K_e = sp.expand(K_s.subs(sub))

print("roots of S:", [complex(z) for z in sp.nroots(sp.Poly(S_e, r))])
print("disc_r S  =", sp.discriminant(sp.Poly(S_e, r)),
      " (nonzero => six distinct roots, genus two)")

# ------------------------------------------- Lemma I.4: solve the 11x11 system
a = sp.symbols('a0:6')
c = sp.symbols('c0:5')
A_poly = sum(a[i] * r**i for i in range(6))
ident = sp.expand(2 * sp.diff(A_poly, r) * S_e - A_poly * sp.diff(S_e, r)
                  + 2 * S_e * sum(c[k] * r**k for k in range(5)) - 2 * N_e)
eqs = sp.Poly(ident, r).all_coeffs()
sol = sp.solve(eqs, list(a) + list(c), dict=True)[0]
Am_poly = sp.expand(A_poly.subs(sol))
ck = [sp.nsimplify(sol[c[k]]) for k in range(5)]
print("\nreduction coefficients c_k (exact rationals):")
for k, v in enumerate(ck):
    print(f"   c_{k} = {v}  = {float(v): .10f}")
print("residual of the polynomial identity:",
      sp.simplify(ident.subs(sol)))

# ------------------------------------------------------- numeric callables
f_S = sp.lambdify(r, S_e, 'mpmath')
f_N = sp.lambdify(r, N_e, 'mpmath')
f_K = sp.lambdify(r, K_e, 'mpmath')
f_Am = sp.lambdify(r, Am_poly, 'mpmath')
Efl, Mfl = mp.mpf(7) / 5, mp.mpf(1)
ckf = [mp.mpf(sp.Rational(v).p) / mp.mpf(sp.Rational(v).q) for v in ck]

r0, r1 = mp.mpf(5), mp.mpf(9)
print(f"\nS(r0)={f_S(r0)}, S(r1)={f_S(r1)}  (both must be > 0)")

Q = lambda g, x: mp.quad(g, [r0, x])

# dF/dm = N/S^{3/2};  F = K/S^{1/2}
dFdm = lambda x: f_N(x) / f_S(x)**mp.mpf(1.5)
U = lambda k, x: Q(lambda t: t**k / mp.sqrt(f_S(t)), x)
alpha0 = f_Am(r0) / mp.sqrt(f_S(r0))
A_m = lambda x: f_Am(x) / mp.sqrt(f_S(x)) - alpha0 + sum(
    ckf[k] * U(k, x) for k in range(5))

# consistency: d(A_m)/dr must equal dF/dm  (Lemma I.4 in force)
xt = mp.mpf(7)
print("\ncheck  d(A_m)/dr - dF/dm  at r=7 :",
      mp.nstr(mp.diff(A_m, xt) - dFdm(xt), 8))
print("check  A_m(r0)                    :", mp.nstr(A_m(r0), 8))

# advanced-time clock  v = E U_3 + r + 2m ln(r-2m)
v_of = lambda x: Efl * U(3, x) + x + 2 * Mfl * mp.log(x - 2 * Mfl)

# ------------------------------------------------------------------ LHS
LHS = mp.quad(lambda t: dFdm(t) * v_of(t), [r0, r1])
print("\nLHS  int dF/dm * v dr           =", mp.nstr(LHS, 20))

# ------------------------------------------------------------------ RHS
I_poly = Q(lambda t: f_Am(t) * t**3 / f_S(t), r1)
L_2m = Q(lambda t: f_Am(t) * t / ((t - 2 * Mfl) * mp.sqrt(f_S(t))), r1)
Dk = [Q(lambda t, k=k: U(k, t) / (t - 2 * Mfl), r1) for k in range(5)]
W3k = [Q(lambda t, k=k: (U(3, t) * t**k - U(k, t) * t**3) / mp.sqrt(f_S(t)), r1)
       for k in range(5)]
Uv = [U(k, r1) for k in range(6)]
Am1, v1 = A_m(r1), v_of(r1)

term1 = -Efl * (I_poly - alpha0 * Uv[3]
                + sum(ckf[k] * (Uv[k] * Uv[3] - W3k[k]) / 2 for k in range(5)))
term2 = -(L_2m - alpha0 * (v1 - Efl * Uv[3])
          + sum(ckf[k] * (r1 * Uv[k] - Uv[k + 1] + 2 * Mfl * Dk[k])
                for k in range(5)))
RHS = Am1 * v1 + term1 + term2

print("RHS  eq:vaidya-full assembled   =", mp.nstr(RHS, 20))
print("\nRHS - LHS                       =", mp.nstr(RHS - LHS, 8))
print("relative                        =", mp.nstr(abs((RHS - LHS) / LHS), 8))
