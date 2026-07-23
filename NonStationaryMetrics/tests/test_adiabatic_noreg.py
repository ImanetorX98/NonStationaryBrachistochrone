# -*- coding: utf-8 -*-
"""
Anti-regression oracle for the first-order adiabatic correction (main11 referee).

Why this test exists
--------------------
The earlier validation was CIRCULAR: the "true non-autonomous geodesic" was
integrated from the SAME normalized equations used to derive the correction, both
omitting the dilation term -alpha*P_r that accompanies the time-dependent canonical
transformation P_r = p_r/A of Thakurta-Kerr. Two reorganizations of the same
(surrogate) system agree to O(eps^2), giving a false slope-2 pass.

This test builds an INDEPENDENT ground truth by integrating the ORIGINAL Hamiltonian
H = A * Hbar(r, p_r/A, p_phi/A ; Ehat/A) in the UN-transformed canonical variables
(r, p_r), where Hamilton's equations are elementary and unambiguous. It then asserts:

  TK (t and tau branches):
    * original-variable flow == normalized flow WITH -alpha*P_r      (< 1e-10)
    * corrected source  S_D = int (Theta H + P_r H_Pr) dl  ->  slope ~ 2
    * old source        S   = int  Theta H            dl  ->  slope ~ 1  (the bug)

  Vaidya (v branch): the slow datum m is a metric modulus, NOT a conformal
  rescaling of momenta, so there is NO -alpha*P_r term:
    * canonical source  Theta = m d_m                     ->  slope ~ 2
    * spuriously adding  + P_r H_pr                        ->  slope ~ 1  (breaks it)

If a future edit reintroduces the omission (or wrongly adds the term to Vaidya),
the corresponding slope flips and the test fails.

Run:  pytest -q tests/test_adiabatic_noreg.py     (or  python tests/test_adiabatic_noreg.py)
"""
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp, cumulative_trapezoid as ct
from scipy.optimize import brentq
from scipy.interpolate import interp1d

np.seterr(all="ignore")
M, A_SPIN = 1.0, 0.9
_r, _pr, _E, _J = sp.symbols("r pr E J_")


# ----------------------------------------------------------------------------- TK
def _tk_hamiltonian(branch):
    """Normalized (unit-cost) Thakurta-Kerr optical Hamiltonian Hbar(r,P_r,J_eff;E_eff).
    branch='eta' (t-family, cyclic constant -1) or 'tau' (proper-time momentum shift)."""
    f = 1 - 2 * M / _r
    Dl = _r**2 - 2 * M * _r + A_SPIN**2
    b = 2 * M * A_SPIN / _r
    v = 1 - f / _E**2
    Pcap = _r**2 + A_SPIN**2 + 2 * M * A_SPIN**2 / _r
    if branch == "eta":
        Pb = Pcap + b**2 / _E**2
        return _J * b * v / Pb + sp.sqrt(Dl * v / Pb) * sp.sqrt((Dl / _r**2) * _pr**2 + _J**2 / Pb) - 1
    # tau branch: shifted azimuthal momentum ptil = J_eff - b/E_eff (A0=1)
    ptil = _J - b / _E
    Pb = Pcap + b**2 / _E**2
    return ptil * (b * v / Pb) + sp.sqrt(Dl * v / Pb) * sp.sqrt((Dl / _r**2) * _pr**2 + ptil**2 / Pb) - f / _E


def _lam(expr):
    return sp.lambdify((_r, _pr, _E, _J), expr, "numpy")


class _TK:
    def __init__(self, branch, Ehat, J0, r0):
        H = _tk_hamiltonian(branch)
        G = sp.diff(H, _J) / sp.diff(H, _pr)
        self.Ehat, self.J0, self.r0 = Ehat, J0, r0
        self.Hn = _lam(H)
        self.Hp = _lam(sp.diff(H, _pr))
        self.Hr = _lam(sp.diff(H, _r))
        self.HJ = _lam(sp.diff(H, _J))
        self.HE = _lam(sp.diff(H, _E))
        self.Gpr = _lam(sp.diff(G, _pr))
        self.GE = _lam(sp.diff(G, _E))
        self.GJ = _lam(sp.diff(G, _J))

    def ingoing_pr(self, rv, E, Jv):
        pg = np.linspace(-80, 80, 3001)
        Hv = self.Hn(rv, pg, E, Jv)
        roots = [brentq(lambda p: self.Hn(rv, p, E, Jv), pg[i], pg[i + 1])
                 for i in range(len(pg) - 1)
                 if np.isfinite(Hv[i]) and np.isfinite(Hv[i + 1]) and Hv[i] * Hv[i + 1] < 0]
        ing = [p for p in roots if self.Hp(rv, p, E, Jv) < 0]
        return min(ing) if ing else np.nan

    def _event(self):
        ev = lambda l, y: y[1]
        ev.terminal, ev.direction = True, 1
        return ev

    def flow_original(self, eps):
        """Ground truth: original H = A*Hbar in ORIGINAL vars (r, p_r_phys)."""
        E0, J0, r0 = self.Ehat, self.J0, self.r0
        def rhs(s, y):
            r_, pphys, _ = y
            A = np.exp(eps * s); Ee = E0 / A; Je = J0 / A; Pr = pphys / A
            return [self.Hp(r_, Pr, Ee, Je), -A * self.Hr(r_, Pr, Ee, Je), self.HJ(r_, Pr, Ee, Je)]
        ev = lambda s, y: self.Hp(y[0], y[1] / np.exp(eps * s), E0 / np.exp(eps * s), J0 / np.exp(eps * s))
        ev.terminal, ev.direction = True, 1
        so = solve_ivp(rhs, [0, 300], [r0, self.ingoing_pr(r0, E0, J0), 0.0],
                       rtol=1e-12, atol=1e-14, max_step=0.005, dense_output=True, events=ev)
        s = np.linspace(0, so.t[-1], 12000)
        Y = so.sol(s)
        return Y[0], Y[2]

    def flow_normalized(self, eps, damping):
        """Normalized vars; damping=True adds the -alpha*P_r generator term."""
        E0, J0, r0 = self.Ehat, self.J0, self.r0
        def rhs(l, y):
            r_, Pr, _ = y
            A = np.exp(eps * l); E = E0 / A; Jv = J0 / A
            return [self.Hp(r_, Pr, E, Jv),
                    -self.Hr(r_, Pr, E, Jv) - (eps * Pr if damping else 0.0),
                    self.HJ(r_, Pr, E, Jv)]
        so = solve_ivp(rhs, [0, 300], [r0, self.ingoing_pr(r0, E0, J0), 0.0],
                       rtol=1e-12, atol=1e-14, max_step=0.005, dense_output=True, events=self._event())
        l = np.linspace(0, so.t[-1], 12000)
        Y = so.sol(l)
        return Y[0], Y[2]

    def delta_phi(self, corrected):
        """Adiabatic correction dphi(r). corrected=True uses S_D (+P_r H_Pr); False uses old S."""
        E0, J0, r0 = self.Ehat, self.J0, self.r0
        so = solve_ivp(lambda l, y: [self.Hp(y[0], y[1], E0, J0), -self.Hr(y[0], y[1], E0, J0),
                                     self.HJ(y[0], y[1], E0, J0)],
                       [0, 300], [r0, self.ingoing_pr(r0, E0, J0), 0.0],
                       rtol=1e-12, atol=1e-14, max_step=0.005, dense_output=True, events=self._event())
        L = np.linspace(0, so.t[-1], 12000)
        rF, prF, phF = so.sol(L)
        phi0 = interp1d(rF, phF, fill_value="extrapolate", bounds_error=False)
        lam_r = interp1d(rF, L, fill_value="extrapolate", bounds_error=False)
        rg = np.linspace(r0, rF.min() + 0.4, 5000)
        prg = interp1d(rF, prF, fill_value="extrapolate", bounds_error=False)(rg)
        ThetaH = E0 * self.HE(rg, prg, E0, J0) + J0 * self.HJ(rg, prg, E0, J0)
        ThetaG = E0 * self.GE(rg, prg, E0, J0) + J0 * self.GJ(rg, prg, E0, J0)
        PrHpr = prg * self.Hp(rg, prg, E0, J0)
        src = ThetaH + PrHpr if corrected else ThetaH
        S = ct(src * np.gradient(lam_r(rg), rg), rg, initial=0)
        integ = self.Gpr(rg, prg, E0, J0) * (lam_r(rg) * ThetaH - S) / self.Hp(rg, prg, E0, J0) - lam_r(rg) * ThetaG
        xc = interp1d(rg, ct(integ, rg, initial=0), fill_value="extrapolate", bounds_error=False)
        return phi0, xc


# ------------------------------------------------------------------------- Vaidya
def _vaidya():
    r2, p2, m2 = sp.symbols("r pr m", positive=True)
    E2, J2 = sp.symbols("E J_", positive=True)
    f = 1 - 2 * m2 / r2
    w = E2**2 - f
    H = p2 * (f - E2**2) - 1 + sp.sqrt(w) * sp.sqrt(E2**2 * p2**2 + J2**2 / r2**2)
    G = sp.diff(H, J2) / sp.diff(H, p2)
    L = lambda e: sp.lambdify((r2, p2, m2, E2, J2), e, "numpy")
    return dict(Hn=L(H), Hp=L(sp.diff(H, p2)), Hr=L(sp.diff(H, r2)), HJ=L(sp.diff(H, J2)),
                Hm=L(sp.diff(H, m2)), Gpr=L(sp.diff(G, p2)), Gm=L(sp.diff(G, m2)))


def _vaidya_slope(spurious):
    V = _vaidya()
    E0, J0, r0, m0 = 1.4, 6.0, 12.0, 1.0

    def ingoing(rv, mv):
        pg = np.linspace(-80, 80, 4001)
        Hv = V["Hn"](rv, pg, mv, E0, J0)
        roots = [brentq(lambda p: V["Hn"](rv, p, mv, E0, J0), pg[i], pg[i + 1])
                 for i in range(len(pg) - 1)
                 if np.isfinite(Hv[i]) and np.isfinite(Hv[i + 1]) and Hv[i] * Hv[i + 1] < 0]
        ing = [p for p in roots if V["Hp"](rv, p, mv, E0, J0) < 0]
        return min(ing) if ing else np.nan

    ev = lambda l, y: y[1]
    ev.terminal, ev.direction = True, 1

    def flow(eps):  # canonical: p_r genuine (no rescaling), m(v) = m0 e^{eps l}
        def rhs(l, y):
            rv, pv, _ = y; mv = m0 * np.exp(eps * l)
            return [V["Hp"](rv, pv, mv, E0, J0), -V["Hr"](rv, pv, mv, E0, J0), V["HJ"](rv, pv, mv, E0, J0)]
        so = solve_ivp(rhs, [0, 300], [r0, ingoing(r0, m0), 0.0], rtol=1e-12, atol=1e-14,
                       max_step=0.005, dense_output=True, events=ev)
        l = np.linspace(0, so.t[-1], 12000)
        return l, *so.sol(l)

    l0, rF, prF, phF = flow(0.0)
    phi0 = interp1d(rF, phF, fill_value="extrapolate", bounds_error=False)
    lam_r = interp1d(rF, l0, fill_value="extrapolate", bounds_error=False)
    rg = np.linspace(r0, rF.min() + 0.3, 5000)
    prg = interp1d(rF, prF, fill_value="extrapolate", bounds_error=False)(rg)
    ThetaH = m0 * V["Hm"](rg, prg, m0, E0, J0)
    src = ThetaH + (prg * V["Hp"](rg, prg, m0, E0, J0) if spurious else 0.0)
    S = ct(src * np.gradient(lam_r(rg), rg), rg, initial=0)
    integ = V["Gpr"](rg, prg, m0, E0, J0) * (S - lam_r(rg) * ThetaH) / V["Hp"](rg, prg, m0, E0, J0) \
        + lam_r(rg) * m0 * V["Gm"](rg, prg, m0, E0, J0)
    xc = interp1d(rg, ct(integ, rg, initial=0), fill_value="extrapolate", bounds_error=False)
    rc = np.linspace(8.0, 11.0, 1500)
    eps = np.array([2.5e-3, 5e-3, 1e-2, 2e-2])
    res = []
    for e in eps:
        _, rL, _, pL = flow(e)
        p = interp1d(rL, pL, fill_value="extrapolate", bounds_error=False)(rc)
        res.append(np.nanmax(np.abs(p - (phi0(rc) + e * xc(rc)))))
    return float(np.polyfit(np.log(eps), np.log(np.array(res)), 1)[0])


# ------------------------------------------------------------------------- helpers
def _tk_slope(tk, corrected, eps=np.array([1e-3, 2e-3, 4e-3, 8e-3])):
    phi0, xc = tk.delta_phi(corrected)
    rc = np.linspace(8.0, 11.0, 1500)
    res = []
    for e in eps:
        rL, pL = tk.flow_original(e)
        p = interp1d(rL, pL, fill_value="extrapolate", bounds_error=False)(rc)
        res.append(np.nanmax(np.abs(p - (phi0(rc) + e * xc(rc)))))
    return float(np.polyfit(np.log(eps), np.log(np.array(res)), 1)[0])


# TK configs: (branch, Ehat, J0, r0)
_TK_T = _TK("eta", 1.4, 6.0, 12.0)
_TK_TAU = _TK("tau", 1.4, 2.5, 12.0)   # scattering |J| > J_c


# --------------------------------------------------------------------------- tests
def test_tk_original_equals_damped_normalized():
    """The -alpha*P_r term IS the genuine canonical dynamics (not a modeling choice)."""
    for tk in (_TK_T, _TK_TAU):
        rO, phO = tk.flow_original(0.02)
        rD, phD = tk.flow_normalized(0.02, damping=True)
        rc = np.linspace(8.0, 11.0, 1500)
        fO = interp1d(rO, phO, fill_value="extrapolate", bounds_error=False)
        fD = interp1d(rD, phD, fill_value="extrapolate", bounds_error=False)
        assert np.nanmax(np.abs(fO(rc) - fD(rc))) < 1e-10


def test_tk_t_branch_slopes():
    assert _tk_slope(_TK_T, corrected=True) > 1.7      # S_D closes to O(eps^2)
    assert _tk_slope(_TK_T, corrected=False) < 1.3     # old S leaves O(eps) (the bug)


def test_tk_tau_branch_slopes():
    assert _tk_slope(_TK_TAU, corrected=True) > 1.7
    assert _tk_slope(_TK_TAU, corrected=False) < 1.3


def test_vaidya_no_dilation_term():
    """m is a metric modulus: Theta = m d_m is complete; adding P_r H_pr must BREAK it."""
    assert _vaidya_slope(spurious=False) > 1.8         # canonical source closes
    assert _vaidya_slope(spurious=True) < 1.3          # spurious term destroys closure


if __name__ == "__main__":
    print("TK t-branch    : corrected slope=%.3f  old slope=%.3f"
          % (_tk_slope(_TK_T, True), _tk_slope(_TK_T, False)))
    print("TK tau-branch  : corrected slope=%.3f  old slope=%.3f"
          % (_tk_slope(_TK_TAU, True), _tk_slope(_TK_TAU, False)))
    print("Vaidya         : m d_m slope=%.3f      spurious+P_rH_pr slope=%.3f"
          % (_vaidya_slope(False), _vaidya_slope(True)))
    print("ground truth   : original == damped-normalized (see test_tk_original_equals_damped_normalized)")
    print("OK")
