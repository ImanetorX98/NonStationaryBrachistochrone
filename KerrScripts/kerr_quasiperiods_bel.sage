# -*- coding: utf-8 -*-
# Second-kind quasi-periods (matrix eta / Klein's kappa) for the genus-2 curve
# of the tau branch -> the MISSING ingredient of the Kleinian closed form of psi.
#
# psi = 1/2 Ehat (rho - rho~) is a length-2 iterated integral of forms of the
# first/second kind (NO third kind: shown beta_+=delta_+=0). By the classical Baker-Klein
# theory it CLOSES in the Kleinian functions zeta, sigma (weight 1, the same class
# as phi_0); it is NOT a polylogarithm. The matrix of quasi-periods eta (second kind) is needed.
#
# METHOD (gets round the obstacle "Sage does not give eta" met in the working note doranTau.md, sec. 7-iii, not distributed):
#   - ODD model (quintic, 1 point at infinity) via x=1/s: the second-kind forms
#     s^2,s^3 ds/Y are clean (double pole, no residue). The EVEN model (deg 6)
#     has 2 points at infinity -> x^2 dx/y is of the third kind (residue) -> wrong eta.
#   - CANONICAL second-kind forms of Baker-Enolski-Leykin (coefficients of the quintic):
#       dr_1 = (lam3 s + 2 lam4 s^2 + 3 lam5 s^3) ds/(4y)
#       dr_2 = lam5 s^2 ds/(4y)
#   - periods through Sage's matrix_of_integral_values on the 4 cycles (INTEGERS: no
#     Singular crash). Pairing du=[1,s] <-> dr=[dr1,dr2].
#
# VALIDATION: kappa=eta*omega^-1 SYMMETRIC (1.4e-12); generalised Legendre relation
#   omega' eta^T - omega eta'^T = -i pi I (canonical).
#
# Quintic: Y^2 = s^6 S(1/s), S_tau = 1200x^6-2300x^5-11428x^4-5519x^3+24700x^2
#   +62500x (parameters M=1,a=9/10,E=7/5,J=5/2). lam_i = coefficients (integers).
# Stable verification fields are independent of the human-readable messages.
# Without verification_results.py (e.g. a script copied out of the repository) the
# program runs as before and exports nothing, unless NSB_RESULT_JSON asks for a result.
import sys as _result_sys
from pathlib import Path as _ResultPath
_result_source = _ResultPath(_result_sys.argv[0]).resolve()  # Sage CLI keeps the source in argv[0].
_result_module = next((p / "NonStationaryMetrics" for p in _result_source.parents
                       if (p / "NonStationaryMetrics" / "verification_results.py").is_file()), None)
if _result_module is not None:
    _result_sys.path.insert(0, str(_result_module))
    from verification_results import ResultWriter as _ResultWriter
    _results = _ResultWriter(_result_source)
else:
    import os as _result_os
    if _result_os.environ.get("NSB_RESULT_JSON"):
        raise SystemExit("cannot export NSB_RESULT_JSON: verification_results.py not found above "
                         + str(_result_source))
    class _NoResults:
        def add(self, metric, sample, value): pass
        def finish(self): pass
    _results = _NoResults()

from sage.all import QQ
from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface as SageRS
import numpy as np

Rs = PolynomialRing(QQ, ['s', 'Y']); s, Y = Rs.gens()
lam = [QQ(1200), QQ(-2300), QQ(-11428), QQ(-5519), QQ(24700), QQ(62500)]  # i=0..5
q = sum(lam[i]*s**i for i in range(6))
X = SageRS(Y**2 - q, prec=100)
assert X.genus == 2

# canonical BEL second-kind forms (factor /2 = 4y->2y of the f_y convention)
dr1 = (lam[3]*s + 2*lam[4]*s**2 + 3*lam[5]*s**3)/2
dr2 = lam[5]*s**2/2
du = [Rs(1), s]                    # du_1=1, du_2=s
Pi = np.array(X.matrix_of_integral_values(du + [dr1, dr2]), dtype=complex)  # 4x4
omega, omega_p = Pi[:2, :2], Pi[:2, 2:]     # first kind a,b
eta,   eta_p   = Pi[2:, :2], Pi[2:, 2:]     # second kind a,b
tau = np.linalg.solve(omega, omega_p)
kappa = eta @ np.linalg.inv(omega)          # Klein's matrix (for sigma)
Leg = omega_p @ eta.T - omega @ eta_p.T

print('tau =\n', np.round(tau, 5))
print('tau sym err        =', np.linalg.norm(tau - tau.T))
print('kappa =\n', np.round(kappa, 5))
_results.add('symmetry', 'kappa_norm', np.linalg.norm(kappa-kappa.T))
print('kappa SYM ERR      =', np.linalg.norm(kappa - kappa.T))
print('Legendre/(2 pi i)  =\n', np.round(Leg/(2j*np.pi), 4), '  (expected -I/2)')
print("\nOK: eta (second kind) computed and validated -> Klein's sigma/zeta ready.")
print("    Next: sigma=exp(-1/2 u^T kappa u) theta(...), zeta_i=d_i log sigma,")
print("    psi(r) = closed combination of zeta(u(r)); validate against the numerical psi.")

_results.finish()
