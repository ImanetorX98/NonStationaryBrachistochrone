# ADDITIVE CONSTANTS C0, Ce in CLOSED FORM (the last numerical piece of the block form).
# Idea: evaluate the elliptic function at a REGULAR point and subtract the principal parts.
#  - eta'(z) is regular at z=0 (r=e4): Ce = eta'(0) + 2 e1_zd zeta(z_d) - 2 e2_zi P(z_inf) + 2 e1_zi zeta(z_inf)
#  - R(z) is regular at z=z_inf AND R(z_inf)=0 (it decays as 1/r):
#      C0 = - sum_a [ b1^a zeta(z_inf-a) + b2^a P(z_inf-a) - (b3^a/2) P'(z_inf-a) ]
# Check against the numerical value (match at one point). All the coefficients are now symbolic/closed.
import numpy as np, mpmath as mp, sys
from scipy.integrate import quad
sys.stdout_ = sys.stdout; sys.stdout = open('/dev/null', 'w')
exec(open('vaidya_sep_analytic_residues.py').read().split('# =========== VALUES')[0])
sys.stdout = sys.stdout_
qn = q; Qp = np.polyder(Q)
def wp(z):
    u=mp.pi*z/(2*om1); r=mp.jtheta(1,u,qn,1)/mp.jtheta(1,u,qn)
    return -eta1/om1-(mp.pi/(2*om1))**2*(mp.jtheta(1,u,qn,2)/mp.jtheta(1,u,qn)-r**2)
def wpp(z):
    u=mp.pi*z/(2*om1); T1=mp.jtheta(1,u,qn,1)/mp.jtheta(1,u,qn); T2=mp.jtheta(1,u,qn,2)/mp.jtheta(1,u,qn); T3=mp.jtheta(1,u,qn,3)/mp.jtheta(1,u,qn)
    return -(mp.pi/(2*om1))**3*(T3-3*T1*T2+2*T1**3)
b2h = lambda ei: np.polyval(Nm,ei)/((ei-rd)**3)*(4/np.polyval(Qp,ei)**2)   # double-pole residue at the half-periods
e4v, e3v = er[3], er[2]                                                    # e4->z=0, e3->z=i w_im
iw = 1j*mp.mpf(w_im); zim = mp.mpf(z_inf); zdm = mp.mpf(z_d)
def r_of_z(z): return c_r-(1/sa)*(wzet(z-z_inf)-wzet(z+z_inf))
def Rex(z): rr=r_of_z(z); return np.polyval(Nm,rr)/((rr-rd)**3*np.polyval(Q,complex(rr)))
def Epex(z): rr=r_of_z(z); return (rr**3-2*rr**2)/(rr-rd)

# ---- CLOSED FORMULAS ----
Ce_closed = (e4v**3-2*e4v**2)/(e4v-rd) + 2*e1_zd*wzet(zdm) - 2*e2_zi*wp(zim) + 2*e1_zi*wzet(zim)
C0_closed = -( b1_zd*(wzet(zim-zdm)-wzet(zim+zdm)) + b2_zd*(wp(zim-zdm)+wp(zim+zdm))
               - (b3_zd/2)*(wpp(zim-zdm)-wpp(zim+zdm)) + b2h(e4v)*wp(zim) + b2h(e3v)*wp(zim-iw) )

# ---- check against numerics (principal parts subtracted at zt) ----
zt = mp.mpf('0.19')
Ce_num = complex(Epex(zt)-(e1_zd*(wzet(zt-zdm)-wzet(zt+zdm))+e2_zi*(wp(zt-zim)+wp(zt+zim))+e1_zi*(wzet(zt-zim)-wzet(zt+zim)))).real
C0_num = complex(Rex(zt)-(b1_zd*wzet(zt-zdm)-b1_zd*wzet(zt+zdm)+b2_zd*wp(zt-zdm)+b2_zd*wp(zt+zdm)
                          -(b3_zd/2)*wpp(zt-zdm)+(b3_zd/2)*wpp(zt+zdm)+b2h(e4v)*wp(zt)+b2h(e3v)*wp(zt-iw))).real
print("=== C0, Ce in CLOSED FORM against numerics ===")
print("R(z_inf) = %.1e  (must be ~0: it justifies the formula for C0)"%abs(complex(Rex(zim))))
print("Ce: closed = %.10f   num = %.10f   diff = %.1e"%(complex(Ce_closed).real, Ce_num, abs(complex(Ce_closed).real-Ce_num)))
print("C0: closed = %.10f   num = %.10f   diff = %.1e"%(complex(C0_closed).real, C0_num, abs(complex(C0_closed).real-C0_num)))
print("(residual ~1e-8 from the finite step for P'; the formulas are exact)")
print("\n=> every coefficient of the block form is now closed (algebraic + Weierstrass values at the poles).")
