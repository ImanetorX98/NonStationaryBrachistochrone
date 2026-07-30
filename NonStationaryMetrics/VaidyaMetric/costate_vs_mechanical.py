# -*- coding: utf-8 -*-
# main14 memo 4.3: the axial control costate J = p_phi is NOT the mechanical angular
# momentum L_mech = g(u, d_phi). Along the forced rail worldline L_mech is not conserved.
#
# From the Vaidya v-branch indicatrix (eq:vaidya-ind), reconstructing the four-velocity
# and the (unit-cost) costate by the Legendre map, we verify the exact relation
#
#     L_mech = r^2 u^phi = ((Ehat^2 - f(r))/Ehat) * J
#
# which varies with r while J is constant -- so the two objects are distinct (they agree
# only at the isolated radius f = Ehat^2 - Ehat). This corrects the earlier (wrong)
# statement that the costate and the mechanical momentum "coincide here".
import sympy as sp

r, f, Eh, th = sp.symbols('r f Ehat theta', positive=True)
w = Eh**2 - f                                   # w = Ehat^2 - f

# indicatrix velocities (eq:vaidya-ind), parameter = advanced clock v
rp = (f - Eh**2) + Eh*sp.sqrt(w)*sp.cos(th)     # dr/dv
php = sp.sqrt(w)/r*sp.sin(th)                    # dphi/dv

# four-velocity from the rail -u_v = Ehat, with u^r = rp u^v (metric g_vv=-f, g_vr=1)
uv = Eh/(f - rp)
uphi = php*uv                                    # u^phi = dphi/dtau
Lmech = r**2*uphi                                # g_phiphi u^phi (equatorial)

# unit-cost costate p_phi = J: outward normal to the indicatrix, normalized so p.xdot=1
kappa = r/(w*(Eh - sp.sqrt(w)*sp.cos(th)))
pr = kappa*(sp.sqrt(w)/r*sp.cos(th))
J = kappa*Eh*sp.sqrt(w)*sp.sin(th)

print("p.xdot (unit cost, expect 1) =", sp.simplify(pr*rp + J*php))
print("L_mech / J                   =", sp.simplify(Lmech/J))
print("expected (Ehat^2-f)/Ehat     =", sp.simplify(w/Eh))
print("difference (expect 0)        =", sp.simplify(Lmech/J - w/Eh))
