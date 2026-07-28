# -*- coding: utf-8 -*-
# Generic special-function closed form of the Vaidya v-branch off-shell correction.
# Source is Theta = m d_m (NO conformal dilation letter), so the inner letter is
#   Sigma(r) = int (Theta H / H_pr) dr = int (m H_m / H_pr) dr   (NOT the radial action).
# On the frozen shell H_v=0:  p_r0 = sqrt(S_V)/den,  S_V = r*Emu*Q2(a=0) (genus-2),
# Emu=(E^2-1)r+2M.  Off-shell wrap  Phi = -int (G_pr/H_pr) Sigma dr,  (G_pr/H_pr)=A_V/sqrt(S_V).
# STEP 1 here: build S_V, p_r0, verify kernel=A_V/sqrt(S_V), reduce inner letter
#   m H_m/H_pr = (rational)/sqrt(S_V) + third-kind at the Schwarzschild horizon r=2m.
import sympy as sp

r,pr,m,E,J=sp.symbols('r p_r m E J',positive=True)
f=1-2*m/r; w=E**2-f
Rad=sp.sqrt(E**2*pr**2+J**2/r**2)
Hv=pr*(f-E**2)-1+sp.sqrt(w)*Rad
# on-shell p_r: sqrt(w)*Rad = 1-pr(f-E^2) => square
lhs2=w*(E**2*pr**2+J**2/r**2); rhs=1-pr*(f-E**2)
pr2sol=sp.solve(sp.expand(lhs2-rhs**2),pr**2)[0]
pr2sol=sp.simplify(pr2sol)
print("Vaidya on-shell p_r^2 =",sp.factor(pr2sol))
# curve: S_V = r*Emu*Q2(a=0); verify p_r0^2 = S_V/den^2 for some polynomial den
Emu=(E**2-1)*r+2*m
Q2v=(2*E**2*J**2*m*r-E**2*J**2*r**2+E**2*r**4+4*J**2*m**2-4*J**2*m*r+J**2*r**2)  # Q2 at a=0
S_V=sp.expand(r*Emu*Q2v)
den2=sp.simplify(S_V/pr2sol)
print("\nS_V / p_r0^2 = den^2 =",sp.factor(den2))
print("  => den =",sp.sqrt(sp.factor(den2)) if sp.sqrt(den2).is_polynomial else sp.factor(den2))

# H_m and H_pr on shell (use p_r0)
Hm=sp.diff(Hv,m); Hpr=sp.diff(Hv,pr)
# on shell substitute sqrt(w)*Rad = rho := 1 - pr(f-E^2)
rho=1-pr*(f-E**2)
Hm_on=sp.simplify(Hm.subs(sp.sqrt(w)*Rad, rho).rewrite(sp.Pow))
# safer: express Rad = rho/sqrt(w) on shell
Rad_on=rho/sp.sqrt(w)
Hm_on=sp.simplify(pr*sp.diff(f,m)+ (sp.diff(sp.sqrt(w),m))*Rad_on + sp.sqrt(w)*0)  # Rad has no m
Hpr_on=sp.simplify((f-E**2)+sp.sqrt(w)*sp.diff(Rad,pr).subs(Rad,Rad_on))
ThetaH_on=sp.simplify(m*Hm_on)
inner_integrand=sp.simplify(ThetaH_on/Hpr_on)   # = Theta H / H_pr on shell (fn of r, pr)
print("\ninner letter integrand Theta H/H_pr (on shell, before p_r0 sub):")
print("  ",inner_integrand)
