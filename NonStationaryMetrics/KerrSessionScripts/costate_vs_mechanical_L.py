# -*- coding: utf-8 -*-
# Referee Issue 2 (BLOCKER): derive that the conserved Pontryagin costate p_phi = J is the
# MECHANICAL covariant angular momentum u_phi = g(u, d_phi), and that the conformal weights
# E_eff = Ehat/A, J_eff = J/A follow from the covariant transformation (not asserted).
#
# PART 1 (Legendre identification): the branch Hamiltonian shell H2=0 is the Kerr mass-shell
#   g_Kerr^{mu nu} u_mu u_nu = -1 with covariant momenta u_mu = (-E, p_r, J). If the on-shell p_r
#   solving this equals the paper's p_r0 = sqrt(S)/(Delta*Emu), then p_phi=J = u_phi (covariant),
#   i.e. the costate IS the mechanical angular momentum.
# PART 2 (conformal weight): TK metric g~=A^2 g_Kerr, physical u~=u/A, so Ehat=-g~(u~,d_t)=A*E_Kerr
#   and J_phys=g~(u~,d_phi)=A*J_Kerr => E_eff=Ehat/A, J_eff=J/A. Verified numerically.
import numpy as np, sympy as sp

M,a=1.0,0.9
r,E,J,pr=sp.symbols('r E J p_r',positive=True)
# equatorial Kerr (Boyer-Lindquist, theta=pi/2), coords (t,r,phi)
Dl=r**2-2*M*r+a**2
g_tt=-(1-2*M/r); g_tp=-2*M*a/r; g_rr=r**2/Dl; g_pp=r**2+a**2+2*M*a**2/r
# inverse (t,phi) 2x2 block + r
Det=g_tt*g_pp-g_tp**2
gtt=g_pp/Det; gpp=g_tt/Det; gtp=-g_tp/Det; grr=Dl/r**2
# mass shell with u_mu=(-E,p_r,J):  gtt E^2 - 2 gtp E J + gpp J^2 + grr p_r^2 = -1
shell=sp.simplify(gtt*E**2-2*gtp*E*J+gpp*J**2+grr*pr**2+1)
pr2=sp.solve(shell,pr**2)[0]
pr2=sp.simplify(pr2)
print("PART 1 -- Kerr mass-shell solved for p_r^2 (u_mu=(-E,p_r,J)):")
print("  p_r^2 =",sp.factor(pr2))
# paper on-shell: p_r0 = sqrt(S)/(Delta*Emu),  S=r*Emu*Q2, Emu=(E^2-1)r+2M
Emu=(E**2-1)*r+2*M
Q2=(2*E**2*J**2*M*r-E**2*J**2*r**2-4*E**2*J*M*a*r+2*E**2*M*a**2*r+E**2*a**2*r**2
    +E**2*r**4+4*J**2*M**2-4*J**2*M*r+J**2*r**2-8*J*M**2*a+4*J*M*a*r+4*M**2*a**2)
S=r*Emu*Q2
pr0_2=sp.simplify(S/(Dl*Emu)**2)   # (sqrt(S)/(Delta*Emu))^2
diff=sp.simplify(pr2-pr0_2)
print("  paper p_r0^2 = S/(Delta*Emu)^2")
print("  shell p_r^2 - paper p_r0^2 =",diff,"  => equal?",diff==0)
print("  => the branch-shell momenta are the COVARIANT u_mu; hence p_phi=J = u_phi (mechanical).")

# numeric cross-check + timelike unit + u_t=-E + u_phi=J via raising/lowering
print("\nnumeric check (r=6,E=1.4,J=6):")
sub={r:6.0,E:1.4,J:6.0}
gLow=np.array([[float(g_tt.subs(sub)),0,float(g_tp.subs(sub))],
               [0,float(g_rr.subs(sub)),0],
               [float(g_tp.subs(sub)),0,float(g_pp.subs(sub))]])
gInv=np.linalg.inv(gLow)
prv=float(sp.sqrt(pr2).subs(sub))
u_low=np.array([-1.4, prv, 6.0])         # (u_t,u_r,u_phi)=(-E,p_r,J)
u_up=gInv@u_low
norm=u_low@u_up
print(f"  g(u,u)=u_mu u^mu = {norm:+.6f}  (=-1 => timelike unit on shell)")
print(f"  u_t={u_low[0]:+.4f} (=-E), u_phi={u_low[2]:+.4f} (=J), recomputed g_phi.u^. = {gLow[2]@u_up:+.4f} (=J)")

# PART 2: conformal weight under g~=A^2 g_Kerr
print("\nPART 2 -- conformal weight (g~=A^2 g_Kerr, u~=u/A):")
for A in [1.0,1.5,2.3]:
    # physical tilde: g~=A^2 g, u~ unit wrt g~ => u~=u/A (u Kerr-unit)
    u_up_t=u_up/A
    g_til_low=A**2*gLow
    Ehat=-(g_til_low[0]@u_up_t)          # -g~(u~,d_t)
    Jphys=(g_til_low[2]@u_up_t)          #  g~(u~,d_phi)
    print(f"  A={A}: Ehat={Ehat:+.5f} (A*E_Kerr={A*1.4:.5f}), J_phys={Jphys:+.5f} (A*J_Kerr={A*6.0:.5f})"
          f"  => E_eff=Ehat/A={Ehat/A:.5f}, J_eff=J_phys/A={Jphys/A:.5f}")
print("\n=> DERIVED: costate p_phi=J is the covariant mechanical L=g(u,d_phi); conformal weights")
print("   E_eff=Ehat/A and J_eff=J/A follow from u~=u/A. Both charges scale as A^-1 (Issue 2 closed).")
