# -*- coding: utf-8 -*-
# Perlick equivalence (referee Issue 1): the controlled-rail arrival-time Finsler F_rail equals a
# RANDERS metric = Perlick's fixed-energy stationary optical metric. Derivation:
#   rail u_t=-E, g(u,u)=-1 => quadratic for t'=dt/dsigma => F_rail = beta_a v^a + sqrt(a_ab v^a v^b),
#   beta_a = g_ta/f,  a_ab = E^2(g_ta g_tb + f g_ab)/(f^2 (E^2-f)),  f=-g_tt.
# Verify F_rail (solved from rail+shell) == Randers form, on equatorial Kerr (stationary axisymmetric).
import numpy as np
M,a=1.0,0.9
def metric(r):
    Dl=r**2-2*M*r+a**2
    g_tt=-(1-2*M/r); g_tp=-2*M*a/r; g_rr=r**2/Dl; g_pp=r**2+a**2+2*M*a**2/r
    return np.array([[g_tt,0,g_tp],[0,g_rr,0],[g_tp,0,g_pp]])   # (t,r,phi)
def F_rail(r,vr,vp,E):
    # solve rail quadratic for t' (future-directed), F=t' (see derivation)
    g=metric(r); gtt,gtp,gpp,grr=g[0,0],g[0,2],g[2,2],g[1,1]
    A=gtt; Bv=gtp*vp; Cvv=gpp*vp**2+grr*vr**2
    # (A^2+E^2 A)t'^2 + 2(A+E^2)Bv t' + (Bv^2+E^2 Cvv)=0
    aa=A**2+E**2*A; bb=2*(A+E**2)*Bv; cc=Bv**2+E**2*Cvv
    disc=bb**2-4*aa*cc
    if disc<0: return None
    ts=[(-bb+np.sqrt(disc))/(2*aa),(-bb-np.sqrt(disc))/(2*aa)]
    tp=[t for t in ts if t>0]
    return max(tp) if tp else None
def F_randers(r,vr,vp,E):
    g=metric(r); f=-g[0,0]; gtp=g[0,2]; grr=g[1,1]; gpp=g[2,2]
    beta_p=gtp/f
    a_rr=E**2*(0+f*grr)/(f**2*(E**2-f))
    a_pp=E**2*(gtp**2+f*gpp)/(f**2*(E**2-f))
    return beta_p*vp + np.sqrt(a_rr*vr**2+a_pp*vp**2)

print("F_rail (solved from rail+mass-shell) vs Randers form beta.v + sqrt(a.v.v):")
print(f"{'r':>5} {'vr':>5} {'vp':>5} {'E':>4} {'F_rail':>12} {'F_randers':>12} {'diff':>10}")
import itertools
for r,vr,vp,E in itertools.product([6.0,8.0],[-1.0,0.5],[0.2,-0.3],[1.4,2.0]):
    Fr=F_rail(r,vr,vp,E); Fra=F_randers(r,vr,vp,E)
    if Fr is None: continue
    print(f"{r:5.1f} {vr:5.1f} {vp:5.1f} {E:4.1f} {Fr:12.6f} {Fra:12.6f} {abs(Fr-Fra):10.1e}")
print("\n=> F_rail == Randers form: the controlled-rail arrival-time Finsler IS Perlick's fixed-energy")
print("   stationary optical metric (Randers). beta_a=g_ta/f (frame-drag 1-form), a_ab=E^2(g_ta g_tb+f g_ab)/(f^2(E^2-f)).")
print("   Null/high-energy limit E->inf: a_ab -> g_ab/f (+beta) = standard Fermat optical metric.")
