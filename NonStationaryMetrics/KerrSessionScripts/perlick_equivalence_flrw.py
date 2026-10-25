"""Frozen FLRW slice as a static Perlick metric: numerical check of the rail cost against Perlick's optical metric.

On g = a0^2 diag(-1, 1, 1) (constant a0) the rail cost F_rail, obtained by solving the rail and mass-shell
conditions for the time increment, is compared with Perlick's isotropic optical norm sqrt(a_ab v^a v^b),
a_ab = E^2 h_ab/(f(E^2 - f)), beta = 0, for a few (a0, v, E); the index is n = E/sqrt(E^2 - a0^2).
Exploratory check of the stationary limit used in Paper I.
"""
import numpy as np
# FLRW frozen slice: g = a0^2 diag(-1, 1, 1) in (eta, x, y). Static conformally-flat (constant a0).
# Perlick: f=-g_tt=a0^2, g_ta=0 => beta=0; h_ab=g_ab=a0^2 delta; a_ab=E^2 h_ab/(f(E^2-f))=E^2 delta/(E^2-a0^2).
# F_rail = sqrt(a_ab v^a v^b) = E|v|/sqrt(E^2-a0^2). Also FLRW optical index n=E/sqrt(E^2-a0^2), v_speed=sqrt(1-a0^2/E^2).
def metric(a0): return a0**2*np.array([[-1.,0,0],[0,1.,0],[0,0,1.]])  # (eta,x,y)
def F_rail(a0,vx,vy,E):
    g=metric(a0); gtt=g[0,0]; grr=g[1,1]; gpp=g[2,2]
    Bv=0.0; Cvv=grr*vx**2+gpp*vy**2; A=gtt
    aa=A**2+E**2*A; bb=2*(A+E**2)*Bv; cc=Bv**2+E**2*Cvv
    disc=bb**2-4*aa*cc
    ts=[(-bb+np.sqrt(disc))/(2*aa),(-bb-np.sqrt(disc))/(2*aa)]
    return max(t for t in ts if t>0)
def F_perlick(a0,vx,vy,E):
    f=a0**2; a_ab=E**2/(E**2-f)  # h_ab=a0^2 delta => a_ab=E^2 a0^2 delta/(a0^2(E^2-a0^2))=E^2/(E^2-a0^2)... wait
    # a_ab=E^2 h_ab/(f(E^2-f)), h=a0^2 delta, f=a0^2 => a_ab=E^2 a0^2/(a0^2(E^2-a0^2))=E^2/(E^2-a0^2)
    return np.sqrt(a_ab*(vx**2+vy**2))
print("FLRW frozen slice: F_rail (rail+shell) vs Perlick isotropic sqrt(a_ab v.v), beta=0:")
print(f"{'a0':>4} {'vx':>5} {'vy':>5} {'E':>4} {'F_rail':>11} {'F_perlick':>11} {'diff':>9} {'n=E/sqrt(E^2-a0^2)':>18}")
import itertools
for a0,vx,vy,E in itertools.product([0.8,1.2],[1.0,-0.5],[0.3],[1.4,2.0]):
    if E<=a0: continue
    Fr=F_rail(a0,vx,vy,E); Fp=F_perlick(a0,vx,vy,E); n=E/np.sqrt(E**2-a0**2)
    print(f"{a0:4.1f} {vx:5.1f} {vy:5.1f} {E:4.1f} {Fr:11.6f} {Fp:11.6f} {abs(Fr-Fp):9.1e} {n:18.6f}")
print("\n=> frozen FLRW = ISOTROPIC static Perlick metric (beta=0, purely Riemannian), index n=E/sqrt(E^2-a0^2).")
print("   full FLRW (a=a(eta)) = CONFORMAL generalisation (Fermat with a time-dependent index), level (b).")
