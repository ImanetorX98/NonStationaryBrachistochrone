# -*- coding: utf-8 -*-
"""Close pairs of fixed-endpoint tau-brachistochrones (Paper II, sec:bvp), at 50 digits.

Static Schwarzschild, M = 1, E = 6/5, r_A = r_B = 6, tau-branch optical index
N^2 = f/(E^2 - f).  With b = E^2 - 1 and J^2 = N^2 r_m^2 = r_m^2 (r_m - 2)/(b r_m + 2) at the
periapsis r_m, the swept angle and the optical cost are

    Phi(r_m)  = 2 int_{r_m}^{6} J dr / (r sqrt(f) sqrt(N^2 r^2 - J^2)),
    cost(r_m) = 2 int_{r_m}^{6} N^2 r dr / (sqrt(f) sqrt(N^2 r^2 - J^2)),

and the (r - r_m)^{-1/2} singularity is removed analytically (r = r_m + u^2, the factor
N^2 r^2 - J^2 divided exactly by r - r_m).  Reproduced here:

  (1) aperture 2.038982616263322: two roots J = 1.13366..., 1.13402..., 3.6e-4 apart,
      closer than the J-grid step 2.7e-3 of bvp_estremi_fissi.py, so a sign-change scan
      misses both;
  (2) aperture 2.038982626262098, between the three-sample vertex and the true maximum of
      Phi: roots J = 1.1338407..., 1.1338447..., on the same side of the sampled node;
  (3) aperture 2 x 1.01945: two solutions r_min = 2.585956..., 2.611228..., whose costs differ
      by 1.8e-6, the deeper one being the more expensive.
"""
import sys
import mpmath as m

m.mp.dps = 50
b, r0 = m.mpf('0.44'), m.mpf(6)
FAIL = []


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not ok:
        FAIL.append(name)


def J_of(rm):
    return m.sqrt(rm**2*(rm - 2)/(b*rm + 2))


def qdiv(r, rm):          # (N^2 r^2 - J^2)/(r - r_m) (b r + 2)(b r_m + 2) / (b r + 2)(b r_m + 2)
    return (b*rm*r*(rm + r - 2) + 2*(r*r + r*rm + rm*rm) - 4*(r + rm))/((b*r + 2)*(b*rm + 2))


def angle(rm):
    j = J_of(rm)
    return 2*m.quad(lambda u: 2*j/((rm + u*u)*m.sqrt((1 - 2/(rm + u*u))*qdiv(rm + u*u, rm))),
                    [0, m.sqrt(r0 - rm)])


def cost(rm):
    def g(u):
        r = rm + u*u; f = 1 - 2/r
        return 2*(f/(b + 1 - f))*r/m.sqrt(f*qdiv(r, rm))
    return 2*m.quad(g, [0, m.sqrt(r0 - rm)])


# the two forms of the angle agree: direct quadrature away from r_m as a cross-check
rm_t = m.mpf('2.6')
direct = 2*m.quad(lambda r: J_of(rm_t)/(r*m.sqrt(1 - 2/r)*m.sqrt((r - 2)*r*r/(b*r + 2) - J_of(rm_t)**2)),
                  [rm_t, rm_t + m.mpf('1e-6'), rm_t + m.mpf('0.01'), r0])
check("regularised angle equals the direct integral", abs(direct - angle(rm_t)) < m.mpf(10)**-8,
      f"diff {float(abs(direct - angle(rm_t))):.1e}")
peak = m.findroot(lambda r: m.diff(angle, r), m.mpf('2.598487'))
print(f"  maximum of Phi: r_m = {m.nstr(peak, 20)}, Phi = {m.nstr(angle(peak), 20)}")


def pair(target, seeds):
    out = []
    for s in seeds:
        rm = m.findroot(lambda r: angle(r) - target, tuple(map(m.mpf, s)), solver='secant')
        out.append(rm)
    return out


print("(1) aperture 2.038982616263322")
p1 = pair(m.mpf('2.038982616263322'), [('2.59834', '2.59835'), ('2.59862', '2.59863')])
J1 = [J_of(x) for x in p1]
print(f"      J = {m.nstr(J1[0], 12)}, {m.nstr(J1[1], 12)}; gap {m.nstr(J1[1] - J1[0], 3)}")
check("J = 1.13366, 1.13402, gap 3.6e-4 < grid step 2.7e-3",
      abs(J1[0] - m.mpf('1.13366')) < 1e-5 and abs(J1[1] - m.mpf('1.13402')) < 1e-5
      and J1[1] - J1[0] < m.mpf('2.7e-3'))
print("(2) aperture 2.038982626262098")
p2 = pair(m.mpf('2.038982626262098'), [('2.598484', '2.598486'), ('2.598488', '2.598490')])
J2 = [J_of(x) for x in p2]
print(f"      J = {m.nstr(J2[0], 12)}, {m.nstr(J2[1], 12)}")
check("J = 1.1338407, 1.1338447",
      abs(J2[0] - m.mpf('1.1338407')) < 1e-7 and abs(J2[1] - m.mpf('1.1338447')) < 1e-7)
print("(3) aperture 2 x 1.01945")
p3 = pair(2*m.mpf('1.01945'), [('2.585', '2.587'), ('2.610', '2.612')])
c3 = [cost(x) for x in p3]
print(f"      r_min = {m.nstr(p3[0], 12)}, {m.nstr(p3[1], 12)}; costs {m.nstr(c3[0], 15)}, {m.nstr(c3[1], 15)}")
check("r_min = 2.585956, 2.611228",
      abs(p3[0] - m.mpf('2.585956')) < 1e-6 and abs(p3[1] - m.mpf('2.611228')) < 1e-6)
check("costs differ by 1.8e-6, the deeper one more expensive",
      abs((c3[0] - c3[1]) - m.mpf('1.8e-6')) < m.mpf('0.05e-6'), f"diff {m.nstr(c3[0] - c3[1], 4)}")

print()
if FAIL:
    print(f"{len(FAIL)} FAILED:", *FAIL, sep="\n  ")
    sys.exit(1)
print("ALL CHECKS PASSED")
