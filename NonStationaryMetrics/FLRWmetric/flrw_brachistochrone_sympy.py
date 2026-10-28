# -*- coding: utf-8 -*-
"""
t- and tau-brachistochrones in flat FLRW (k=0) -- SymPy derivation.

Formalism: covariantly constrained worldline (as in KerrMetric/doranTau.md),
extended to the non-stationary case through the conformal time eta.

Metric (c=1, signature -+++, motion in the plane z=0, comoving coordinates x,y):

    ds^2 = a(eta)^2 [ -deta^2 + dx^2 + dy^2 ]

d_eta is a CONFORMAL Killing field of the physical metric (a genuine Killing field of
ghat = eta_Mink). Rail constraint proposed in the note: Ehat = -u_eta held fixed.

Lagrangians (generic parameter lam, ' = d/dlam):

  tau branch:  L_tau = (1 + mu*Ehat) * Lam + mu * g_{eta nu} xdot^nu
                   = (1 + mu*Ehat) * Lam - mu * a^2 * eta'
  t branch:  L_t   = a * eta' + mu * ( Ehat*Lam - a^2 * eta' )
             (minimises T_t = int dt = int a deta, same constraint)

  with Lam = sqrt(-g_{ab} x'^a x'^b) = a * sqrt(eta'^2 - x'^2 - y'^2) = dtau/dlam.

Checks performed:
  V1  constraint  =>  gamma = Ehat/a  (E_phys = Ehat/a: redshift of the rail)
                      v(eta) = sqrt(1 - a^2/Ehat^2)  (local physical speed)
  V2  x,y cyclic => p_x, p_y conserved; p_y/p_x = y'/x' => comoving straight line
  V3  straight-line ansatz + mu(eta) solved from the x equation => the eta equation
      holds identically (gauge identity) -- consistent worldline
  V4  rail force: f_eta = a'/a = aH (conformal Hubble rate), f.u = 0,
      du_eta/dtau = 0 (constraint maintained)
  V5  reduction of the functionals on the solution:
        T_tau = (1/Ehat) * int a(eta)^2 deta ,   T_t = int a(eta) deta
      both monotone in the arrival eta => the tau branch and the t branch give the
      SAME curve (comoving straight line) in flat FLRW
  V6  t branch: same constraint, same straight line; mu_t(eta) = C_t/(Ehat^2 v)
  V7  limits: Minkowski a=1 (v = sqrt(1-1/Ehat^2), gamma=Ehat);
      de Sitter a = -1/(H eta): freezing at a=Ehat, maximum comoving distance
      in closed form
"""

import sympy as sp

lam = sp.Symbol('lam', real=True)
Eh = sp.Symbol('Ehat', positive=True)      # Ehat = -u_eta (conformal constraint)
H0 = sp.Symbol('H', positive=True)

eta = sp.Function('eta', real=True)(lam)
x = sp.Function('x', real=True)(lam)
y = sp.Function('y', real=True)(lam)
mu = sp.Function('mu', real=True)(lam)
a = sp.Function('a', positive=True)

etp, xp, yp = eta.diff(lam), x.diff(lam), y.diff(lam)

A = a(eta)
Lam = A * sp.sqrt(etp**2 - xp**2 - yp**2)   # dtau/dlam

def EL(L, q):
    return (sp.diff(L, q.diff(lam)).diff(lam) - sp.diff(L, q)).doit()

print("=" * 72)
print("TAU BRANCH:  L_tau = (1 + mu*Ehat)*Lam - mu*a^2*eta'")
print("=" * 72)

L_tau = (1 + mu * Eh) * Lam - mu * A**2 * etp

eq_mu = sp.diff(L_tau, mu)                 # constraint: Ehat*Lam - a^2*eta' = 0
eq_eta = EL(L_tau, eta)
eq_x = EL(L_tau, x)
eq_y = EL(L_tau, y)

print("\n[V1] constraint (dL/dmu = 0):")
print("     ", sp.simplify(eq_mu), " = 0")
print("      =>  a^2 * eta'/Lam = Ehat   that is  -u_eta = Ehat")

# gauge lam = eta: eta(lam) -> lam
gauge = [(etp, 1), (eta.diff(lam, 2), 0), (eta, lam)]
vx = x.diff(lam)   # now x' = dx/deta
vy = y.diff(lam)

con_g = sp.simplify(eq_mu.subs(gauge))
# constraint in the eta gauge: Ehat*a*sqrt(1 - v^2) = a^2  =>  v^2 = 1 - a^2/Ehat^2
v2_sol = sp.solve(con_g, vx**2 + vy**2)
print("\n[V1] gauge lam=eta:  v^2 = x'^2 + y'^2 =", sp.simplify(v2_sol[0]))
print("      => v(eta) = sqrt(1 - a^2/Ehat^2);  motion exists only if a < Ehat")
print("      gamma = dt/dtau = a*eta'/... :")
# u^t = dt/dtau = a * (deta/dtau); from the constraint a^2 deta/dtau = Ehat
print("      from the constraint: a^2 * deta/dtau = Ehat  =>  gamma = a*deta/dtau"
      " = Ehat/a")
print("      E_phys = gamma = Ehat/a  -> REDSHIFT of the rail (proportional to 1/a)")

print("\n[V2] x,y cyclic: p_x = dL/dx' conserved")
p_x = sp.diff(L_tau, xp)
p_y = sp.diff(L_tau, yp)
print("      p_x =", sp.simplify(p_x))
print("      eq_x is d(p_x)/dlam = 0:", sp.simplify(eq_x - p_x.diff(lam)) == 0)
print("      p_y/p_x = y'/x' = const  =>  comoving STRAIGHT LINE (without loss of generality: y=0)")

# ---- V3: straight-line ansatz y=0, x'(eta) = v(eta); solve for mu(eta) from eq_x ----
etav = sp.Symbol('etav', real=True)        # eta as the variable
av = a(etav)
v = sp.sqrt(1 - av**2 / Eh**2)             # speed from the constraint
muf = sp.Function('mu_sol', real=True)(etav)

# p_x in this gauge, on the ansatz:
p_x_g = p_x.subs(gauge).subs([(vy, 0), (y, 0)])
p_x_ans = p_x_g.subs([(vx, v), (x.diff(lam), v)]).subs(lam, etav) \
               .subs(mu.subs(lam, etav), muf)
p_x_ans = p_x_ans.subs(sp.Symbol('lam'), etav)
p_x_ans = sp.simplify(p_x_ans)
print("\n[V3] p_x on the straight-line ansatz:", p_x_ans)

C = sp.Symbol('C', real=True)              # C = -p_x, constant of motion
mu_sol = sp.solve(sp.Eq(p_x_ans, -C), muf)[0]
mu_sol = sp.simplify(mu_sol)
print("      mu(eta) solved from p_x = -C:")
print("      mu =", mu_sol)
print("      (1 + mu*Ehat) =", sp.simplify(1 + mu_sol * Eh), "  [= C/(Ehat*v)]")

# eta equation on the ansatz: it must vanish identically (gauge identity)
subs_ans = [(mu, muf.subs(etav, eta)), (yp, 0), (y, 0)]
eq_eta_ans = eq_eta.subs(subs_ans)
eq_eta_ans = eq_eta_ans.subs(gauge)
eq_eta_ans = eq_eta_ans.subs(vx, v.subs(etav, lam)).doit()
eq_eta_ans = eq_eta_ans.subs(x.diff(lam), v.subs(etav, lam)).doit()
eq_eta_ans = eq_eta_ans.subs(muf.subs(etav, lam), mu_sol.subs(etav, lam))
eq_eta_ans = eq_eta_ans.doit()
# the only remaining Subs are Subs(Derivative(a(w),w), w, lam) = a'(lam)
for s_ in eq_eta_ans.atoms(sp.Subs):
    eq_eta_ans = eq_eta_ans.xreplace({s_: a(lam).diff(lam)})
eq_eta_ans = sp.simplify(eq_eta_ans)
print("\n[V3] residual of the eta equation on the ansatz (expected 0):", eq_eta_ans)

# ---- V4: rail force from the geodesic residual ----
print("\n[V4] rail force  f_mu = du_mu/dtau - (1/2) d_mu g_ab u^a u^b")
tau_ = sp.Symbol('tau', real=True)
# on the solution: u_eta = -Ehat (const), u_x = Ehat*v(eta), dtau = a^2/Ehat deta
u_eta_l = -Eh
u_x_l = Eh * v
detadtau = Eh / av**2
du_eta = sp.diff(u_eta_l, etav) * detadtau            # = 0
du_x = sp.diff(u_x_l, etav) * detadtau
# u^a u_a e (1/2) d_eta g_ab u^a u^b = (a'/a) g_ab u^a u^b = -a'/a
u_up_eta = Eh / av**2
u_up_x = u_x_l / av**2
norm = sp.simplify(-av**2 * u_up_eta**2 + av**2 * u_up_x**2)
print("      u.u =", norm, " (expected -1)")
half_deta_g = (av.diff(etav) / av) * norm             # (1/2) d_eta g_ab u^a u^b
f_eta = sp.simplify(du_eta - half_deta_g)
f_x = sp.simplify(du_x - 0)
print("      f_eta =", f_eta, "   [expected a'/a = aH, conformal Hubble rate]")
print("      f_x   =", sp.simplify(f_x), "   [= -a'/(a v)]")
fdotu = sp.simplify((-1 / av**2) * f_eta * u_eta_l + (1 / av**2) * f_x * u_x_l)
print("      f.u   =", fdotu, " (expected 0: the rail does no proper work)")
print("      du_eta/dtau =", du_eta, " (constraint -u_eta = Ehat maintained)")

# ---- V5: reduction of the functionals ----
print("\n[V5] functionals on the solution (gauge lam = eta):")
Lam_g = (A * sp.sqrt(etp**2 - xp**2 - yp**2)).subs(gauge)
Lam_sol = Lam_g.subs([(vy, 0), (y, 0), (vx, v.subs(etav, lam))])
Lam_sol = sp.simplify(Lam_sol.subs(x.diff(lam), v.subs(etav, lam)))
print("      dtau/deta =", Lam_sol, "  =>  T_tau = (1/Ehat) int a^2 deta")
print("      dt/deta   = a           =>  T_t   =        int a   deta")
print("      both increasing in eta_arrival => same brachistochrone")
print("      (comoving straight line); they differ only in the VALUE of the time.")

print()
print("=" * 72)
print("T BRANCH:  L_t = a*eta' + mu*(Ehat*Lam - a^2*eta')")
print("=" * 72)

L_t = A * etp + mu * (Eh * Lam - A**2 * etp)

eq_mu_t = sp.diff(L_t, mu)
eq_x_t = EL(L_t, x)
p_x_t = sp.diff(L_t, xp)
p_y_t = sp.diff(L_t, yp)

print("\n[V6] identical constraint:", sp.simplify(eq_mu_t - eq_mu) == 0)
print("      p_x^t =", sp.simplify(p_x_t))
print("      p_y^t/p_x^t = y'/x' = const => same comoving STRAIGHT LINE")

p_x_t_g = p_x_t.subs(gauge).subs([(vy, 0), (y, 0)])
p_x_t_ans = p_x_t_g.subs([(vx, v), (x.diff(lam), v)]).subs(lam, etav) \
                   .subs(mu.subs(lam, etav), muf)
p_x_t_ans = sp.simplify(p_x_t_ans)
Ct = sp.Symbol('C_t', real=True)
mu_t_sol = sp.simplify(sp.solve(sp.Eq(p_x_t_ans, -Ct), muf)[0])
print("      mu_t(eta) =", mu_t_sol, "  [= C_t/(Ehat^2 v)]")
print("      => in FLRW k=0 the t and tau branches: same curve, different multipliers")

print()
print("=" * 72)
print("LIMITS AND CONCRETE CASES")
print("=" * 72)

print("\n[V7a] Minkowski a=1:")
v_mink = v.subs(av, 1)
print("      v =", sp.simplify(v_mink), "  gamma = Ehat  (standard, OK)")

print("\n[V7b] de Sitter  a(eta) = -1/(H eta),  eta in (-oo, 0):")
a_dS = -1 / (H0 * etav)
v_dS = sp.sqrt(1 - a_dS**2 / Eh**2)
print("      v(eta) =", v_dS)
eta_f = sp.solve(sp.Eq(a_dS, Eh), etav)[0]
print("      freezing (v=0) at  a=Ehat:  eta_f =", eta_f,
      "  t_f = ln(Ehat)/H")
# substitution u = -Ehat*H*eta (u>0):  eta0=-1/H -> u=Ehat,  eta_f -> u=1
eta0 = -1 / H0                                        # a(eta0) = 1
u = sp.Symbol('u', positive=True)
Dx = sp.integrate(sp.sqrt(1 - 1 / u**2), (u, 1, Eh)) / (Eh * H0)
Dx = sp.simplify(Dx)
if isinstance(Dx, sp.Piecewise):        # physical branch Ehat > 1
    Dx = Dx.args[-1].expr
print("      maximum comoving distance (from a=1):")
print("      Dx_max =", Dx)
print("      comoving Hubble horizon at eta0: r_H = 1/(a H) = 1/H")
ratio = sp.simplify(Dx * H0)
print("      Dx_max * H =", ratio, "  -> compare with r_H*H = 1")
print("      limit Ehat -> oo (null):",
      sp.limit(ratio, Eh, sp.oo))
print("      Ehat -> oo: Dx*H -> 1 : the null brachistochrone reaches")
print("      exactly the horizon (conformal consistency, Kovner-Perlick)")

print("\n[V7c] de Sitter travel times (from eta0=-1/H to eta1):")
eta1 = sp.Symbol('eta1', negative=True)
T_t_dS = sp.integrate(a_dS, (etav, eta0, eta1))
T_tau_dS = sp.integrate(a_dS**2 / Eh, (etav, eta0, eta1))
print("      T_t   =", sp.simplify(T_t_dS))
print("              (real: eta1<0, log(eta1)=log(-eta1)+i*pi;"
      " T_t = ln(a(eta1))/H, as t = ln(a)/H in de Sitter)")
print("      T_tau =", sp.simplify(T_tau_dS))
print("              T_tau(eta1 -> eta_f=-1/(Ehat*H)) =",
      sp.simplify(T_tau_dS.subs(eta1, -1 / (Eh * H0))),
      " FINITE: freezing happens at finite tau")

print("\n[V8] alternative PHYSICAL rail (constraint -u.n = gamma0, n comoving):")
print("      u^t = gamma0 const => v = sqrt(1-1/gamma0^2) constant,")
print("      dx/deta = v const: still a straight line, BUT no freezing;")
print("      E_conf = -u_eta = a*gamma0 GROWS with a: the rail must")
print("      pump energy against the redshift. The two rails (Ehat vs")
print("      gamma0) coincide only if a=const (Minkowski).")

print("\nDONE.")
