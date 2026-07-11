# -*- coding: utf-8 -*-
"""
Brachistocrone t e tau in FLRW piatta (k=0) — derivazione sympy.

Formalismo: worldline vincolato covariante (come KerrMetric/doranTau.md),
esteso al caso non stazionario via tempo conforme eta.

Metrica (c=1, segnatura -+++, moto nel piano z=0, coordinate comoventi x,y):

    ds^2 = a(eta)^2 [ -deta^2 + dx^2 + dy^2 ]

d_eta e' Killing CONFORME della metrica fisica (Killing genuino di
ghat = eta_Mink). Vincolo rotaia proposto dal doc: Ehat = -u_eta fissata.

Lagrangiane (parametro generico lam, ' = d/dlam):

  ramo tau:  L_tau = (1 + mu*Ehat) * Lam + mu * g_{eta nu} xdot^nu
                   = (1 + mu*Ehat) * Lam - mu * a^2 * eta'
  ramo t:    L_t   = a * eta' + mu * ( Ehat*Lam - a^2 * eta' )
             (minimizza T_t = int dt = int a deta, stesso vincolo)

  con Lam = sqrt(-g_{ab} x'^a x'^b) = a * sqrt(eta'^2 - x'^2 - y'^2) = dtau/dlam.

Verifiche eseguite:
  V1  vincolo  =>  gamma = Ehat/a  (E_phys = Ehat/a: redshift della rotaia)
                   v(eta) = sqrt(1 - a^2/Ehat^2)  (velocita' locale fisica)
  V2  x,y cicliche => p_x, p_y conservati; p_y/p_x = y'/x' => retta comovente
  V3  ansatz retta + mu(eta) risolto dall'eq. di x => eq. di eta soddisfatta
      identicamente (identita' di gauge) — worldline consistente
  V4  forza-rotaia: f_eta = a'/a = aH (conformal Hubble), f.u = 0,
      du_eta/dtau = 0 (vincolo mantenuto)
  V5  riduzione funzionali sulla soluzione:
        T_tau = (1/Ehat) * int a(eta)^2 deta ,   T_t = int a(eta) deta
      entrambi monotoni nell'eta d'arrivo => ramo tau e ramo t danno la
      STESSA curva (retta comovente) in FLRW piatta
  V6  ramo t: stesso vincolo, stessa retta; mu_t(eta) = C_t/(Ehat^2 v)
  V7  limiti: Minkowski a=1 (v = sqrt(1-1/Ehat^2), gamma=Ehat);
      de Sitter a = -1/(H eta): congelamento a a=Ehat, distanza comovente
      massima in forma chiusa
"""

import sympy as sp

lam = sp.Symbol('lam', real=True)
Eh = sp.Symbol('Ehat', positive=True)      # Ehat = -u_eta (vincolo conforme)
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
print("RAMO TAU:  L_tau = (1 + mu*Ehat)*Lam - mu*a^2*eta'")
print("=" * 72)

L_tau = (1 + mu * Eh) * Lam - mu * A**2 * etp

eq_mu = sp.diff(L_tau, mu)                 # vincolo: Ehat*Lam - a^2*eta' = 0
eq_eta = EL(L_tau, eta)
eq_x = EL(L_tau, x)
eq_y = EL(L_tau, y)

print("\n[V1] vincolo (dL/dmu = 0):")
print("     ", sp.simplify(eq_mu), " = 0")
print("      =>  a^2 * eta'/Lam = Ehat   cioe'  -u_eta = Ehat")

# gauge lam = eta: eta(lam) -> lam
gauge = [(etp, 1), (eta.diff(lam, 2), 0), (eta, lam)]
vx = x.diff(lam)   # ora x' = dx/deta
vy = y.diff(lam)

con_g = sp.simplify(eq_mu.subs(gauge))
# vincolo in gauge eta: Ehat*a*sqrt(1 - v^2) = a^2  =>  v^2 = 1 - a^2/Ehat^2
v2_sol = sp.solve(con_g, vx**2 + vy**2)
print("\n[V1] gauge lam=eta:  v^2 = x'^2 + y'^2 =", sp.simplify(v2_sol[0]))
print("      => v(eta) = sqrt(1 - a^2/Ehat^2);  moto esiste solo se a < Ehat")
print("      gamma = dt/dtau = a*eta'/... :")
# u^t = dt/dtau = a * (deta/dtau); dal vincolo a^2 deta/dtau = Ehat
print("      dal vincolo: a^2 * deta/dtau = Ehat  =>  gamma = a*deta/dtau"
      " = Ehat/a")
print("      E_phys = gamma = Ehat/a  -> REDSHIFT della rotaia (prop 1/a)")

print("\n[V2] x,y cicliche: p_x = dL/dx' conservato")
p_x = sp.diff(L_tau, xp)
p_y = sp.diff(L_tau, yp)
print("      p_x =", sp.simplify(p_x))
print("      eq_x e' d(p_x)/dlam = 0:", sp.simplify(eq_x - p_x.diff(lam)) == 0)
print("      p_y/p_x = y'/x' = cost  =>  RETTA comovente (senza perdita: y=0)")

# ---- V3: ansatz retta y=0, x'(eta) = v(eta), risolvo mu(eta) da eq_x ----
etav = sp.Symbol('etav', real=True)        # eta come variabile
av = a(etav)
v = sp.sqrt(1 - av**2 / Eh**2)             # velocita' dal vincolo
muf = sp.Function('mu_sol', real=True)(etav)

# p_x in gauge, sull'ansatz:
p_x_g = p_x.subs(gauge).subs([(vy, 0), (y, 0)])
p_x_ans = p_x_g.subs([(vx, v), (x.diff(lam), v)]).subs(lam, etav) \
               .subs(mu.subs(lam, etav), muf)
p_x_ans = p_x_ans.subs(sp.Symbol('lam'), etav)
p_x_ans = sp.simplify(p_x_ans)
print("\n[V3] p_x sull'ansatz retta:", p_x_ans)

C = sp.Symbol('C', real=True)              # C = -p_x costante del moto
mu_sol = sp.solve(sp.Eq(p_x_ans, -C), muf)[0]
mu_sol = sp.simplify(mu_sol)
print("      mu(eta) risolto da p_x = -C:")
print("      mu =", mu_sol)
print("      (1 + mu*Ehat) =", sp.simplify(1 + mu_sol * Eh), "  [= C/(Ehat*v)]")

# eq. di eta sull'ansatz: deve annullarsi identicamente (identita' di gauge)
subs_ans = [(mu, muf.subs(etav, eta)), (yp, 0), (y, 0)]
eq_eta_ans = eq_eta.subs(subs_ans)
eq_eta_ans = eq_eta_ans.subs(gauge)
eq_eta_ans = eq_eta_ans.subs(vx, v.subs(etav, lam)).doit()
eq_eta_ans = eq_eta_ans.subs(x.diff(lam), v.subs(etav, lam)).doit()
eq_eta_ans = eq_eta_ans.subs(muf.subs(etav, lam), mu_sol.subs(etav, lam))
eq_eta_ans = eq_eta_ans.doit()
# gli unici Subs residui sono Subs(Derivative(a(w),w), w, lam) = a'(lam)
for s_ in eq_eta_ans.atoms(sp.Subs):
    eq_eta_ans = eq_eta_ans.xreplace({s_: a(lam).diff(lam)})
eq_eta_ans = sp.simplify(eq_eta_ans)
print("\n[V3] residuo eq. di eta sull'ansatz (atteso 0):", eq_eta_ans)

# ---- V4: forza-rotaia dal residuo geodetico ----
print("\n[V4] forza-rotaia  f_mu = du_mu/dtau - (1/2) d_mu g_ab u^a u^b")
tau_ = sp.Symbol('tau', real=True)
# sulla soluzione: u_eta = -Ehat (cost), u_x = Ehat*v(eta), dtau = a^2/Ehat deta
u_eta_l = -Eh
u_x_l = Eh * v
detadtau = Eh / av**2
du_eta = sp.diff(u_eta_l, etav) * detadtau            # = 0
du_x = sp.diff(u_x_l, etav) * detadtau
# u^a u_a e (1/2) d_eta g_ab u^a u^b = (a'/a) g_ab u^a u^b = -a'/a
u_up_eta = Eh / av**2
u_up_x = u_x_l / av**2
norm = sp.simplify(-av**2 * u_up_eta**2 + av**2 * u_up_x**2)
print("      u.u =", norm, " (atteso -1)")
half_deta_g = (av.diff(etav) / av) * norm             # (1/2) d_eta g_ab u^a u^b
f_eta = sp.simplify(du_eta - half_deta_g)
f_x = sp.simplify(du_x - 0)
print("      f_eta =", f_eta, "   [atteso a'/a = aH, conformal Hubble]")
print("      f_x   =", sp.simplify(f_x), "   [= -a'/(a v)]")
fdotu = sp.simplify((-1 / av**2) * f_eta * u_eta_l + (1 / av**2) * f_x * u_x_l)
print("      f.u   =", fdotu, " (atteso 0: rotaia non compie lavoro proprio)")
print("      du_eta/dtau =", du_eta, " (vincolo -u_eta = Ehat mantenuto)")

# ---- V5: riduzione dei funzionali ----
print("\n[V5] funzionali sulla soluzione (gauge lam = eta):")
Lam_g = (A * sp.sqrt(etp**2 - xp**2 - yp**2)).subs(gauge)
Lam_sol = Lam_g.subs([(vy, 0), (y, 0), (vx, v.subs(etav, lam))])
Lam_sol = sp.simplify(Lam_sol.subs(x.diff(lam), v.subs(etav, lam)))
print("      dtau/deta =", Lam_sol, "  =>  T_tau = (1/Ehat) int a^2 deta")
print("      dt/deta   = a           =>  T_t   =        int a   deta")
print("      entrambi crescenti in eta_arrivo => stessa brachistocrona")
print("      (retta comovente); differiscono solo nel VALORE del tempo.")

print()
print("=" * 72)
print("RAMO T:  L_t = a*eta' + mu*(Ehat*Lam - a^2*eta')")
print("=" * 72)

L_t = A * etp + mu * (Eh * Lam - A**2 * etp)

eq_mu_t = sp.diff(L_t, mu)
eq_x_t = EL(L_t, x)
p_x_t = sp.diff(L_t, xp)
p_y_t = sp.diff(L_t, yp)

print("\n[V6] vincolo identico:", sp.simplify(eq_mu_t - eq_mu) == 0)
print("      p_x^t =", sp.simplify(p_x_t))
print("      p_y^t/p_x^t = y'/x' = cost => stessa RETTA comovente")

p_x_t_g = p_x_t.subs(gauge).subs([(vy, 0), (y, 0)])
p_x_t_ans = p_x_t_g.subs([(vx, v), (x.diff(lam), v)]).subs(lam, etav) \
                   .subs(mu.subs(lam, etav), muf)
p_x_t_ans = sp.simplify(p_x_t_ans)
Ct = sp.Symbol('C_t', real=True)
mu_t_sol = sp.simplify(sp.solve(sp.Eq(p_x_t_ans, -Ct), muf)[0])
print("      mu_t(eta) =", mu_t_sol, "  [= C_t/(Ehat^2 v)]")
print("      => in FLRW k=0 rami t e tau: stessa curva, moltiplicatori diversi")

print()
print("=" * 72)
print("LIMITI E CASI CONCRETI")
print("=" * 72)

print("\n[V7a] Minkowski a=1:")
v_mink = v.subs(av, 1)
print("      v =", sp.simplify(v_mink), "  gamma = Ehat  (standard, OK)")

print("\n[V7b] de Sitter  a(eta) = -1/(H eta),  eta in (-oo, 0):")
a_dS = -1 / (H0 * etav)
v_dS = sp.sqrt(1 - a_dS**2 / Eh**2)
print("      v(eta) =", v_dS)
eta_f = sp.solve(sp.Eq(a_dS, Eh), etav)[0]
print("      congelamento (v=0) a  a=Ehat:  eta_f =", eta_f,
      "  t_f = ln(Ehat)/H")
# sostituzione u = -Ehat*H*eta (u>0):  eta0=-1/H -> u=Ehat,  eta_f -> u=1
eta0 = -1 / H0                                        # a(eta0) = 1
u = sp.Symbol('u', positive=True)
Dx = sp.integrate(sp.sqrt(1 - 1 / u**2), (u, 1, Eh)) / (Eh * H0)
Dx = sp.simplify(Dx)
if isinstance(Dx, sp.Piecewise):        # ramo fisico Ehat > 1
    Dx = Dx.args[-1].expr
print("      distanza comovente massima (da a=1):")
print("      Dx_max =", Dx)
print("      orizzonte comovente di Hubble a eta0: r_H = 1/(a H) = 1/H")
ratio = sp.simplify(Dx * H0)
print("      Dx_max * H =", ratio, "  -> confronto con r_H*H = 1")
print("      limite Ehat -> oo (nullo):",
      sp.limit(ratio, Eh, sp.oo))
print("      Ehat -> oo: Dx*H -> 1 : la brachistocrona nulla raggiunge")
print("      esattamente l'orizzonte (consistenza conforme, Kovner-Perlick)")

print("\n[V7c] tempi di percorrenza de Sitter (da eta0=-1/H a eta1):")
eta1 = sp.Symbol('eta1', negative=True)
T_t_dS = sp.integrate(a_dS, (etav, eta0, eta1))
T_tau_dS = sp.integrate(a_dS**2 / Eh, (etav, eta0, eta1))
print("      T_t   =", sp.simplify(T_t_dS))
print("              (reale: eta1<0, log(eta1)=log(-eta1)+i*pi;"
      " T_t = ln(a(eta1))/H, come t = ln(a)/H in de Sitter)")
print("      T_tau =", sp.simplify(T_tau_dS))
print("              T_tau(eta1 -> eta_f=-1/(Ehat*H)) =",
      sp.simplify(T_tau_dS.subs(eta1, -1 / (Eh * H0))),
      " FINITO: il congelamento avviene a tau finito")

print("\n[V8] rotaia FISICA alternativa (vincolo -u.n = gamma0, n comovente):")
print("      u^t = gamma0 cost => v = sqrt(1-1/gamma0^2) costante,")
print("      dx/deta = v cost: ancora retta, MA nessun congelamento;")
print("      E_conf = -u_eta = a*gamma0 CRESCE con a: la rotaia deve")
print("      pompare energia contro il redshift. Le due rotaie (Ehat vs")
print("      gamma0) coincidono solo se a=cost (Minkowski).")

print("\nFATTO.")
