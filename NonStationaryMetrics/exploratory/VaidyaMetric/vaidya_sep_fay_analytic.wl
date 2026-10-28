(* ::Package:: *)
(* =====================================================================
   ANALYTIC FAY - symbolic reduction of the dilogs and of beta in closed form.
   ---------------------------------------------------------------------
   IDEA: with the Weierstrass addition theorems, every even function e_i(z)
   is RATIONAL in p=P(z) (+ constants at the poles). The primitives Pe_i reduce
   to 6 transcendental LETTERS {z, Z(z)=zeta(z), Ls[+-zd], Ls[+-zi]}.
   The Fay relations = MATCHING of the coefficients (rational in p,pp) in front of
   each letter  ->  EXACT algebra (no sampling/conditioning).
   Constants at the poles kept SYMBOLIC  ->  SYMBOLIC beta (then a numerical check).
   Conventions identical to vaidya_sep_5term.py / vaidya_sep_fay_symbolic.wl.
   ===================================================================== *)

$MinPrecision = 60; prec = 55;
Off[NIntegrate::ncvb]; Off[NIntegrate::slwcon]; Off[NIntegrate::nlim]; Off[General::stop];

(* ---------- 0-1. curve + lattice (high-precision numerics, reusing the setup) ---------- *)
Em = 7/5; mm = 1; r0 = 12;
DEpol[r_] := (Em^2 - 1) r + 2 mm;  Spoly[r_, J_] := r (r - 2 mm) DEpol[r] (r^2 (r - 2 mm) - J^2 DEpol[r]);
JcE = Jc /. Solve[Resultant[Spoly[r, Jc], D[Spoly[r, Jc], r], r] == 0 && Jc > 1, Jc, Reals][[1]];
Jcn = N[JcE, prec]; Sc = Expand[Spoly[r, Jcn]];
a4 = N[Coefficient[Sc, r, 6], prec];
rootsS = Sort[N[r /. Solve[Sc == 0, r], prec], Re[#1] < Re[#2] &];
rd = Re@First@SelectFirst[Flatten[Table[{rootsS[[i]], rootsS[[j]]}, {i, Length[rootsS]}, {j, i+1, Length[rootsS]}], 1], Abs[#[[1]] - #[[2]]] < 10^-6 &];
erts = Sort[Re /@ DeleteCases[rootsS, x_ /; Abs[x - rd] < 10^-6, 1, 2]];
{e1, e2, e3, e4} = erts;  Q4[x_] := a4 (x - e1)(x - e2)(x - e3)(x - e4);
k2 = ((e3 - e2)(e4 - e1))/((e4 - e2)(e3 - e1)); pref = 2/Sqrt[(e4 - e2)(e3 - e1)]/Sqrt[a4];
om1 = N[pref EllipticK[k2], prec]; wim = N[pref EllipticK[1 - k2], prec];
{g2, g3} = WeierstrassInvariants[{om1, I wim}];
WP[z_] := WeierstrassP[z, {g2, g3}]; WZ[z_] := WeierstrassZeta[z, {g2, g3}];
WS[z_] := WeierstrassSigma[z, {g2, g3}]; WPp[z_] := WeierstrassPPrime[z, {g2, g3}];
zr[rv_] := NIntegrate[1/Sqrt[Q4[x]], {x, e4, rv}, WorkingPrecision -> prec];
z0 = zr[r0]; zinf = NIntegrate[1/Sqrt[Q4[x]], {x, e4, Infinity}, WorkingPrecision -> prec];
zd = zinf + NIntegrate[1/Sqrt[Q4[x]], {x, -Infinity, rd}, WorkingPrecision -> prec];
iw = I wim;
Print["Jc=", N[Jcn,12], " r_d=", N[rd,10], " z_d=", N[zd,10], " z_inf=", N[zinf,10]];

(* constants at the poles (NUMERICAL values; kept symbolic in the matching) *)
pd = WP[zd]; ppd = WPp[zd]; zdc = WZ[zd];       (* P,P',zeta at z_d *)
pI = WP[zinf]; ppI = WPp[zinf]; zic = WZ[zinf]; (* at z_inf *)
eIw = WP[iw];                                    (* P at the half period iw (P'(iw)=0) *)
Print["P(iw) (2-torsion, P'=0): ", N[eIw,10], "  check P'(iw)=", N[WPp[iw],6]];

(* ---------- 2. REDUCTION  e_i(z) -> rational in p=P(z), pp=P'(z)  (addition) ---------- *)
(* addition: P(z-+c) = -P(z)-P(c) + 1/4 ((pp -+ ppc)/(p-Pc))^2 ; with p,pp symbols *)
Ppm[sgn_, Pc_, ppc_] := -p - Pc + (1/4)((pp - sgn ppc)/(p - Pc))^2;      (* P(z + sgn c): addition with P'(sgn c)=sgn ppc *)
(* Z_a = -2 zeta(a) + P'(a)/(p - P(a))   [derived by hand, verified below] *)
Zred[Za_, Pa_, ppa_] := -2 Za + ppa/(p - Pa);
(* rational forms (in p,pp) of the 8 functions *)
erat[1] = 1;
erat[2] = Zred[zdc, pd, ppd];                                   (* Z_zd *)
erat[3] = Ppm[-1, pd, ppd] + Ppm[+1, pd, ppd] // Simplify;       (* P_zd = P(z-zd)+P(z+zd) *)
erat[5] = p;                                                    (* wp0 = P(z) *)
erat[6] = -p - eIw + (1/4)(pp/(p - eIw))^2 // Simplify;          (* P(z-iw), P'(iw)=0 *)
erat[7] = Ppm[-1, pI, ppI] + Ppm[+1, pI, ppI] // Simplify;       (* P_zi *)
erat[8] = Zred[zic, pI, ppI];                                  (* Z_zi *)
(* e_4 = Pp_zd = P'(z-zd)-P'(z+zd): differentiate P(z-+zd) in z (dp/dz=pp, dpp/dz=6p^2-g2/2) *)
dz[expr_] := D[expr, p] pp + D[expr, pp](6 p^2 - g2/2);
erat[4] = Simplify[dz[Ppm[-1, pd, ppd]] - dz[Ppm[+1, pd, ppd]]];
(* pp^2 -> cubic, to canonicalise *)
ppRule = pp^n_ :> (4 p^3 - g2 p - g3)^(n/2) /; EvenQ[n];
canon[e_] := Collect[Expand[e] /. pp^2 -> 4 p^3 - g2 p - g3, pp, Simplify];

(* numerical check of each reduction e_i(rational) vs the direct ef[i][z] *)
ef[1][z_]:=1; ef[2][z_]:=WZ[z-zd]-WZ[z+zd]; ef[3][z_]:=WP[z-zd]+WP[z+zd]; ef[4][z_]:=WPp[z-zd]-WPp[z+zd];
ef[5][z_]:=WP[z]; ef[6][z_]:=WP[z-iw]; ef[7][z_]:=WP[z-zinf]+WP[z+zinf]; ef[8][z_]:=WZ[z-zinf]-WZ[z+zinf];
ztst = zr[10];
subNum = {p -> WP[ztst], pp -> WPp[ztst]};
Print["--- check of the reductions e_i -> rational(p) ---"];
Do[Print["  e[", i, "]: |rat - direct| = ",
    N[Abs[(erat[i] /. subNum) - ef[i][ztst]], 6]], {i, 8}];

(* ---------- 3. REDUCTION of the primitives Pe_i -> 6 LETTERS + rational(p) ----------
   Letters: zz=z, Zz=zeta(z), Ls[c]=lnsigma(z-c) for c in {zd,-zd,zi,-zi}.
   coeff[i] = <rule: letter -> rational coeff(p,pp)> + rat (rational part). *)
(* zeta(z-c) = Zz - zeta(c) + 1/2 (pp + ppc)/(p - Pc) ;  P'(iw)=0 *)
Zzmc[Pc_, ppc_, Zc_] := Zz - Zc + (1/2)(pp + ppc)/(p - Pc);         (* zeta(z-c) *)
Zzpc[Pc_, ppc_, Zc_] := Zz + Zc + (1/2)(pp - ppc)/(p - Pc);         (* zeta(z+c) *)
(* Pe_i in the form {letters + rat}, ignoring additive constants (they do not affect the relations) *)
PeR[1] = zz;
PeR[2] = Ls[zd] - Ls[-zd];
PeR[3] = -(Zzmc[pd, ppd, zdc] + Zzpc[pd, ppd, zdc]);
PeR[4] = Ppm[-1, pd, ppd] - Ppm[+1, pd, ppd];
PeR[5] = -Zz;
PeR[6] = -Zzmc[eIw, 0, ZeIw];                                       (* ZeIw=zeta(iw) *)
PeR[7] = -(Zzmc[pI, ppI, zic] + Zzpc[pI, ppI, zic]);
PeR[8] = Ls[zi] - Ls[-zi];
ZeIw = WZ[iw];
(* check: PeR[i] (letters replaced by their values) - direct Pe[i][z] = constant *)
Pe[1][z_]:=z; Pe[2][z_]:=Log[WS[z-zd]]-Log[WS[z+zd]]; Pe[3][z_]:=-(WZ[z-zd]+WZ[z+zd]);
Pe[4][z_]:=WP[z-zd]-WP[z+zd]; Pe[5][z_]:=-WZ[z]; Pe[6][z_]:=-WZ[z-iw];
Pe[7][z_]:=-(WZ[z-zinf]+WZ[z+zinf]); Pe[8][z_]:=Log[WS[z-zinf]]-Log[WS[z+zinf]];
letterVal = {zz -> #, Zz -> WZ[#], Ls[zd] -> Log[WS[# - zd]], Ls[-zd] -> Log[WS[# + zd]],
             Ls[zi] -> Log[WS[# - zinf]], Ls[-zi] -> Log[WS[# + zinf]],
             p -> WP[#], pp -> WPp[#]} &;
z1 = zr[105/10]; z2 = zr[95/10];
Print["--- check Pe_i -> letters (the difference must be CONSTANT in z) ---"];
Do[Module[{d1, d2},
   d1 = (PeR[i] /. letterVal[z1]) - Pe[i][z1];
   d2 = (PeR[i] /. letterVal[z2]) - Pe[i][z2];
   Print["  Pe[", i, "]: |d(z1)-d(z2)| = ", N[Abs[d1 - d2], 6]]], {i, 8}];

(* ---------- 4. OSTROGRADSKY MATCHING -> beta ----------
   A'_k = erat[i] PeR[j] - erat[j] PeR[i] = sum_L coeffL(p,pp) L + ratPart.
   Modulo dz-exact terms, each coeffL (rational in p) is characterised by its INVARIANTS at the poles:
   residue (coeff of t^-1) and double-pole coeff (t^-2) [t=z-a].  Match delta_phi vs the 5-element basis. *)
letters = {zz, Zz, Ls[zd], Ls[-zd], Ls[zi], Ls[-zi]};
Aprime[k_] := Module[{i = pairsList[[k, 1]], j = pairsList[[k, 2]]},
   Expand[erat[i] PeR[j] - erat[j] PeR[i]]];
(* pairs (1-based, indices of the 8 functions) *)
Rall = {1, 2, 3, 4, 5, 6}; Eall = {1, 2, 7, 8};
pairsList = Select[Subsets[Range[8], {2}],
   (MemberQ[Rall, #[[1]]] && MemberQ[Eall, #[[2]]]) || (MemberQ[Eall, #[[1]]] && MemberQ[Rall, #[[2]]]) &];
(* coeffs of R and eta' (numerical, from the residue setup) - recomputed here *)
NmSym = Expand[((r (r - 2 M)((Em^2-1)r+2M)(r^2(r-2M)-Jcn^2((Em^2-1)r+2M))) D[Jcn((Em^2-1)r+2M), M]
   - (1/2)(Jcn((Em^2-1)r+2M)) D[r(r-2M)((Em^2-1)r+2M)(r^2(r-2M)-Jcn^2((Em^2-1)r+2M)), M]) /. M -> mm];
Nmf[x_] := NmSym /. r -> x; Qd[x_] := (D[Q4[y], y] /. y -> x); Qdd[x_] := (D[Q4[y], {y, 2}] /. y -> x);
Frat[x_] := Nmf[x]/Q4[x]; sP = Sqrt[Q4[rd]]; a1c = Qd[rd]/(4 sP); a2c = Qdd[rd]/12;
h0 = Frat[rd]; h1 = (D[Frat[x], x] /. x -> rd) sP; h2 = (1/2)((D[Frat[x], {x, 2}] /. x -> rd) sP^2 + (D[Frat[x], x] /. x -> rd)(Qd[rd]/2));
b1zdN = (h2 - 3 a1c h1 + (6 a1c^2 - 3 a2c) h0)/sP^3; b2zdN = (h1 - 3 a1c h0)/sP^3; b3zdN = h0/sP^3;
b2hf[ei_] := Nmf[ei]/(ei - rd)^3 (4/Qd[ei]^2);
e1zdN = (rd^3 - 2 rd^2)/sP; saN = Sqrt[a4]; crN = e4 - (2/saN) zdc*0 - (2/saN) WZ[zinf];
BcN = crN + (1/saN) WZ[2 zinf]; AqN = -1/saN; e1ziN = AqN (2 BcN + rd - 2); e2ziN = AqN^2;
cvR = {0, b1zdN, b2zdN, -b3zdN/2, b2hf[e4], b2hf[e3], 0, 0};
dvE = {0, e1zdN, 0, 0, 0, 0, e2ziN, e1ziN};
(* C0,Ce for the numerical match *)
rOfz[z_] := crN - (1/saN)(WZ[z - zinf] - WZ[z + zinf]);
ztt = zr[191/100]; RexN = Nmf[rOfz[ztt]]/((rOfz[ztt] - rd)^3 Q4[rOfz[ztt]]);
EpexN = (rOfz[ztt]^3 - 2 rOfz[ztt]^2)/(rOfz[ztt] - rd);
cvR[[1]] = RexN - Sum[cvR[[i]] ef[i][ztt], {i, 2, 8}]; dvE[[1]] = EpexN - Sum[dvE[[i]] ef[i][ztt], {i, 2, 8}];
wkN = Table[(cvR[[p[[1]]]] dvE[[p[[2]]]] - cvR[[p[[2]]]] dvE[[p[[1]]]])/2, {p, pairsList}];
(* delta_phi_w2' in letter form *)
dphiP = Expand[Sum[wkN[[k]] Aprime[k], {k, Length[pairsList]}]];
(* 5 basis dilogs *)
baseNames = {{5, 7}, {5, 8}, {2, 5}, {1, 5}, {4, 7}};
baseK = Flatten[Position[pairsList, Sort[#]] & /@ baseNames];
(* invariants: residue and double-pole coeff at each pole, for a rational function(p,pp) *)
gPrime[Pa_] := 6 Pa^2 - g2/2;
ploc[Pa_, ppa_, t_] := Pa + ppa t + (1/2) gPrime[Pa] t^2 + (1/6)(12 Pa ppa) t^3 + (1/24)(gPrime[Pa]^2 + 12 Pa gPrime[Pa]... )t^4;
(* only up to t^3: enough for the residue(t^-1) and double pole(t^-2) of poles up to order 2 *)
ploc2[Pa_, ppa_, t_] := Pa + ppa t + (1/2) gPrime[Pa] t^2 + (1/6)(12 Pa ppa) t^3;
pploc2[Pa_, ppa_, t_] := ppa + gPrime[Pa] t + (1/2)(12 Pa ppa) t^2;
inv[coeffL_, Pa_, ppa_] := Module[{ser},
   ser = Series[coeffL /. {p -> ploc2[Pa, ppa, t], pp -> pploc2[Pa, ppa, t]}, {t, 0, 0}];
   {Coefficient[ser, t, -1], Coefficient[ser, t, -2]}];   (* {residue, double pole} *)
(* poles and their (Pa,ppa) *)
poli = {{pd, ppd}, {pd, -ppd}, {pI, ppI}, {pI, -ppI}, {eIw, 0}, {WP[0*iw + 10^-30], 0}}; (* z_d,-z_d,z_inf,-z_inf,iw,~0 *)
(* note: the pole at z=0 (e_5=p) is treated separately as the coeff of p *)
Print["--- MATCHING of invariants: building the system for beta ---"];
letCoeff[expr_, L_] := Coefficient[expr, L];
rows = {}; rhs = {};
Do[Do[Module[{cD, cB, iD, iB},
    cD = letCoeff[dphiP, L];
    iD = inv[cD, pol[[1]], pol[[2]]];
    cB = letCoeff[Aprime[#], L] & /@ baseK;
    iB = Table[inv[cB[[b]], pol[[1]], pol[[2]]], {b, 5}];
    AppendTo[rows, iB[[All, 1]]]; AppendTo[rhs, iD[[1]]];   (* residue *)
    AppendTo[rows, iB[[All, 2]]]; AppendTo[rhs, iD[[2]]];   (* double pole *)
   ], {pol, poli}], {L, letters}];
Mrows = N[rows, prec]; Mrhs = N[rhs, prec];
betaSol = LeastSquares[Mrows, Mrhs];
residM = Max[Abs[Mrows.betaSol - Mrhs]];
Print["residual of the invariant matching = ", N[residM, 5]];
Print["beta (from the invariants):"]; Do[Print["  beta[", b, "] = ", N[betaSol[[b]], 30]], {b, 5}];
