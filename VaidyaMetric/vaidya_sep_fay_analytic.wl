(* ::Package:: *)
(* =====================================================================
   FAY ANALITICA - riduzione simbolica dei dilog e beta in forma chiusa.
   ---------------------------------------------------------------------
   IDEA: con i teoremi di addizione di Weierstrass, ogni funzione pari e_i(z)
   e' RAZIONALE in p=P(z) (+ costanti ai poli). Le primitive Pe_i si riducono
   a 6 LETTERE trascendenti {z, Z(z)=zeta(z), Ls[+-zd], Ls[+-zi]}.
   Le relazioni di Fay = MATCHING dei coefficienti (razionali in p,pp) davanti
   a ciascuna lettera  ->  algebra ESATTA (niente sampling/condizionamento).
   Costanti ai poli tenute SIMBOLICHE  ->  beta SIMBOLICI (poi verifica numerica).
   Convenzioni identiche a vaidya_sep_5term.py / vaidya_sep_fay_symbolic.wl.
   ===================================================================== *)

$MinPrecision = 60; prec = 55;
Off[NIntegrate::ncvb]; Off[NIntegrate::slwcon]; Off[NIntegrate::nlim]; Off[General::stop];

(* ---------- 0-1. curva + reticolo (numerico high-prec, riuso setup) ---------- *)
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

(* costanti ai poli (valori NUMERICI; le terremo simboliche nel matching) *)
pd = WP[zd]; ppd = WPp[zd]; zdc = WZ[zd];       (* P,P',zeta a z_d *)
pI = WP[zinf]; ppI = WPp[zinf]; zic = WZ[zinf]; (* a z_inf *)
eIw = WP[iw];                                    (* P a semiperiodo iw (P'(iw)=0) *)
Print["P(iw) (2-torsione, P'=0): ", N[eIw,10], "  check P'(iw)=", N[WPp[iw],6]];

(* ---------- 2. RIDUZIONE  e_i(z) -> razionale in p=P(z), pp=P'(z)  (addizione) ---------- *)
(* addizione: P(z-+c) = -P(z)-P(c) + 1/4 ((pp -+ ppc)/(p-Pc))^2 ; con p,pp simboli *)
Ppm[sgn_, Pc_, ppc_] := -p - Pc + (1/4)((pp - sgn ppc)/(p - Pc))^2;      (* P(z + sgn c): addizione con P'(sgn c)=sgn ppc *)
(* Z_a = -2 zeta(a) + P'(a)/(p - P(a))   [derivato a mano, verificato sotto] *)
Zred[Za_, Pa_, ppa_] := -2 Za + ppa/(p - Pa);
(* forme razionali (in p,pp) delle 8 funzioni *)
erat[1] = 1;
erat[2] = Zred[zdc, pd, ppd];                                   (* Z_zd *)
erat[3] = Ppm[-1, pd, ppd] + Ppm[+1, pd, ppd] // Simplify;       (* P_zd = P(z-zd)+P(z+zd) *)
erat[5] = p;                                                    (* wp0 = P(z) *)
erat[6] = -p - eIw + (1/4)(pp/(p - eIw))^2 // Simplify;          (* P(z-iw), P'(iw)=0 *)
erat[7] = Ppm[-1, pI, ppI] + Ppm[+1, pI, ppI] // Simplify;       (* P_zi *)
erat[8] = Zred[zic, pI, ppI];                                  (* Z_zi *)
(* e_4 = Pp_zd = P'(z-zd)-P'(z+zd): derivo P(z-+zd) in z (dp/dz=pp, dpp/dz=6p^2-g2/2) *)
dz[expr_] := D[expr, p] pp + D[expr, pp](6 p^2 - g2/2);
erat[4] = Simplify[dz[Ppm[-1, pd, ppd]] - dz[Ppm[+1, pd, ppd]]];
(* pp^2 -> cubica per canonizzare *)
ppRule = pp^n_ :> (4 p^3 - g2 p - g3)^(n/2) /; EvenQ[n];
canon[e_] := Collect[Expand[e] /. pp^2 -> 4 p^3 - g2 p - g3, pp, Simplify];

(* verifica numerica di ogni riduzione e_i(razionale) vs ef[i][z] diretto *)
ef[1][z_]:=1; ef[2][z_]:=WZ[z-zd]-WZ[z+zd]; ef[3][z_]:=WP[z-zd]+WP[z+zd]; ef[4][z_]:=WPp[z-zd]-WPp[z+zd];
ef[5][z_]:=WP[z]; ef[6][z_]:=WP[z-iw]; ef[7][z_]:=WP[z-zinf]+WP[z+zinf]; ef[8][z_]:=WZ[z-zinf]-WZ[z+zinf];
ztst = zr[10];
subNum = {p -> WP[ztst], pp -> WPp[ztst]};
Print["--- verifica riduzioni e_i -> razionale(p) ---"];
Do[Print["  e[", i, "]: |rat - diretto| = ",
    N[Abs[(erat[i] /. subNum) - ef[i][ztst]], 6]], {i, 8}];

(* ---------- 3. RIDUZIONE primitive Pe_i -> 6 LETTERE + razionale(p) ----------
   Lettere: zz=z, Zz=zeta(z), Ls[c]=lnsigma(z-c) per c in {zd,-zd,zi,-zi}.
   coeff[i] = <regola: lettera -> coeff razionale(p,pp)> + rat (parte razionale). *)
(* zeta(z-c) = Zz - zeta(c) + 1/2 (pp + ppc)/(p - Pc) ;  P'(iw)=0 *)
Zzmc[Pc_, ppc_, Zc_] := Zz - Zc + (1/2)(pp + ppc)/(p - Pc);         (* zeta(z-c) *)
Zzpc[Pc_, ppc_, Zc_] := Zz + Zc + (1/2)(pp - ppc)/(p - Pc);         (* zeta(z+c) *)
(* Pe_i in forma {lettere + rat}, ignorando costanti additive (non toccano le relazioni) *)
PeR[1] = zz;
PeR[2] = Ls[zd] - Ls[-zd];
PeR[3] = -(Zzmc[pd, ppd, zdc] + Zzpc[pd, ppd, zdc]);
PeR[4] = Ppm[-1, pd, ppd] - Ppm[+1, pd, ppd];
PeR[5] = -Zz;
PeR[6] = -Zzmc[eIw, 0, ZeIw];                                       (* ZeIw=zeta(iw) *)
PeR[7] = -(Zzmc[pI, ppI, zic] + Zzpc[pI, ppI, zic]);
PeR[8] = Ls[zi] - Ls[-zi];
ZeIw = WZ[iw];
(* verifica: PeR[i] (con lettere sostituite dai valori) - Pe[i][z] diretto = costante *)
Pe[1][z_]:=z; Pe[2][z_]:=Log[WS[z-zd]]-Log[WS[z+zd]]; Pe[3][z_]:=-(WZ[z-zd]+WZ[z+zd]);
Pe[4][z_]:=WP[z-zd]-WP[z+zd]; Pe[5][z_]:=-WZ[z]; Pe[6][z_]:=-WZ[z-iw];
Pe[7][z_]:=-(WZ[z-zinf]+WZ[z+zinf]); Pe[8][z_]:=Log[WS[z-zinf]]-Log[WS[z+zinf]];
letterVal = {zz -> #, Zz -> WZ[#], Ls[zd] -> Log[WS[# - zd]], Ls[-zd] -> Log[WS[# + zd]],
             Ls[zi] -> Log[WS[# - zinf]], Ls[-zi] -> Log[WS[# + zinf]],
             p -> WP[#], pp -> WPp[#]} &;
z1 = zr[105/10]; z2 = zr[95/10];
Print["--- verifica Pe_i -> lettere (differenza deve essere COSTANTE in z) ---"];
Do[Module[{d1, d2},
   d1 = (PeR[i] /. letterVal[z1]) - Pe[i][z1];
   d2 = (PeR[i] /. letterVal[z2]) - Pe[i][z2];
   Print["  Pe[", i, "]: |d(z1)-d(z2)| = ", N[Abs[d1 - d2], 6]]], {i, 8}];
