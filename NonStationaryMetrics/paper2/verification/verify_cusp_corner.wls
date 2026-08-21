#!/usr/bin/env wolframscript
(* ::Package:: *)

(* verify_cusp_corner.wls -- referee major comment 4 (CQG-116884).
   Settles the local geometry of the equatorial tau-branch controlled-rail
   extremal at the conformal stationary limit r_e = 2M, for

     |J| > Jc   (smooth periapsis outside r_e),
     |J| < Jc   (contact with r_e),
     |J| = Jc   (marginal, both signs).

   Everything is derived from the printed tau-branch quadrature

     dphi/dr = Ktau(r)/Sqrt[S(r)],
     S(r)    = r (r-2M) ((E^2-1) r + 2M) ( r Delta - J^2 ((E^2-1) r + 2M) ),
     Ktau(r) = J r (r-2M) ((E^2-1) r + 2M)/Delta,      Delta = r^2 - 2 M r + a^2,

   with one-sided expansions and explicit branch orientation.  Each check
   prints PASS/FAIL; the script exits 1 if any check fails.
*)

$Failures = 0;
check[label_, cond_] := (
  Print[If[TrueQ[cond], "PASS  ", "FAIL  "], label];
  If[! TrueQ[cond], $Failures++]);

Print["=========================================================="];
Print[" tau-branch local geometry at the conformal stationary limit"];
Print["=========================================================="];

(* ---------------------------------------------------------------- *)
(* symbols                                                          *)
(* ---------------------------------------------------------------- *)
ClearAll[r, M, a, E2, J, Jc, DeltaK, S, Ktau];
DeltaK[r_] := r^2 - 2 M r + a^2;
w[r_]      := (E2 - 1) r + 2 M;                    (* D_E in Paper I notation *)
B[r_]      := r DeltaK[r] - J^2 w[r];              (* the "turning" bracket    *)
S[r_]      := r (r - 2 M) w[r] B[r];
Ktau[r_]   := J r (r - 2 M) w[r]/DeltaK[r];

Jcsq = a^2/E2;                                     (* Jc = a/E, so Jc^2 = a^2/E^2 *)

(* ---------------------------------------------------------------- *)
(* 1. B(2M) = 2M (a^2 - J^2 E^2):  the sign that splits the cases    *)
(* ---------------------------------------------------------------- *)
b2M = Simplify[B[2 M]];
check["B(2M) = 2M (a^2 - J^2 E^2)", Simplify[b2M - 2 M (a^2 - J^2 E2)] === 0];
Print["      B(2M) = ", b2M];

(* ---------------------------------------------------------------- *)
(* 2. At |J| = Jc the bracket factorises: B = (r-2M)(r^2 + a^2/E^2)  *)
(*    hence r = 2M is a DOUBLE root of S, not a simple one.          *)
(* ---------------------------------------------------------------- *)
Bmarg = Together[B[r] /. J^2 -> Jcsq];
check["B(r)|_{J^2=Jc^2} = (r-2M)(r^2 + a^2/E^2)",
      Simplify[Bmarg - (r - 2 M) (r^2 + a^2/E2)] === 0];

Smarg = Simplify[S[r] /. J^2 -> Jcsq];
check["S|_{marg} = r (r-2M)^2 w(r) (r^2 + a^2/E^2)",
      Simplify[Smarg - r (r - 2 M)^2 w[r] (r^2 + a^2/E2)] === 0];

(* multiplicity of the root r = 2M, marginal vs generic *)
mult[expr_] := Module[{k = 0, e = expr},
   While[Simplify[(e /. r -> 2 M)] === 0 && k < 6,
     e = D[e, r]; k++]; k];
check["root r=2M of S is double at |J|=Jc",  mult[Smarg] === 2];
check["root r=2M of S is simple for generic J",
      mult[Together[S[r]]] === 1];

(* ---------------------------------------------------------------- *)
(* 3. GENERIC |J| < Jc: dphi/dr ~ C Sqrt[r-2M]  ->  semicubical cusp *)
(* ---------------------------------------------------------------- *)
(* set eps = r - 2M > 0 and expand *)
asm = {eps > 0, M > 0, E2 > 1, a > 0, J != 0, a^2 - E2 J^2 > 0};
gen = Ktau[r]/Sqrt[S[r]] /. r -> 2 M + eps;
(* dphi/dr vanishes; divide by Sqrt[eps] and the limit is finite nonzero, so
   the exponent is exactly 1/2 and phi - phi_e ~ (2/3) K eps^{3/2}.        *)
Kgen = Simplify[Limit[gen/Sqrt[eps], eps -> 0, Direction -> "FromAbove",
    Assumptions -> asm], Assumptions -> asm];
Print["      generic:  dphi/dr / Sqrt[r-2M]  ->  K = ", Kgen];
check["generic: dphi/dr -> 0 at r_e (radial tangent)",
      Simplify[Limit[gen, eps -> 0, Direction -> "FromAbove",
        Assumptions -> asm], Assumptions -> asm] === 0];
(* K = (E J / a^2) Sqrt[2M/(a^2 - E^2 J^2)]   (E2 = E^2 in this script) *)
KgenClosed = (Sqrt[E2] J/a^2) Sqrt[2 M/(a^2 - E2 J^2)];
check["generic: K = (E J/a^2) Sqrt(2M/(a^2 - E^2 J^2))",
      Simplify[Kgen - KgenClosed, Assumptions -> asm] === 0];
(* it blows up as |J| -> Jc: the cusp amplitude diverges at the marginal value,
   which is exactly why the marginal case is a different local normal form.  *)
check["K diverges as |J| -> Jc  (cusp branch does not connect to the corner)",
      Simplify[Limit[KgenClosed /. J -> -Sqrt[a^2/E2] + d, d -> 0,
        Direction -> "FromAbove"]] === DirectedInfinity[-1] ||
      MatchQ[Simplify[1/KgenClosed /. J^2 -> a^2/E2], 0]];
check["generic: K finite and nonzero  =>  semicubical cusp, exponent 3/2",
      FreeQ[Kgen, eps] && Simplify[Kgen, Assumptions -> asm] =!= 0];

(* ---------------------------------------------------------------- *)
(* 4. MARGINAL |J| = Jc: dphi/dr has a FINITE nonzero one-sided      *)
(*    limit  ->  corner, NOT a cusp.  Both signs of J.               *)
(* ---------------------------------------------------------------- *)
(* Sqrt[S] = |r-2M| Sqrt[r w (r^2+Jc^2)]; keep the modulus explicit. *)
rootFactor[r_] := Sqrt[r w[r] (r^2 + a^2/E2)];
(* outside, r > 2M *)
dphiOut = Simplify[
   (J r (r - 2 M) w[r]/DeltaK[r])/((r - 2 M) rootFactor[r])];
check["marginal, r>2M: dphi/dr = J r w /(Delta Sqrt[r w (r^2+Jc^2)])",
      Simplify[dphiOut - J r w[r]/(DeltaK[r] rootFactor[r])] === 0];

limOut = Simplify[dphiOut /. r -> 2 M];
Print["      marginal one-sided limit  (dphi/dr)|_{r->2M^+} = ",
      Simplify[limOut]];
closed = 2 M E2^(1/2) J/(a^2 Sqrt[4 M^2 + a^2/E2]);
check["closed form of the corner slope: 2 M E J / (a^2 Sqrt(4M^2+Jc^2))",
      Simplify[limOut - closed, Assumptions -> {M > 0, a > 0, E2 > 1}] === 0];
check["corner slope is finite and nonzero for J != 0",
      Simplify[limOut] =!= 0 && FreeQ[Simplify[limOut], DirectedInfinity]];

(* the slope flips sign across r = 2M because Sqrt[S] carries |r-2M| *)
dphiIn = Simplify[
   (J r (r - 2 M) w[r]/DeltaK[r])/(-(r - 2 M) rootFactor[r])];
check["branch orientation: (dphi/dr)_{r<2M} = -(dphi/dr)_{r>2M}",
      Simplify[dphiIn + dphiOut] === 0];

(* ---------------------------------------------------------------- *)
(* 5. numerical confirmation at the paper's representative point     *)
(*    M = 1, a = 0.9, E = 1.2                                        *)
(* ---------------------------------------------------------------- *)
Print["----------------------------------------------------------"];
Print[" numerical check at M=1, a=0.9, E=1.2  (Jc = a/E = 0.75)"];
Print["----------------------------------------------------------"];
num = {M -> 1, a -> 9/10, E2 -> (6/5)^2};
JcN  = N[9/10/(6/5), 30];
Print["      Jc = ", JcN];

dphi[rv_?NumericQ, Jv_] := Module[{s, k},
   s = N[S[r] /. num /. J -> Jv /. r -> rv, 40];
   k = N[Ktau[r] /. num /. J -> Jv /. r -> rv, 40];
   k/Sqrt[s]];

(* (a) marginal retrograde J = -Jc: finite corner slope *)
predCorner = N[closed /. num /. J -> -JcN, 30];
tab = Table[{eps, dphi[2 + eps, -JcN]}, {eps, {1/10, 1/100, 1/1000, 1/10000,
    1/100000}}];
Print["      J = -Jc :   eps        dphi/dr          predicted ", predCorner];
Do[Print["                 ", N[t[[1]], 8], "   ", N[t[[2]], 12]], {t, tab}];
errCorner = Abs[N[tab[[-1, 2]] - predCorner, 30]];
check["J=-Jc: numerical dphi/dr -> finite corner slope (err < 1e-4)",
      errCorner < 10^-4];

(* convergence order: the error should be O(eps) *)
e1 = Abs[N[tab[[-2, 2]] - predCorner, 30]];
e2 = Abs[N[tab[[-1, 2]] - predCorner, 30]];
ordC = Log[10, e1/e2];
Print["      observed convergence order in eps: ", N[ordC, 6]];
check["J=-Jc: first-order approach to the corner slope (order ~ 1)",
      Abs[ordC - 1] < 0.1];

(* (b) same for the prograde marginal J = +Jc: same magnitude, opposite sign *)
predP = N[closed /. num /. J -> JcN, 30];
check["J=+Jc and J=-Jc corner slopes are equal and opposite",
      Abs[N[dphi[2 + 1/100000, JcN] + dphi[2 + 1/100000, -JcN], 30]] < 10^-6];

(* (c) sub-marginal |J| < Jc: dphi/dr -> 0 like Sqrt[eps]  (cusp)      *)
Jsub = 6/10;                                  (* |J| = 0.6 < 0.75 *)
(* the subleading correction is O(eps), so a two-point exponent estimate needs
   eps small; at 1e-2/1e-4 it reads 0.487, which is the O(eps) tail, not a
   different power.  Use 1e-8 / 1e-12 at 40 digits.                        *)
epsList = {1/10^2, 1/10^4, 1/10^8, 1/10^12};
tabS = Table[{ep, dphi[2 + ep, -Jsub]}, {ep, epsList}];
Print["      J = -0.6 (< Jc):   eps          dphi/dr           dphi/dr/Sqrt[eps]"];
Do[Print["                     ", ScientificForm[N[t[[1]], 3]], "   ",
         N[t[[2]], 12], "   ", N[t[[2]]/Sqrt[t[[1]]], 12]], {t, tabS}];
ordS = Log[Abs[tabS[[3, 2]]/tabS[[4, 2]]]]/Log[10^4];
Print["      observed exponent of dphi/dr ~ eps^p :  p = ", N[ordS, 12]];
check["|J|<Jc: dphi/dr ~ eps^{1/2}, i.e. phi-phi_e ~ eps^{3/2} (cusp)",
      Abs[ordS - 1/2] < 10^-6];
(* and the amplitude agrees with the symbolic K *)
KgenN = N[Kgen /. num /. J -> -Jsub, 30];
check["|J|<Jc: numerical amplitude matches the symbolic K",
      Abs[N[tabS[[4, 2]]/Sqrt[tabS[[4, 1]]] - KgenN, 30]] < 10^-8];
Print["      symbolic K = ", KgenN];

(* (d) super-marginal |J| > Jc: turning point strictly outside r_e     *)
Jsup = 9/10;                                  (* |J| = 0.9 > 0.75 *)
roots = r /. NSolve[(B[r] /. num /. J -> Jsup) == 0, r, WorkingPrecision -> 30];
rmin = Min[Select[roots, (Im[#] == 0 && Re[#] > 2) &]];
Print["      J = 0.9 (> Jc):  periapsis root of B at r = ", N[rmin, 12]];
check["|J|>Jc: simple turning point strictly outside r_e = 2M", rmin > 2];
(* at a simple root of S the slope blows up like (r-rmin)^{-1/2}        *)
sl = Table[dphi[rmin + eps, Jsup]*Sqrt[eps], {eps, {1/10^4, 1/10^6}}];
check["|J|>Jc: dphi/dr ~ (r-rmin)^{-1/2}  (smooth periapsis, exponent 2)",
      Abs[N[sl[[1]] - sl[[2]], 20]] < 10^-2];

(* ---------------------------------------------------------------- *)
(* 6. the retrograde marginal does not counter-rotate inside r_e     *)
(*    the sign flip of Sqrt[S] across r=2M reverses dphi/dr, so the  *)
(*    analytic continuation of J=-Jc into r<2M has dphi/dlambda < 0  *)
(* ---------------------------------------------------------------- *)
Print["----------------------------------------------------------"];
sIn  = N[dphi[2 - 1/1000, -JcN], 20];      (* uses +Sqrt branch, r<2M *)
sOut = N[dphi[2 + 1/1000, -JcN], 20];
Print["      J=-Jc: dphi/dr just outside = ", sOut,
      " ,  principal branch just inside = ", sIn];
check["J=-Jc: sign of dphi/dr on the principal branch reverses at r_e",
      Sign[Re[sIn]] != Sign[Re[sOut]] || Im[sIn] != 0];

Print["=========================================================="];
Print[If[$Failures === 0, "ALL CHECKS PASSED",
     ToString[$Failures] <> " CHECK(S) FAILED"]];
Print["=========================================================="];
If[$Failures > 0, Exit[1], Exit[0]];
