#!/usr/bin/env wolframscript
(* ::Package:: *)

(* verify_two_branch_optics.wls
   -------------------------------------------------------------------------
   WHICH optical surface, and for WHICH clock.

   This script exists because an earlier draft attributed a coordinate-time
   result to the proper-time branch.  The two branches carry DIFFERENT optical
   metrics and, as it turns out, OPPOSITE global behaviour, so the distinction
   is not bookkeeping.

   Static case g = -f dt^2 + h, rail -g(u,d/dt) = Ehat.  Then f u^t = Ehat and
   |v|_h^2 = (Ehat^2 - f)/f, whence per unit h-length:

       dt   = Lambda_t   |dx|_h ,   Lambda_t^2   = Ehat^2/(f (Ehat^2 - f))
       dtau = Lambda_tau |dx|_h ,   Lambda_tau^2 = f/(Ehat^2 - f)

   Lambda_t is Perlick's fixed-energy optical factor: it measures ARRIVAL
   (coordinate) time.  Lambda_tau measures PROPER time.  They differ by
   Ehat^2/f^2, which is not a constant, so the two surfaces are not homothetic
   and their curvatures are unrelated.
*)

$Failures = 0;
check[l_, c_] := (Print[If[TrueQ[c], "PASS  ", "FAIL  "], l];
   If[! TrueQ[c], $Failures++]);

Print["=========================================================="];
Print[" the two optical surfaces, and which clock each measures"];
Print["=========================================================="];

ClearAll[r, ph, M, Eh, f];
ff[r_] := 1 - 2 M/r;

(* --- which clock ------------------------------------------------------- *)
ut = Eh/f; vh = Sqrt[(Eh^2 - f)/f];
check["Lambda_t^2 = Ehat^2/(f(Ehat^2-f)) is dt per unit h-length (Perlick)",
      Simplify[(ut/vh)^2 - Eh^2/(f (Eh^2 - f))] === 0];
check["Lambda_tau^2 = f/(Ehat^2-f) is dtau per unit h-length",
      Simplify[(1/vh)^2 - f/(Eh^2 - f)] === 0];
check["they differ by Ehat^2/f^2, not by a constant: not homothetic",
      Simplify[(Eh^2/(f (Eh^2 - f)))/(f/(Eh^2 - f)) - Eh^2/f^2] === 0];

(* --- Gauss curvature of each ------------------------------------------- *)
gauss[Lam2_] := Module[{aa, ai, Ch, Rm, Rc, co},
   co = {r, ph};
   aa = Lam2 {{1/ff[r], 0}, {0, r^2}};
   ai = Simplify[Inverse[aa]];
   Ch = Table[Simplify[(1/2) Sum[ai[[i, l]] (D[aa[[l, j]], co[[k]]] +
        D[aa[[l, k]], co[[j]]] - D[aa[[j, k]], co[[l]]]), {l, 2}]],
      {i, 2}, {j, 2}, {k, 2}];
   Rm = Table[Simplify[D[Ch[[i, j, l]], co[[k]]] - D[Ch[[i, j, k]], co[[l]]] +
       Sum[Ch[[i, k, s]] Ch[[s, j, l]] - Ch[[i, l, s]] Ch[[s, j, k]], {s, 2}]],
      {i, 2}, {j, 2}, {k, 2}, {l, 2}];
   Rc = Table[Simplify[Sum[Rm[[s, i, s, j]], {s, 2}]], {i, 2}, {j, 2}];
   Simplify[(1/2) Sum[ai[[i, j]] Rc[[i, j]], {i, 2}, {j, 2}]]];

Kt = gauss[Eh^2/(ff[r] (Eh^2 - ff[r]))];
Ktau = gauss[ff[r]/(Eh^2 - ff[r])];
Print["----------------------------------------------------------"];
Print["   K_t   = ", Kt];
Print["   K_tau = ", Ktau];
check["the two curvatures are different functions", Simplify[Kt - Ktau] =!= 0];

(* --- far field --------------------------------------------------------- *)
at = Simplify[Limit[r^3 Kt, r -> Infinity]];
atau = Simplify[Limit[r^3 Ktau, r -> Infinity]];
Print["   r^3 K_t   -> ", at];
Print["   r^3 K_tau -> ", atau];
check["arrival branch: r^3 K_t -> -M(2Ehat^2-3)/Ehat^2, sign set by Ehat",
      Simplify[at + M (2 Eh^2 - 3)/Eh^2] === 0];
check["proper-time branch: r^3 K_tau -> +M, POSITIVE for every Ehat",
      Simplify[atau - M] === 0];
Print["   => the far field focuses on the proper-time branch for EVERY rail"];
Print["      energy, and defocuses on the arrival branch above Ehat^2 = 3/2."];

(* --- sign over the exterior -------------------------------------------- *)
Print["----------------------------------------------------------"];
sgn[K_, lab_] := Module[{neg, tot},
   neg = Select[Flatten[Table[{rv, ev, N[K /. {M -> 1, r -> rv, Eh -> ev}, 20]},
        {rv, 21/10, 100, 1/2}, {ev, 101/100, 20, 1/20}], 1], #[[3]] < 0 &];
   tot = Length[Flatten[Table[{rv, ev}, {rv, 21/10, 100, 1/2},
        {ev, 101/100, 20, 1/20}], 1]];
   Print["   ", lab, ": ", Length[neg], " negative samples out of ", tot];
   Length[neg]/tot // N];
fT = sgn[Kt, "K_t  "]; fTau = sgn[Ktau, "K_tau"];
check["K_tau is positive on essentially the whole exterior", fTau < 1/1000];

(* --- completeness ------------------------------------------------------ *)
Print["----------------------------------------------------------"];
dsT[r_] := Sqrt[Eh^2/(ff[r] (Eh^2 - ff[r]))]/Sqrt[ff[r]];
dsTau[r_] := Sqrt[ff[r]/(Eh^2 - ff[r])]/Sqrt[ff[r]];
poleT = Simplify[Limit[(r - 2 M) dsT[r], r -> 2 M, Direction -> "FromAbove"],
   Assumptions -> {M > 0, Eh > 1}];
limTau = Simplify[Limit[dsTau[r], r -> 2 M, Direction -> "FromAbove"],
   Assumptions -> {M > 0, Eh > 1}];
Print["   inner end:  (r-2M) ds_t/dr -> ", poleT,
      " (simple pole, log-divergent)"];
Print["               ds_tau/dr      -> ", limTau, " (finite)"];
check["arrival branch: simple pole at r=2M with residue 2M => complete",
      Simplify[poleT - 2 M] === 0];
check["proper-time branch: finite limit 1/Ehat => FINITE distance, incomplete",
      Simplify[limTau - 1/Eh] === 0];
num = {M -> 1, Eh -> 14/10};
dT[e_] := NIntegrate[Evaluate[dsT[r] /. num], {r, 2 + e, 3}, WorkingPrecision -> 30];
dTau[e_] := NIntegrate[Evaluate[dsTau[r] /. num], {r, 2 + e, 3}, WorkingPrecision -> 30];
Print["   distance from r=3M to r=2M+eps  (M=1, Ehat=1.4):"];
Do[Print["     eps=", ScientificForm[N[e, 2]], "    t: ", N[dT[e], 8],
    "    tau: ", N[dTau[e], 8]], {e, {1/10^2, 1/10^4, 1/10^6, 1/10^8}}];
check["the proper-time distance to r=2M converges (finite, ~0.7523)",
      Abs[N[dTau[1/10^8] - dTau[1/10^6], 30]] < 10^-5];
check["the arrival-time distance diverges (grows by ~9.21 per two decades)",
      N[dT[1/10^8] - dT[1/10^6], 30] > 9];

Print["----------------------------------------------------------"];
Print["   CONCLUSION.  Cartan-Hadamard applies to the ARRIVAL-time surface"];
Print["   (complete, K < 0 for Ehat^2 >= 3/2) and to that one only.  The"];
Print["   proper-time surface is positively curved and incomplete at r = 2M:"];
Print["   it focuses, and a falling worldline reaches the horizon in finite"];
Print["   proper time, as it must.  Any no-conjugate-point or global"];
Print["   minimality statement therefore belongs to the arrival branch."];

Print["=========================================================="];
Print[If[$Failures === 0, "ALL CHECKS PASSED",
     ToString[$Failures] <> " CHECK(S) FAILED"]];
Print["=========================================================="];
If[$Failures > 0, Exit[1], Exit[0]];
