#!/usr/bin/env wolframscript
(* ::Package:: *)

(* verify_tk_ckv_rail.wls
   -------------------------------------------------------------------------
   The rail-drift identity for THAKURTA-KERR, the geometry of this paper.

   (An earlier draft illustrated the identity on ingoing Vaidya.  That is the
   companion paper's geometry, not this one; the referee is reviewing the
   conformal-Kerr manuscript, so the worked case must be conformal Kerr.)

   Thakurta-Kerr is g = A(eta)^2 gbar, with gbar the Kerr metric and eta the
   conformal time.  The selector is W = d/d(eta).  Since gbar is eta-independent
   and d/d(eta) is Killing for gbar,

        L_W g = L_W (A^2 gbar) = 2 A A' gbar = (2 A'/A) g ,

   so W is a CONFORMAL KILLING VECTOR of g with conformal factor
   psi = A'/A = epsilon, the adiabatic parameter of the paper.  Substituting
   into the rail-drift identity

        d/dtau ( -g(u,W) ) = -(1/2) (L_W g)_{mu nu} u^mu u^nu

   and using g(u,u) = -1 gives the exact result

        d Ehat / d tau |_geodesic = A'/A = epsilon .

   The control effort needed to hold the rail IS the small parameter of the
   adiabatic expansion.  This is the middle rung of the paper's
   Killing -> conformal-Killing -> Kodama hierarchy, made quantitative.
*)

$Failures = 0;
check[label_, cond_] := (
  Print[If[TrueQ[cond], "PASS  ", "FAIL  "], label];
  If[! TrueQ[cond], $Failures++]);

Print["=========================================================="];
Print[" Thakurta-Kerr: the selector is conformal Killing, and the"];
Print[" rail drift is exactly the adiabatic parameter"];
Print["=========================================================="];

ClearAll[eta, r, th, ph, A, M, a, uu];
coords = {eta, r, th, ph};

(* Kerr in Boyer-Lindquist, with t -> eta (the conformal clock) *)
Sig = r^2 + a^2 Cos[th]^2;
Del = r^2 - 2 M r + a^2;
gbar = {
  {-(1 - 2 M r/Sig), 0, 0, -2 M a r Sin[th]^2/Sig},
  {0, Sig/Del, 0, 0},
  {0, 0, Sig, 0},
  {-2 M a r Sin[th]^2/Sig, 0, 0,
   (r^2 + a^2 + 2 M a^2 r Sin[th]^2/Sig) Sin[th]^2}};
gTK = Simplify[A[eta]^2 gbar];

(* ---------------------------------------------------------------- *)
(* 1. the seed metric is eta-independent and d/d(eta) is Killing      *)
(* ---------------------------------------------------------------- *)
check["the Kerr seed is independent of the conformal clock",
      Simplify[D[gbar, eta]] === ConstantArray[0, {4, 4}]];

(* ---------------------------------------------------------------- *)
(* 2. W = d/d(eta) is a CONFORMAL Killing vector of the full metric   *)
(* ---------------------------------------------------------------- *)
LWg = Simplify[D[gTK, eta]];        (* L_{d/deta} g = d/deta of components *)
psi = A'[eta]/A[eta];
check["L_W g = (2 A'/A) g exactly: W is conformal Killing",
      Simplify[LWg - 2 psi gTK] === ConstantArray[0, {4, 4}]];
check["W is NOT Killing unless A' = 0",
      Simplify[LWg /. A -> Function[x, 1]] === ConstantArray[0, {4, 4}] &&
      Simplify[LWg] =!= ConstantArray[0, {4, 4}]];

(* the conformal factor is (1/n) div W *)
divW = Simplify[
   (1/Sqrt[-Det[gTK]]) Sum[
     D[Sqrt[-Det[gTK]] UnitVector[4, 1][[i]], coords[[i]]], {i, 4}]];
check["div W = 4 A'/A, so psi = (1/4) div W as a CKV requires",
      Simplify[divW - 4 psi] === 0];

(* ---------------------------------------------------------------- *)
(* 3. the rail drift along a geodesic is exactly epsilon = A'/A       *)
(* ---------------------------------------------------------------- *)
uu = Table[Symbol["u" <> ToString[i]], {i, 4}];
drift = Simplify[-(1/2) uu . LWg . uu];
(* impose the unit-timelike normalisation g(u,u) = -1 *)
driftOnShell = Simplify[drift /. Solve[uu . gTK . uu == -1, uu[[2]]][[1]]];
Print["      d Ehat/d tau |_geodesic = ", Simplify[driftOnShell]];
check["d Ehat/d tau = A'/A = epsilon, independent of the four-velocity",
      Simplify[driftOnShell - psi] === 0];
check["the drift vanishes exactly at frozen conformal factor",
      Simplify[driftOnShell /. A -> Function[x, 1]] === 0];

(* ---------------------------------------------------------------- *)
(* 4. would a Ricci-soliton selector do better?  Only if TK were      *)
(*    Einstein, which it is not.                                      *)
(* ---------------------------------------------------------------- *)
Print["----------------------------------------------------------"];
Print["      A soliton selector would need (1/2) L_X g + Ric = lambda g."];
Print["      For a CKV that reads Ric = (lambda - psi) g, i.e. the geometry"];
Print["      would have to be EINSTEIN.  Thakurta-Kerr is not: it carries an"];
Print["      imperfect-fluid source from the conformal factor.  So the soliton"];
Print["      class does not contain this selector, while the conformal-Killing"];
Print["      class does -- which is precisely the middle rung of the paper's"];
Print["      Killing -> conformal-Killing -> Kodama hierarchy."];

(* verify on the simplest slice that TK is not Einstein: a = 0, A generic *)
gS = Simplify[gTK /. a -> 0];
gSinv = Simplify[Inverse[gS]];
Chr = Table[
   Simplify[(1/2) Sum[gSinv[[i, l]] (D[gS[[l, j]], coords[[k]]] +
        D[gS[[l, k]], coords[[j]]] - D[gS[[j, k]], coords[[l]]]), {l, 4}]],
   {i, 4}, {j, 4}, {k, 4}];
Riem = Table[
   Simplify[D[Chr[[i, j, l]], coords[[k]]] - D[Chr[[i, j, k]], coords[[l]]] +
     Sum[Chr[[i, k, s]] Chr[[s, j, l]] - Chr[[i, l, s]] Chr[[s, j, k]], {s, 4}]],
   {i, 4}, {j, 4}, {k, 4}, {l, 4}];
RicS = Table[Simplify[Sum[Riem[[s, i, s, j]], {s, 4}]], {i, 4}, {j, 4}];
RS = Simplify[Sum[gSinv[[i, j]] RicS[[i, j]], {i, 4}, {j, 4}]];
einsteinDefect = Simplify[RicS - (RS/4) gS];
check["conformal Schwarzschild (a=0) is NOT Einstein for generic A(eta)",
      Simplify[einsteinDefect] =!= ConstantArray[0, {4, 4}]];
check["it IS Ricci-flat when A is constant (the seed)",
      Simplify[RicS /. A -> Function[x, 1]] === ConstantArray[0, {4, 4}]];

Print["=========================================================="];
Print[If[$Failures === 0, "ALL CHECKS PASSED",
     ToString[$Failures] <> " CHECK(S) FAILED"]];
Print["=========================================================="];
If[$Failures > 0, Exit[1], Exit[0]];
