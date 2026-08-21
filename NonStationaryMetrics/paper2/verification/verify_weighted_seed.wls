#!/usr/bin/env wolframscript
(* ::Package:: *)

(* verify_weighted_seed.wls
   -------------------------------------------------------------------------
   Where the weighted-soliton framework meets the controlled rail, exactly.

   The cited preprint on curvature-operator eigenvalues works on WEIGHTED
   manifolds e^{-2 f^2} g_{mu nu}, and its section 5 studies the 4-velocity as
   a conformal Killing field, arriving at div(v) = 0 (its Eq. 32).  Its
   Theorem T8 states that a strictly parallel weight, nabla_mu F = 0, is
   necessary and sufficient for the Ricci solitons to be EINSTEIN solitons,
   whence constant sectional curvature and space forms (T7, C1).

   Reading it makes the relation to Thakurta-Kerr precise, and it is closer
   than expected:

     with the weight chosen as e^{-2 f^2} = 1/A(eta)^2, i.e. f^2 = ln A,
     the weighted metric IS the Kerr seed, and W = d/d(eta) is a genuine
     KILLING field of it -- div_seed(W) = 0, matching their Eq. (32).

   So their framework describes exactly the FROZEN level of this construction.
   What separates the two is equally precise:

     their central hypothesis is that the weight is strictly parallel,
     nabla_mu F = 0.  Ours is f^2 = ln A(eta), which is NOT parallel whenever
     A runs -- and the whole content of this paper is what happens then.

   Their T8 then says their solitons are Einstein; we verify Thakurta-Kerr is
   not Einstein at the unweighted level, so the two statements are consistent
   and describe complementary regimes.
*)

$Failures = 0;
check[label_, cond_] := (
  Print[If[TrueQ[cond], "PASS  ", "FAIL  "], label];
  If[! TrueQ[cond], $Failures++]);

Print["=========================================================="];
Print[" the weighted seed: where the two frameworks coincide"];
Print["=========================================================="];

ClearAll[eta, r, th, ph, A, M, a];
coords = {eta, r, th, ph};
Sig = r^2 + a^2 Cos[th]^2;
Del = r^2 - 2 M r + a^2;
gbar = {
  {-(1 - 2 M r/Sig), 0, 0, -2 M a r Sin[th]^2/Sig},
  {0, Sig/Del, 0, 0},
  {0, 0, Sig, 0},
  {-2 M a r Sin[th]^2/Sig, 0, 0,
   (r^2 + a^2 + 2 M a^2 r Sin[th]^2/Sig) Sin[th]^2}};
gTK = Simplify[A[eta]^2 gbar];
W = {1, 0, 0, 0};

(* ---------------------------------------------------------------- *)
(* 1. W is conformal Killing for g, but KILLING for the weighted seed *)
(* ---------------------------------------------------------------- *)
psi = A'[eta]/A[eta];
LWg = Simplify[D[gTK, eta]];
check["W is conformal Killing for the full metric: L_W g = 2 psi g",
      Simplify[LWg - 2 psi gTK] === ConstantArray[0, {4, 4}]];

(* the weighted metric with e^{-2 f^2} = 1/A^2 is the seed *)
gWeighted = Simplify[gTK/A[eta]^2];
check["the weighted metric e^{-2f^2} g with f^2 = ln A IS the Kerr seed",
      Simplify[gWeighted - gbar] === ConstantArray[0, {4, 4}]];
check["W is a genuine KILLING field of the weighted seed: L_W gbar = 0",
      Simplify[D[gWeighted, eta]] === ConstantArray[0, {4, 4}]];

(* their Eq. (32): the divergence of the 4-velocity vanishes.  For the seed *)
divSeed = Simplify[
   (1/Sqrt[-Det[gbar]]) Sum[
     D[Sqrt[-Det[gbar]] W[[i]], coords[[i]]], {i, 4}]];
Print["      div_seed(W) = ", divSeed, "   (their Eq. 32)"];
check["div of W with respect to the weighted seed vanishes, as their Eq. (32)",
      Simplify[divSeed] === 0];

(* whereas for the full metric it does not *)
divFull = Simplify[
   (1/Sqrt[-Det[gTK]]) Sum[
     D[Sqrt[-Det[gTK]] W[[i]], coords[[i]]], {i, 4}]];
Print["      div_g(W)    = ", Simplify[divFull]];
check["div of W with respect to the full metric is 4 psi, nonzero when A runs",
      Simplify[divFull - 4 psi] === 0];

(* ---------------------------------------------------------------- *)
(* 2. their hypothesis is a strictly parallel weight; ours is not     *)
(* ---------------------------------------------------------------- *)
Print["----------------------------------------------------------"];
Print["      their Theorem T8: nabla_mu F = 0 (strictly parallel weight) is"];
Print["      necessary and sufficient for the solitons to be Einstein."];
(* our weight is f^2 = log A(eta): its gradient is nonzero when A runs *)
fsq = Log[A[eta]];
gradf = Simplify[Table[D[fsq, coords[[i]]], {i, 4}]];
Print["      our weight f^2 = log A(eta),  nabla_mu f^2 = ", gradf];
check["our weight is NOT strictly parallel unless A is constant",
      Simplify[gradf] =!= ConstantArray[0, 4] &&
      Simplify[gradf /. A -> Function[x, 1]] === ConstantArray[0, 4]];
check["the gradient is exactly the conformal factor psi",
      Simplify[gradf[[1]] - psi] === 0];

Print["      => their hypothesis holds precisely at frozen A, which is this"];
Print["         paper's frozen limit; the running weight is the whole subject"];
Print["         of the non-stationary theory."];

(* ---------------------------------------------------------------- *)
(* 3. consistency: their T8 forces Einstein; TK is not Einstein       *)
(* ---------------------------------------------------------------- *)
Print["----------------------------------------------------------"];
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
check["conformal Schwarzschild is NOT Einstein for generic A",
      Simplify[RicS - (RS/4) gS] =!= ConstantArray[0, {4, 4}]];
check["it IS Ricci-flat at frozen A -- the seed, where their theory applies",
      Simplify[RicS /. A -> Function[x, 1]] === ConstantArray[0, {4, 4}]];
Print["      Their C1 (the Schwarzschild metrized soliton has vanishing"];
Print["      sectional curvature eigenvalue) is therefore a statement about"];
Print["      the seed, and is consistent with what we find: the seed is"];
Print["      Ricci-flat, the conformally rescaled spacetime is not."];

Print["=========================================================="];
Print[If[$Failures === 0, "ALL CHECKS PASSED",
     ToString[$Failures] <> " CHECK(S) FAILED"]];
Print["=========================================================="];
If[$Failures > 0, Exit[1], Exit[0]];
