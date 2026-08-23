#!/usr/bin/env wolframscript
(* ::Package:: *)

(* verify_vaidya_tube.wls
   -------------------------------------------------------------------------
   The Kodama selector on ingoing Vaidya, and the signature of the apparent
   horizon tube.

   Ashtekar and Krishnan (PRD 68, 104030, Def. 1) require a dynamical horizon
   to be a SPACE-LIKE three-manifold foliated by marginally trapped surfaces;
   they note the timelike case is excluded only when the dominant energy
   condition holds.  For ingoing Vaidya the criterion is elementary and is
   checked here, together with the Kodama identities Paper I relies on.
*)
$Failures = 0;
check[l_, c_] := (Print[If[TrueQ[c], "PASS  ", "FAIL  "], l]; If[! TrueQ[c], $Failures++]);

crd = {v, r, th, ph};
f  = 1 - 2 m[v]/r;
g4 = {{-f, 1, 0, 0}, {1, 0, 0, 0}, {0, 0, r^2, 0}, {0, 0, 0, r^2 Sin[th]^2}};
K4 = {1, 0, 0, 0};
sq = Simplify[Sqrt[-Det[g4]], Assumptions -> {r > 0, 0 < th < Pi}];

check["K = d/dv is divergence-free",
  Simplify[(1/sq) Sum[D[sq K4[[i]], crd[[i]]], {i, 4}]] === 0];
check["g(K,K) = -f, so K is timelike exactly outside r = 2 m(v)",
  Simplify[(K4 . g4 . K4) + f] === 0];
MS = Simplify[(r/2) (1 - Sum[Inverse[g4][[i,j]] D[r, crd[[i]]] D[r, crd[[j]]],
        {i,4},{j,4}])];
check["the Kodama charge is the Misner-Sharp mass m(v)", Simplify[MS - m[v]] === 0];

(* the apparent-horizon tube r = 2 m(v) *)
T = {1, 2 m'[v], 0, 0};
nrm = Simplify[(T . g4 . T) /. r -> 2 m[v]];
Print["      g(T,T) on the tube = ", nrm];
check["g(T,T) = 4 m'(v) on the tube", Simplify[nrm - 4 m'[v]] === 0];
Print["      m' > 0 accretion  -> spacelike -> dynamical horizon (Def. 1)"];
Print["      m' = 0            -> null      -> isolated / Killing horizon"];
Print["      m' < 0 evaporation-> timelike  -> NOT a dynamical horizon"];

Print[If[$Failures === 0, "ALL CHECKS PASSED",
     ToString[$Failures] <> " CHECK(S) FAILED"]];
If[$Failures > 0, Exit[1], Exit[0]];
