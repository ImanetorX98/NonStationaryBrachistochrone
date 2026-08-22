#!/usr/bin/env wolframscript
(* ::Package:: *)

(* verify_conformal_killing_horizon.wls
   -------------------------------------------------------------------------
   Is the Thakurta-Kerr conformal stationary limit surface a conformal Killing
   horizon?

   Sultana and Dyer (Gen. Relativ. Gravit. 37 (2005) 1349, Sec. 1) set out the
   terminology this paper uses and supply the criterion.  For a spacetime with a
   conformal Killing field xi, the CONFORMAL STATIONARY LIMIT SURFACE is the
   hypersurface Sigma_0 where xi becomes null, N = xi.xi = 0.  If Sigma_0 is
   itself a null surface it is a CONFORMAL KILLING HORIZON, and they state the
   necessary and sufficient condition for that: the ROTATION of the conformal
   Killing congruence must vanish on Sigma_0, the surface being non-degenerate.

   This script applies that criterion to our selector.  Two facts make it easy:

     - twist vanishing is equivalent to xi ^ dxi = 0, and that test is
       CONFORMALLY INVARIANT, since xi_a -> Om^2 xi_a sends
       (Om^2 xi) ^ d(Om^2 xi) = Om^4 (xi ^ dxi).  So the Thakurta-Kerr question
       reduces to the seed Kerr one, for any A(eta);
     - on the equatorial CSLS the surviving component is exactly -a/(2M).

   Conclusion: for a =/= 0 the Thakurta-Kerr CSLS is NOT a conformal Killing
   horizon, and the obstruction is the spin itself.  For a = 0 the twist vanishes
   identically, which is why Sultana and Dyer's own conformal-Schwarzschild
   construction does produce one.
*)

$Failures = 0;
check[label_, cond_] := (
  Print[If[TrueQ[cond], "PASS  ", "FAIL  "], label];
  If[! TrueQ[cond], $Failures++]);

Print["=========================================================="];
Print[" Sultana-Dyer criterion on the Thakurta-Kerr selector"];
Print["=========================================================="];

crd = {t, r, th, ph};
Sig = r^2 + aa^2 Cos[th]^2; Del = r^2 - 2 M r + aa^2;
g = {{-(1 - 2 M r/Sig), 0, 0, -2 M aa r Sin[th]^2/Sig},
     {0, Sig/Del, 0, 0}, {0, 0, Sig, 0},
     {-2 M aa r Sin[th]^2/Sig, 0, 0,
      (r^2 + aa^2 + 2 M aa^2 r Sin[th]^2/Sig) Sin[th]^2}};
xiUp = {1, 0, 0, 0};
xiDn = Simplify[g . xiUp];
dxi  = Table[D[xiDn[[b]], crd[[a]]] - D[xiDn[[a]], crd[[b]]], {a,4},{b,4}];
wedge = Simplify[Table[
    xiDn[[a]] dxi[[b,c]] + xiDn[[b]] dxi[[c,a]] + xiDn[[c]] dxi[[a,b]],
    {a,4},{b,4},{c,4}], Assumptions -> {r > 0, M > 0, aa > 0, 0 < th < Pi}];

N2   = Simplify[xiUp . g . xiUp];
csls = M + Sqrt[M^2 - aa^2 Cos[th]^2];
check["the CSLS is r = M + Sqrt(M^2 - a^2 cos^2 th), where N = xi.xi vanishes",
  Simplify[N2 /. r -> csls, Assumptions -> {M > 0, 0 < aa < M, 0 < th < Pi}] === 0];
check["equatorially the CSLS is r_e = 2M",
  Simplify[(csls /. th -> Pi/2) - 2 M, Assumptions -> {M > 0}] === 0];

nzAll = Select[Flatten[wedge], ! TrueQ[Simplify[# == 0]] &];
check["xi ^ dxi is not identically zero for a =/= 0 (the congruence twists)",
  Length[nzAll] > 0];

onS  = Simplify[wedge /. r -> csls, Assumptions -> {M > 0, 0 < aa < M, 0 < th < Pi}];
nzS  = Select[Flatten[onS], ! TrueQ[Simplify[# == 0]] &];
check["the twist does NOT vanish ON the CSLS",  Length[nzS] > 0];

eqTwist = Simplify[First[nzS] /. th -> Pi/2, Assumptions -> {M > 0, 0 < aa < M}];
Print["      surviving component on the equatorial CSLS: ", eqTwist];
check["equatorially that component is exactly -a/(2M)",
  Simplify[eqTwist + aa/(2 M)] === 0];

w0 = Simplify[wedge /. aa -> 0];
check["for the Schwarzschild seed (a = 0) the twist vanishes identically",
  Length[Select[Flatten[w0], ! TrueQ[Simplify[# == 0]] &]] == 0];

Print[""];
Print["      Reading.  By Sultana and Dyer's own necessary and sufficient"];
Print["      condition, the Thakurta-Kerr conformal stationary limit surface is"];
Print["      NOT a conformal Killing horizon whenever a =/= 0: the rotation of"];
Print["      the conformal Killing congruence is -a/(2M) on the equatorial"];
Print["      surface, nonzero precisely when the seed spins.  The test is"];
Print["      conformally invariant, so no A(eta) can repair it.  This is an"];
Print["      independent, published route to the statement of Sec. 2 that r_e is"];
Print["      a control-domain boundary and not a horizon; and it explains why"];
Print["      the conformal-SCHWARZSCHILD constructions in the same literature"];
Print["      do obtain a conformal Killing horizon, since there the twist is"];
Print["      identically zero."];

Print["=========================================================="];
Print[If[$Failures === 0, "ALL CHECKS PASSED",
     ToString[$Failures] <> " CHECK(S) FAILED"]];
Print["=========================================================="];
If[$Failures > 0, Exit[1], Exit[0]];
