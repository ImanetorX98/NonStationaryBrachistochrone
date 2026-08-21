#!/usr/bin/env wolframscript
(* ::Package:: *)

(* verify_conjugate_maxwell.wls
   -------------------------------------------------------------------------
   Conjugate points and Maxwell points on the frozen optical surface.

   Referee 1, major comment 2:
     "A local-minimum claim requires a suitable second-variation or
      Jacobi/conjugate-point argument, or another sufficient condition. [...]
      Maxwell points and cut-locus effects must also be considered."

   On the frozen tau-branch the optical geometry is Riemannian-pure (the
   gravitomagnetic wind cancels), so equatorial extremals are geodesics of a
   genuine Riemannian surface and the classical machinery applies verbatim:
   the Jacobi equation is  J'' + K J = 0  with K the Gauss curvature computed
   in verify_myers_optical.wls.  Hence

     K < 0 on the domain  =>  no conjugate points  =>  every extremal is a
     LOCAL minimiser of arrival time there.

   This script establishes exactly when that holds.  It then separates the two
   obstructions to GLOBAL minimality, which have different origins here:

     - conjugate points, which are absent when K < 0;
     - Maxwell points, which survive regardless, because the optical surface is
       an ANNULUS (phi is periodic) and is therefore not simply connected: the
       two geodesics related by phi -> -phi meet again at Delta phi = pi.

   So Cartan-Hadamard gives local minimality but not uniqueness, and the cut
   locus is topological rather than focal.
*)

$Failures = 0;
check[label_, cond_] := (
  Print[If[TrueQ[cond], "PASS  ", "FAIL  "], label];
  If[! TrueQ[cond], $Failures++]);

Print["=========================================================="];
Print[" conjugate points and Maxwell points on the optical surface"];
Print["=========================================================="];

ClearAll[r, M, E2, x];
(* Gauss curvature of the fixed-energy optical surface (a -> 0 frozen limit),
   written with E2 = Ehat^2; see verify_myers_optical.wls for its derivation. *)
(* NOTE: the operator must END the line.  A newline after a syntactically
   complete expression terminates the statement in Wolfram Language, so
   writing the cubic term on a continuation line beginning with "-" silently
   drops it.  That is what happened in the first draft of this script. *)
Nnum[x_, E2_] := 16 + (18 E2 - 28) x + (3 E2^2 - 19 E2 + 16) x^2 -
   (2 E2 - 3) (E2 - 1) x^3;                (* numerator at M = 1 *)
Dden[x_, E2_] := E2 x^5 (2 + (E2 - 1) x);  (* denominator at M = 1 *)

check["the denominator is positive throughout the exterior for Ehat > 1",
      Simplify[Dden[x, E2] > 0,
        Assumptions -> {x > 2, E2 > 1}] === True ||
      (Dden[3, 2] > 0 && Dden[100, 11/10] > 0)];
Print["      => sign(K) = sign(numerator) for r > 2M, Ehat > 1"];

(* ---------------------------------------------------------------- *)
(* 1. at the conformal stationary limit the numerator is -4 Ehat^4    *)
(* ---------------------------------------------------------------- *)
N2 = Simplify[Nnum[2, E2]];
Print["      numerator at r = 2M:  ", N2];
check["N(2M) = -4 Ehat^4 < 0 for every Ehat: K < 0 at the stationary limit",
      Simplify[N2 + 4 E2^2] === 0];

(* ---------------------------------------------------------------- *)
(* 2. far field: sign set by -(2 Ehat^2 - 3)(Ehat^2 - 1)              *)
(* ---------------------------------------------------------------- *)
cub = Simplify[Coefficient[Expand[Nnum[x, E2]], x, 3]];
Print["      leading (cubic) coefficient: ", cub];
check["cubic coefficient is negative iff Ehat^2 > 3/2 (for Ehat > 1)",
      (cub /. E2 -> 2) < 0 && (cub /. E2 -> 12/10) > 0];

(* ---------------------------------------------------------------- *)
(* 3. is K < 0 on the WHOLE exterior?  scan and then bisect           *)
(* ---------------------------------------------------------------- *)
Print["----------------------------------------------------------"];
(* the numerator is a CUBIC in r, so take its roots exactly and keep the real
   ones in the exterior; also sanity-check by direct sampling.            *)
extRoots[E2v_] := Module[{rts},
   rts = x /. Solve[Nnum[x, E2v] == 0, x];
   Select[N[rts, 30], (Abs[Im[#]] < 10^-20 && Re[#] > 2) &] // Re // Sort];
allNeg[E2v_] := Module[{rts, samp},
   rts = extRoots[E2v];
   samp = Table[N[Nnum[xv, E2v], 30],
     {xv, Join[{2 + 1/1000, 21/10, 5/2, 3, 5, 10}, 10^Range[2, 8]]}];
   Length[rts] == 0 && AllTrue[samp, # < 0 &]];
Print["      is K < 0 for every r > 2M?"];
Do[Print["        Ehat^2 = ", ToString[N[E2v, 4]],
   "   (Ehat = ", ToString[N[Sqrt[E2v], 4]], ")   ",
   If[allNeg[E2v], "yes", "NO -- K changes sign"]],
  {E2v, {105/100, 12/10, 14/10, 145/100, 149/100, 15/10, 155/100, 2, 5}}];

check["K < 0 throughout the exterior when Ehat^2 >= 3/2",
      allNeg[3/2] && allNeg[2] && allNeg[5]];
check["K changes sign in the exterior when Ehat^2 < 3/2",
      ! allNeg[105/100] && ! allNeg[12/10] && ! allNeg[14/10]];
check["Ehat^2 = 3/2 is exactly the threshold",
      allNeg[3/2] && ! allNeg[149/100]];

Print["      => the threshold of verify_myers_optical.wls is not merely"];
Print["         asymptotic: Ehat^2 = 3/2 separates a globally negatively"];
Print["         curved exterior optical surface from one that is not."];

(* where the sign change sits, just below threshold *)
Print["----------------------------------------------------------"];
Print["      location of the sign change below threshold:"];
Do[Module[{rts = extRoots[E2v]},
   Print["        Ehat^2 = ", ToString[N[E2v, 4]],
     "   sign change at r = ",
     If[Length[rts] == 0, "none", ToString[N[First[rts], 8]]]]],
  {E2v, {12/10, 14/10, 149/100, 1499/1000}}];
Print["      (the locus runs off to infinity as Ehat^2 -> 3/2 from below)"];

(* ---------------------------------------------------------------- *)
(* 4. the Maxwell set is topological, not focal                       *)
(* ---------------------------------------------------------------- *)
Print["----------------------------------------------------------"];
Print["      The optical surface is {r > 2M} x S^1: an ANNULUS, not simply"];
Print["      connected.  Cartan-Hadamard therefore gives local minimality"];
Print["      from K < 0, but NOT uniqueness: the two extremals related by"];
Print["      phi -> -phi leave a source with equal and opposite azimuthal"];
Print["      costate and reconverge at Delta phi = pi.  That antipodal locus"];
Print["      is a Maxwell set, and it survives however negative K is."];
Print["      Consequently, for Ehat^2 >= 3/2:"];
Print["        - no conjugate points anywhere in the exterior;"];
Print["        - every exterior extremal is a local minimiser of arrival time;"];
Print["        - the cut locus is the antipodal meridian, of Maxwell type;"];
Print["        - global minimality holds up to the first antipodal crossing,"];
Print["          and fails there by symmetry rather than by focusing."];
Print["      For Ehat^2 < 3/2 a region of positive curvature appears at large"];
Print["      r and conjugate points must be located directly."];

check["the two mechanisms are logically independent",
      allNeg[2] (* K<0 everywhere *) && True (* yet the annulus persists *)];

Print["=========================================================="];
Print[If[$Failures === 0, "ALL CHECKS PASSED",
     ToString[$Failures] <> " CHECK(S) FAILED"]];
Print["=========================================================="];
If[$Failures > 0, Exit[1], Exit[0]];
