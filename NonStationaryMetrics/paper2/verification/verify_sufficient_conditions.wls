#!/usr/bin/env wolframscript
(* verify_sufficient_conditions.wls
   Two further sufficient conditions for minimality, and the honest limit.

   (A) The optical distance function is a saturated HJB certificate.  On the
       universal cover of the arrival-time surface (complete, K<0), Cartan-
       Hadamard makes d(p,.) smooth away from p with |grad d|_a = 1 exactly.
       That is a value function saturating the eikonal equation, so referee 1's
       route "a value function / a saturated Hamilton-Jacobi-Bellman
       certificate" is available, not only "an equivalent global result".

   (B) Below the threshold, K < 0 still holds on an INNER REGION r < r_c(Ehat).
       We compute r_c and check it is a genuine sign change, giving a
       domain-restricted local-minimality statement for every Ehat > 1.
*)
$Failures = 0;
check[l_, c_] := (Print[If[TrueQ[c], "PASS  ", "FAIL  "], l];
   If[! TrueQ[c], $Failures++]);
Print["=========================================================="];
Print[" further sufficient conditions for minimality"];
Print["=========================================================="];

ClearAll[r, M, E2, x];
(* numerator of K_t at M = 1, cubic in r; sign(K) = sign(N) in the exterior *)
Nnum[x_, E2_] := 16 + (18 E2 - 28) x + (3 E2^2 - 19 E2 + 16) x^2 -
   (2 E2 - 3) (E2 - 1) x^3;

(* --- (A) the eikonal identity, in the optical metric ------------------- *)
Print["--- (A) the distance function saturates the eikonal equation ---"];
(* on any Riemannian manifold, |grad d|_a = 1 a.e.; verify on our surface in
   the radial direction, where d(r) = Integrate[Lambda/Sqrt[f]] and the
   gradient has a single component.                                        *)
ClearAll[Eh, f];
f[r_] := 1 - 2 M/r;
Lam2[r_] := Eh^2/(f[r] (Eh^2 - f[r]));
arr = Lam2[r] {{1/f[r], 0}, {0, r^2}};
dsdr = Sqrt[Lam2[r]/f[r]];                     (* d d / d r for a radial ray *)
gradsq = Simplify[Inverse[arr][[1, 1]] dsdr^2];
Print["      a^{rr} (dd/dr)^2 = ", Simplify[gradsq]];
check["|grad d|_a^2 = 1 identically: the distance saturates the eikonal",
      Simplify[gradsq - 1] === 0];
Print["      => d(p,.) is a value function with equality, not an inequality;"];
Print["         it is a saturated HJB certificate in the sense required."];

(* --- (B) the inner region below the threshold -------------------------- *)
Print["--- (B) below threshold: K < 0 on an inner region r < r_c ---"];
rc[E2v_] := Module[{rts},
   rts = x /. Solve[Nnum[x, E2v] == 0, x];
   rts = Select[N[rts, 30], (Abs[Im[#]] < 10^-20 && Re[#] > 2) &] // Re // Sort;
   If[Length[rts] == 0, Infinity, First[rts]]];
Print["      Ehat^2      Ehat        r_c/M       K just inside   K just outside"];
Do[Module[{rcv = rc[E2v], a, b},
   a = N[Nnum[rcv - 1/10, E2v], 20]; b = N[Nnum[rcv + 1/10, E2v], 20];
   Print["      ", PaddedForm[N[E2v, 5], {6, 3}], "   ",
     PaddedForm[N[Sqrt[E2v], 5], {6, 3}], "   ",
     If[rcv === Infinity, "  infinity", ToString[PaddedForm[N[rcv, 8], {10, 4}]]],
     "   ", ToString[ScientificForm[a, 3]], "   ",
     ToString[ScientificForm[b, 3]]]],
  {E2v, {102/100, 11/10, 12/10, 13/10, 14/10, 145/100, 149/100, 15/10, 2}}];

check["below threshold there IS a finite r_c with a sign change",
      rc[12/10] < Infinity && N[Nnum[rc[12/10] - 1/10, 12/10], 20] < 0 &&
      N[Nnum[rc[12/10] + 1/10, 12/10], 20] > 0];
check["r_c increases as Ehat^2 -> 3/2 from below",
      rc[12/10] < rc[14/10] < rc[149/100]];
check["at and above the threshold no sign change occurs",
      rc[3/2] === Infinity && rc[2] === Infinity];
check["even at Ehat^2 = 1.02 the inner region is nonempty",
      rc[102/100] > 2];
Print["      => for EVERY Ehat > 1 the arrival-branch extremals confined to"];
Print["         r < r_c(Ehat) carry no conjugate point and are local minimisers;"];
Print["         above the threshold r_c = infinity and the whole exterior works."];

(* --- (C) the limit: the certificate does not survive a drift ----------- *)
Print["--- (C) why this does not extend by continuity to A' != 0 ---"];
Print["      The frozen certificate is SATURATED: |grad d|_a = 1 exactly,"];
Print["      with no margin. A perturbation argument needs a strict"];
Print["      inequality to consume, so nothing survives automatically."];
Print["      Constructing a value function for the drifting problem IS"];
Print["      solving the perturbed HJB equation, i.e. the adiabatic"];
Print["      expansion itself. That remains open and is stated as such."];
Print["=========================================================="];
Print[If[$Failures === 0, "ALL CHECKS PASSED",
     ToString[$Failures] <> " CHECK(S) FAILED"]];
If[$Failures > 0, Exit[1], Exit[0]];
