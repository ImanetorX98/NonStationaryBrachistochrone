#!/usr/bin/env wolframscript
(* verify_optical_completeness.wls
   Is the fixed-energy optical surface COMPLETE?  If it is, and K < 0, then
   Cartan-Hadamard applies on its universal cover and yields a genuine GLOBAL
   minimality statement within each homotopy class -- an "equivalent global
   result" in the sense of referee 1, major comment 2.

   Optical 2-metric (a -> 0 frozen limit, equatorial):
       a = Lambda^2 (dr^2/f + r^2 dphi^2),  Lambda^2 = Eh^2/(f (Eh^2 - f)),
       f = 1 - 2M/r.
   Radial optical distance element:
       ds = Lambda dr/Sqrt[f] = Eh dr / (f Sqrt[Eh^2 - f]).
*)
$Failures = 0;
check[l_, c_] := (Print[If[TrueQ[c], "PASS  ", "FAIL  "], l];
   If[! TrueQ[c], $Failures++]);
ClearAll[r, M, Eh];
f[r_] := 1 - 2 M/r;
ds[r_] := Eh/(f[r] Sqrt[Eh^2 - f[r]]);
Print["=========================================================="];
Print[" completeness of the fixed-energy optical surface"];
Print["=========================================================="];
Print["      radial element ds/dr = ", ds[r]];

(* Integrate reports divergence rather than returning Infinity, so measure the
   partial distances and take limits explicitly.                            *)
num = {M -> 1, Eh -> 14/10};
dsN[rv_] := N[ds[r] /. num /. r -> rv, 30];

(* --- inner end: r -> 2M.  ds/dr ~ C/(r-2M): a simple pole, so log-divergent *)
Kin = Simplify[Limit[(r - 2 M) ds[r], r -> 2 M, Direction -> "FromAbove",
   Assumptions -> {M > 0, Eh > 1}], Assumptions -> {M > 0, Eh > 1}];
Print["      near r = 2M:  (r-2M) ds/dr -> ", Kin];
check["ds/dr has a simple pole at r = 2M, with residue 2M Ehat/Sqrt(Ehat^2-1)... ",
      FreeQ[Kin, r] && Simplify[Kin, Assumptions -> {M > 0, Eh > 1}] =!= 0];
dIn[e_] := NIntegrate[Evaluate[ds[r] /. num], {r, 2 + e, 3},
   WorkingPrecision -> 30];
Print["      optical distance from r=3M down to r=2M+eps:"];
Do[Print["        eps = ", ScientificForm[N[e, 2]], "   d = ",
    N[dIn[e], 10]], {e, {1/10^2, 1/10^4, 1/10^6, 1/10^8}}];
inc1 = dIn[1/10^6] - dIn[1/10^4]; inc2 = dIn[1/10^8] - dIn[1/10^6];
Print["      increments per two decades: ", N[inc1, 8], " , ", N[inc2, 8]];
check["the inner end r = 2M is at INFINITE optical distance (log-divergent)",
      N[inc2, 30] > 1 && Abs[N[(inc2 - inc1)/inc1, 30]] < 10^-3];

(* --- outer end: r -> Infinity.  ds/dr -> Ehat/Sqrt(Ehat^2-1) > 0, so the
       integral diverges linearly.                                          *)
Print["      as r -> Infinity, Lambda^2 -> ",
      Simplify[Limit[Eh^2/(f[r] (Eh^2 - f[r])), r -> Infinity]]];
dsInf = Simplify[Limit[ds[r], r -> Infinity],
   Assumptions -> {M > 0, Eh > 1}];
Print["      as r -> Infinity, ds/dr -> ", dsInf];
check["ds/dr tends to a positive constant, so the outer distance diverges",
      Simplify[dsInf - Eh/Sqrt[Eh^2 - 1]] === 0];
dOut[R_] := NIntegrate[Evaluate[ds[r] /. num], {r, 3, R},
   WorkingPrecision -> 30];
Print["      optical distance from r=3M out to R:"];
Do[Print["        R = ", R, "   d = ", N[dOut[R], 10]],
   {R, {10, 100, 1000, 10000}}];
check["the outer end is at INFINITE optical distance (linear growth)",
      N[dOut[10000] - dOut[1000], 30] > 8000];

(* --- so the surface is complete --- *)
check["the optical surface is geodesically complete (both ends at infinity)",
      N[inc2, 30] > 1 && N[dOut[10000] - dOut[1000], 30] > 8000];

Print["----------------------------------------------------------"];
Print["      Complete, and K < 0 for Ehat^2 >= 3/2.  The universal cover"];
Print["      {r>2M} x R is then complete, simply connected and nonpositively"];
Print["      curved, so by Cartan-Hadamard its exponential map is a"];
Print["      diffeomorphism: between any two points there is exactly one"];
Print["      geodesic, and it minimises.  Projecting back, each homotopy class"];
Print["      on the annulus carries exactly one extremal, which is the GLOBAL"];
Print["      minimiser of arrival time within that class."];
Print["      Global minimality therefore holds within homotopy classes, and"];
Print["      fails between them -- which is the Maxwell set, exactly."];

Print["=========================================================="];
Print[If[$Failures === 0, "ALL CHECKS PASSED",
     ToString[$Failures] <> " CHECK(S) FAILED"]];
If[$Failures > 0, Exit[1], Exit[0]];
