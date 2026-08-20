#!/usr/bin/env wolframscript
ok[l_, v_] := Print[If[PossibleZeroQ[Simplify[v]], "  [OK]   ", "  [FAIL] "], l,
   "  ->  ", Simplify[v]];
DE = b r + 2 m;                       (* b = Ehat^2 - 1 > 0 *)
Jc = r^2 (r - 2 m)/DE;                (* the turning function J(r) *)

Print["=== A = J'(r) D_E^2 , the right-hand side ==="];
A = Simplify[D[Jc, r] DE^2];
Print["  A = ", Factor[A]];
ok["A - 2r( b r(r-m) + m(3r-4m) )", Simplify[A - 2 r (b r (r - m) + m (3 r - 4 m))]];
Print["  for r > 2m > 0, b > 0:  r-m > m > 0 and 3r-4m > 2m > 0, so A > 0."];

Print["\n=== the pointwise comparison, squared ==="];
L = b r (r - m) + m (3 r - 4 m);
Rr = (b r + 2 m)^2 (r (r - 2 m));
diff = Factor[Expand[L^2 - Rr]];
Print["  L^2 - (b r + 2m)^2 r(r-2m) = ", diff];
Print["\n  at m = 0 : ", Simplify[diff /. m -> 0], "   (saturation: the flat chord)"];

Print["\n=== positivity for r > 2m > 0, b > 0 ==="];
Print["  substitute r = 2m + s with s > 0, and scale m = 1:"];
d2 = Factor[Expand[diff /. {r -> 2 + s, m -> 1}]];
Print["    ", d2];
Print["  expanded in s with b > 0, s > 0:"];
Print["    ", CoefficientList[Expand[d2], s]];
Print["  all coefficients nonnegative for b > 0 ? ",
  Reduce[And @@ (# >= 0 & /@ CoefficientList[Expand[d2], s]) , b, Reals]];
Print[];
Print["  full statement: ", Reduce[diff > 0 && r > 2 m && m > 0 && b > 0,
   {r, m, b}, Reals] =!= False];
ok["Resolve: the inequality holds on the whole physical domain",
   Boole[Resolve[ForAll[{r, m, b}, r > 2 m && m > 0 && b > 0, diff > 0], Reals]] - 1];
