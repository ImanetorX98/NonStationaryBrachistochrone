#!/usr/bin/env wolframscript
ok[l_, v_] := Print[If[PossibleZeroQ[Simplify[v]], "  [OK]   ", "  [FAIL] "], l,
   "  ->  ", Simplify[v]];
b = Eh^2 - 1;
Cc = r^3 - 2 m r^2 - J2 b r - 2 m J2;         (* C(r) with J2 = J^2 *)
Print["=== when does C have a DOUBLE root? (that is what would make Phi diverge) ==="];
(* eliminate J2 between C = 0 and C' = 0 *)
J2sol = J2 /. First@Solve[D[Cc, r] == 0, J2];
Print["  C'(r) = 0  =>  J^2 = ", Simplify[J2sol], "   (needs r > 4m/3 for J^2 > 0)"];
elim = Simplify[Numerator[Together[Cc /. J2 -> J2sol]]];
Print["  substituting into C = 0 gives, after removing the factor r:"];
loc = Simplify[Cancel[elim/(-2 r)]];
Print["    ", Collect[loc, r], " = 0"];
ok["that locus equals (Eh^2-1) r^2 - m(Eh^2-4) r - 4 m^2",
   Simplify[loc - ((Eh^2 - 1) r^2 - m (Eh^2 - 4) r - 4 m^2)]];

Print["\n=== the positive root of that quadratic is always below 4m/3 ==="];
P[x_] := (Eh^2 - 1) x^2 - m (Eh^2 - 4) x - 4 m^2;
Print["  P(0)     = ", Simplify[P[0]], "   (negative for m > 0)"];
Print["  P(4m/3)  = ", Simplify[P[4 m/3]], "   = (4/9)(Eh^2-1) m^2 > 0 for Eh > 1"];
ok["P(4m/3) - (4/9)(Eh^2-1) m^2", Simplify[P[4 m/3] - (4/9) (Eh^2 - 1) m^2]];
Print["  leading coefficient Eh^2-1 > 0, P(0)<0<P(4m/3)  =>  the unique positive"];
Print["  root lies strictly in (0, 4m/3), where J^2 = (3r^2-4mr)/(Eh^2-1) < 0."];
Print["  Hence NO double root with J^2 > 0: for every physical J the turning point"];
Print["  is a SIMPLE root of C, the integrand behaves as (r-r_t)^(-1/2), and the"];
Print["  apsidal integral Phi converges.  (This is the same fact that puts the"];
Print["  genus degeneration at r_d < 0 in appendix B.)"];
