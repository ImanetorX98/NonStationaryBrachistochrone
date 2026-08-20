#!/usr/bin/env wolframscript
ok[l_, v_] := Print[If[PossibleZeroQ[Simplify[v]], "  [OK]   ", "  [FAIL] "], l,
   "  ->  ", Simplify[v]];
b = Eh^2 - 1; DE = b r + 2 m;
Cc = Expand[r^2 (r - 2 m) - J^2 DE];
S = Expand[r (r - 2 m) DE Cc];

Print["=== the turning condition and the exact division ==="];
J2t = rt^2 (rt - 2 m)/(b rt + 2 m);
ok["C(r_t) = 0 when J^2 = r_t^2(r_t-2m)/D_E(r_t)",
   Simplify[Cc /. {r -> rt, J^2 -> J2t}]];
Qq = r^2 + (rt - 2 m) r + (rt (rt - 2 m) - J2t b);
ok["C(r) - (r - r_t) Q(r)  with Q the displayed quadratic",
   Simplify[Expand[(Cc /. J^2 -> J2t) - (r - rt) Qq]]];
Print["  Q(r) = ", Collect[Qq, r], "   (regular at r = r_t)"];

Print["\n=== regularity of the u-substituted integrand ==="];
Print["  r = r_t + u^2  =>  dr = 2u du,  S = (r-r_t) P = u^2 P,  sqrt(S) = u sqrt(P)"];
Print["  so  G dr = (K/sqrt(S)) dr = 2 K du / sqrt(P):  no singularity at u = 0."];
Psub = Simplify[(r (r - 2 m) DE Qq) /. r -> rt + u^2];
ok["P(0) - r_t(r_t-2m) D_E(r_t) Q(r_t)  (finite, nonzero)",
   Simplify[(Psub /. u -> 0) - rt (rt - 2 m) (b rt + 2 m) (Qq /. r -> rt)]];
Print["  P(0) = ", Simplify[Psub /. u -> 0]];

Print["\n=== the two endpoint limits that force an interior maximum ==="];
Print["  r_t -> 2m+ :  J^2 = r_t^2(r_t-2m)/D_E -> 0, so K = J D_E -> 0"];
ok["J^2 -> 0 as r_t -> 2m", Simplify[Limit[J2t, rt -> 2 m]]];
Print["  the integrand carries an overall factor J, so Phi -> 0: a radial plunge"];
Print["  sweeps no azimuth."];
Print["  r_t -> R- :  the integration range [0, sqrt(R-r_t)] shrinks to a point,"];
Print["  the integrand staying bounded, so Phi -> 0 as well."];
Print[];
Print["  Phi is continuous and strictly positive on (2m,R), and vanishes at both"];
Print["  ends => it attains an interior maximum => Phi'(J) = 0 has a solution."];
