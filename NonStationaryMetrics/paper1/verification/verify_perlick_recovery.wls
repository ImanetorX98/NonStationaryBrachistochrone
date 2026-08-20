#!/usr/bin/env wolframscript
ok[l_, v_] := Print[If[PossibleZeroQ[Simplify[v]], "  [OK]   ", "  [FAIL] "], l,
   "  ->  ", Simplify[v]];
n = 2;
gta = Array[gt, n];
gab = Table[gs[Min[i, j], Max[i, j]], {i, n}, {j, n}];
nu = Array[nn, n];
gtn = gta . nu; gnn = nu . gab . nu;

(* the rail, squared, is a quadratic in t' *)
quad[tp_] := ff (ff - Eh^2) tp^2 + 2 gtn (Eh^2 - ff) tp + (gtn^2 + Eh^2 gnn);

(* eq:perlick-randers *)
beta = gta/ff;
aab = Eh^2 (Outer[Times, gta, gta] + ff gab)/(ff^2 (Eh^2 - ff));
D2 = Simplify[nu . aab . nu];
Frail = beta . nu + sq;                       (* sq stands for Sqrt[D2] *)

Print["=== substitute F_rail into the quadratic, treating sq^2 -> D2 ==="];
sub = Expand[quad[Frail]] /. sq^k_ :> D2^(k/2) /; EvenQ[k];
sub = Expand[sub //. sq^2 -> D2];
Print["  quad(F_rail) = ", Simplify[sub]];
ok["quad(F_rail) identically zero (both the rational and the sq-linear part)",
   Simplify[sub]];
Print["  rational part  : ", Simplify[sub /. sq -> 0]];
Print["  coefficient of sq : ", Simplify[Coefficient[sub, sq]]];

Print["\n=== branch: F_rail is the root with f t' - gt.nu > 0 ==="];
Print["  f F_rail - gt.nu = ", Simplify[ff Frail - gtn], "  = f*sq > 0 for f>0"];
ok["f F_rail - gt.nu - f sq", Simplify[ff Frail - gtn - ff sq]];
Print["  so the unsquared rail  f t' - gt.nu = Ehat Sqrt[...]  picks this root,"];
Print["  the other root giving a negative left-hand side."];

Print["\n=== Perlick identification (threading), symbolic ==="];
om = Array[ww, n];
hab = Table[hh[Min[i, j], Max[i, j]], {i, n}, {j, n}];
gtaT = ff om; gabT = hab - ff Outer[Times, om, om];
ok["beta_a - omega_a", Total[Flatten[Simplify[gtaT/ff - om]]]];
ok["g_ta g_tb + f g_ab - f h_ab",
   Total[Flatten[Simplify[Outer[Times, gtaT, gtaT] + ff gabT - ff hab]]]];
aT = Simplify[Eh^2 (Outer[Times, gtaT, gtaT] + ff gabT)/(ff^2 (Eh^2 - ff))];
ok["a_ab - Ehat^2 h_ab/(f(Ehat^2-f))  [Perlick fixed-energy optical metric]",
   Total[Flatten[Simplify[aT - Eh^2 hab/(ff (Eh^2 - ff))]]]];
ok["Ehat -> infinity limit  a_ab -> h_ab/f  [Fermat-Randers]",
   Total[Flatten[Simplify[Limit[aT, Eh -> Infinity] - hab/ff]]]];
