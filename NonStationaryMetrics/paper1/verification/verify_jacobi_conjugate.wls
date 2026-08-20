#!/usr/bin/env wolframscript
ok[l_, v_] := Print[If[PossibleZeroQ[Simplify[v]], "  [OK]   ", "  [FAIL] "], l,
   "  ->  ", Simplify[v]];
DE = (Eh^2 - 1) r + 2 m;
Cc = r^2 (r - 2 m) - J^2 DE;
S  = Expand[r (r - 2 m) DE Cc];
K  = J DE;                       (* G = K/Sqrt[S] = dphi/dr *)
Nm = Expand[S D[K, m] - K D[S, m]/2];      (* mass source   *)
NJ = Expand[S D[K, J] - K D[S, J]/2];      (* Jacobi source *)

Print["=== is the Jacobi source the same object we computed? ==="];
ok["N_J - r^3 (r-2m)^2 D_E^2", Simplify[Expand[NJ - r^3 (r - 2 m)^2 DE^2]]];

Print["\n=== Euler: G = K/Sqrt[S] has weight 2 - 3 = -1 ==="];
G = K/Sqrt[S];
sc[f_] := f /. {r -> al r, m -> al m, J -> al J};
ok["G homogeneous of degree -1", Simplify[sc[G] - al^(-1) G]];
Print["  Euler:  r dG/dr + m dG/dm + J dG/dJ = -G,  i.e."];
Print["          m dG/dm + J dG/dJ = -(r G)'"];
ok["m dG/dm + J dG/dJ + d(r G)/dr", Simplify[m D[G, m] + J D[G, J] + D[r G, r]]];

Print["\n=== hence a polynomial identity between the two sources ==="];
ok["m N_m + J N_J - (r K S' - 2 S (r K)')/2",
   Simplify[Expand[m Nm + J NJ - (r K D[S, r] - 2 S D[r K, r])/2]]];
Print["  the right-hand side is purely a total derivative: no U_k can appear,"];
Print["  so the reductions must satisfy  m c_k^m + J c_k^J = 0."];

Print["\n=== solve both reductions and check the consequence ==="];
red[src_] := Module[{A, C5, idn, sol},
   A = Sum[Subscript[aa, i] r^i, {i, 0, 5}];
   C5 = Sum[Subscript[cc, k] r^k, {k, 0, 4}];
   idn = Expand[2 D[A, r] S - A D[S, r] + 2 S C5 - 2 src];
   sol = First@Solve[Thread[CoefficientList[idn, r] == 0],
      Join[Table[Subscript[aa, i], {i, 0, 5}], Table[Subscript[cc, k], {k, 0, 4}]]];
   {Together[A /. sol], Table[Together[Subscript[cc, k] /. sol], {k, 0, 4}]}];
{Am, cm} = red[Nm];
{Aj, cj} = red[NJ];
Do[ok["c_" <> ToString[k] <> "^J + (m/J) c_" <> ToString[k] <> "^m",
   Simplify[cj[[k + 1]] + (m/J) cm[[k + 1]]]], {k, 0, 4}];
ok["m A^m + J A^J + r K   (algebraic parts)", Simplify[m Am + J Aj + r K]];
