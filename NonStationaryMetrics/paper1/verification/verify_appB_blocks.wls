#!/usr/bin/env wolframscript
(* Appendix B of Paper I, block form (Lemma I.K): 4 checks for a GENERIC quartic Q4 = a4 r^4 + ... + a0.
   (1) r(z) near the pole z_inf: r = -1/(sqrt(a4)(z - z_inf)) + B + O(z - z_inf), B = -a3/(4 a4);
   (2)-(3) the clock integrand of T_tau = U_3 - 2m U_2 at z_inf: double-pole coefficient 1/a4 and residue
       -(2B + r_d - 2m)/sqrt(a4), with m symbolic;
   (4) the symmetrised integration by parts  Int A'B = AB/2 + Int (A'B - AB')/2.
   Run: wolframscript -file verify_appB_blocks.wls *)
ok[l_, v_] := Print[If[PossibleZeroQ[Simplify[v]], "  [OK]   ", "  [FAIL] "], l,
   "  ->  ", Simplify[v]];
(* expansion of r(z) at the pole z_infinity, for a GENERIC quartic *)
Q4[x_] := a4 x^4 + a3 x^3 + a2 x^2 + a1 x + a0;
rr = -1/(Sqrt[a4] ep) + c0 + c1 ep + c2 ep^2;
res = Series[D[rr, ep]^2 - Q4[rr], {ep, 0, 0}];
sol = Solve[{SeriesCoefficient[res, -3] == 0, SeriesCoefficient[res, -2] == 0}, {c0, c1}];
c0v = c0 /. First[sol];
Print["  r(z) near z_inf :  r = -1/(sqrt(a4) eps) + c0 + O(eps),  c0 = ", Simplify[c0v]];
ok["c0 + a3/(4 a4)", Simplify[c0v + a3/(4 a4)]];
Print["  => B = -a3/(4 a4) = c0, the constant defined in appendix B (an earlier text left it undefined)."];

(* clock residues at z_infinity from T_tau = U_3 - 2m U_2 *)
Print["\n  integrand of T_tau:  (r^3 - 2 m r^2)/(r - rd)"];
integ = Series[(rr^3 - 2 mm rr^2)/(rr - rd) /. {c0 -> c0v, c1 -> (c1 /. First[sol])},
   {ep, 0, 0}];
e2 = SeriesCoefficient[integ, -2];
e1 = SeriesCoefficient[integ, -1];
Print["   double-pole coefficient  = ", Simplify[e2], "   (paper: e_2^{z_inf} = 1/a4)"];
ok["e_2^{z_inf} - 1/a4", Simplify[e2 - 1/a4]];
Print["   residue                  = ", Simplify[e1]];
Print["   paper: -(2B + rd - 2m)/Sqrt[a4]"];
ok["residue + (2 B + rd - 2 m)/Sqrt[a4]  with B = -a3/(4 a4)",
   Simplify[e1 + (2 (-a3/(4 a4)) + rd - 2 mm)/Sqrt[a4]]];
Print["   (an earlier text wrote -2 for -2m, having set m = 1; the manuscript now carries m.)"];

(* the symmetrised by-parts identity behind the block form *)
Print["\n  symmetrised integration by parts:  Int A' B = (1/2) A B + (1/2) Int (A'B - A B')"];
ok["identity (differentiating both sides)",
   Simplify[D[aa[z] bb[z]/2, z] + (aa'[z] bb[z] - aa[z] bb'[z])/2 - aa'[z] bb[z]]];
