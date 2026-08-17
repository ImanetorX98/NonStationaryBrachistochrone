#!/usr/bin/env wolframscript
(* App B residue formulas, verified as a Laurent expansion.
   Generic Q4 and F: nothing Vaidya-specific enters.                      *)
ok[l_, v_] := Print[If[PossibleZeroQ[Simplify[v]], "  [OK]   ", "  [FAIL] "], l,
   "  ->  ", Simplify[v]];

(* r(z) near z_d :  dr/dz = Sqrt[Q4(r)],  r(z_d) = rd  *)
Q[d_] := sd^2 + Q1 d + Q2 d^2/2 + Q3 d^3/6;      (* Q4(rd + delta) *)
n = 5;
del = Sum[c[i] eps^i, {i, 1, n}];
odeq = Series[D[del, eps]^2 - Q[del], {eps, 0, n - 1}];
sol = Solve[Thread[CoefficientList[Normal[odeq], eps][[1 ;; 4]] == 0],
   {c[1], c[2], c[3], c[4]}];
solp = First@Select[sol, (c[1] /. #) === sd &];
delta = Normal[Series[del /. solp, {eps, 0, 4}]];
Print["  delta(eps) = r - rd = ", Collect[delta, eps]];
a1v = Coefficient[delta, eps, 2]/sd;
a2v = Coefficient[delta, eps, 3]/sd;
Print["  a_1 = coeff(eps^2)/sd = ", Simplify[a1v], "   (paper: Q4'(rd)/(4 sd))"];
ok["a_1 - Q4'(rd)/(4 sd)", Simplify[a1v - Q1/(4 sd)]];
Print["  a_2 = coeff(eps^3)/sd = ", Simplify[a2v], "   (paper: Q4''(rd)/12)"];
ok["a_2 - Q4''(rd)/12", Simplify[a2v - Q2/12]];

(* R(z) = F(r)/(r-rd)^3 ; Laurent coefficients b_3, b_2, b_1 *)
Fser = F0 + F1 delta + F2 delta^2/2;
Rser = Series[Fser/delta^3, {eps, 0, 0}];
b3 = SeriesCoefficient[Rser, -3];
b2 = SeriesCoefficient[Rser, -2];
b1 = SeriesCoefficient[Rser, -1];
h0 = F0; h1 = F1 sd; h2 = F2 sd^2/2 + F1 Q1/4;
Print["\n  paper: h_0 = F(rd), h_1 = F'(rd) sd, h_2 = F''(rd) sd^2/2 + F'(rd) Q4'(rd)/4"];
ok["b_3 - h_0/sd^3", Simplify[b3 - h0/sd^3]];
ok["b_2 - (h_1 - 3 a_1 h_0)/sd^3", Simplify[b2 - (h1 - 3 a1v h0)/sd^3]];
ok["b_1 - (h_2 - 3 a_1 h_1 + (6 a_1^2 - 3 a_2) h_0)/sd^3",
   Simplify[b1 - (h2 - 3 a1v h1 + (6 a1v^2 - 3 a2v) h0)/sd^3]];

(* two-torsion poles: r - e is EVEN in eps, so R has no residue there *)
Print["\n  two-torsion pole at a simple root e of Q4:"];
de = Qp eps^2/4 + k3 eps^3 + k4 eps^4 + k5 eps^5;
res = Series[D[de, eps]^2 - (Qp de + Qpp de^2/2 + Qppp de^3/6), {eps, 0, 5}];
k3v = k3 /. First@Solve[SeriesCoefficient[res, 3] == 0, k3];
ok["k_3 = 0 (no eps^3 term)", k3v];
r4 = res /. k3 -> 0;
k4v = k4 /. First@Solve[SeriesCoefficient[r4, 4] == 0, k4];
Print["    k_4 = ", Simplify[k4v], "  (even order: nonzero)"];
r5 = Simplify[SeriesCoefficient[r4 /. k4 -> k4v, 5]];
ok["k_5 = 0 (no eps^5 term)", k5 /. First@Solve[r5 == 0, k5]];
Print["    => r(z) is even about a two-torsion point."];
dele = Qp eps^2/4 + k4v eps^4;
Rte = Series[Nm0/((em - rd0)^3 (Qp dele + Qpp dele^2/2)), {eps, 0, -1}];
ok["b_2^{e_i} - 4 N_m(e)/((e-rd)^3 Q4'(e)^2)",
   Simplify[SeriesCoefficient[Rte, -2] - 4 Nm0/((em - rd0)^3 Qp^2)]];
ok["no residue at the two-torsion pole", Simplify[SeriesCoefficient[Rte, -1]]];
