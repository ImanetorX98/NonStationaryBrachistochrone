#!/usr/bin/env wolframscript
(* ::Package:: *)
(* Evaluation of the t-branch separatrix closed form, Paper II eq. (t-sep-phi),
     phi(r) = Lam0 z(r) + Sum_k lam_k Log[ sigma(z - v_k) / sigma(z + v_k) ],
   in Mathematica's own Weierstrass functions, against direct quadrature, on the five
   arcs quoted in the caption of fig:t-sep.  Everything is computed from the printed
   definitions at 40 digits: J_c^- and r_star from Q2 = Q2' = 0, the residues alpha_k as
   limits (no fit), A, B, g2, g3 from the quartic Q4 = R6/(r - r_star)^2, which has a root
   at r = 0, and z(r) from InverseWeierstrassP on the pair (P, P'), P' fixed by the
   orientation dz = dr/sqrt(Q4).
   M = 1, a = 9/10, E = 6/5. *)

$MinPrecision = 40; wp = 60;
Off[NIntegrate::precw];
M = 1; a = 9/10; EE = 6/5;
Q2[r_, J_] := 2 EE^2 J^2 M r - EE^2 J^2 r^2 - 4 EE^2 J M a r + 2 EE^2 M a^2 r + EE^2 a^2 r^2 +
   EE^2 r^4 + 4 J^2 M^2 - 4 J^2 M r + J^2 r^2 - 8 J M^2 a + 4 J M a r + 4 M^2 a^2;
DE[r_] := (EE^2 - 1) r + 2 M; Dl[r_] := r^2 - 2 M r + a^2;
sol = FindRoot[{Q2[r, J] == 0, D[Q2[r, J], r] == 0}, {{r, 351/100}, {J, -805/100}},
   WorkingPrecision -> wp];
{rs, Jc} = {r, J} /. sol;
Print["J_c^- = ", N[Jc, 20], "   r_star = ", N[rs, 20]];

R6 = Expand[r Q2[r, Jc] DE[r]];
{Q4, rem} = PolynomialQuotientRemainder[R6, (r - rs)^2, r];
Print["remainder of R6/(r-r_star)^2: ", N[Max[Abs[CoefficientList[rem, r]]], 3]];
Kt = r DE[r] (Jc (r - 2 M) + 2 M a)/Dl[r];
G = Together[Kt/(r - rs)];                     (* dphi/dr = G/sqrt(Q4) *)
rp = M + Sqrt[M^2 - a^2]; rm = M - Sqrt[M^2 - a^2];
cs = {rp, rm, rs};
Lpoly = Limit[G, r -> Infinity];
al = Table[Limit[(r - c) G, r -> c], {c, cs}];
Print["no pole at r=0: G(0) finite = ", FreeQ[G /. r -> 0, ComplexInfinity]];
pf = Lpoly + Total[al/(r - cs)];
Print["partial fractions exact: ", N[Max[Abs[Table[(G - pf) /. r -> x, {x, {4, 7, 13}}]]], 3]];

q = CoefficientList[Q4, r];                    (* {q0, q1, q2, q3, q4}, q0 = 0 *)
Print["Q4(0) = ", N[q[[1]], 3]];
{q1, q2, q3, q4} = q[[2 ;; 5]];
AA = q1/4; BB = q2/12;
cub = Expand[(q1 s^3 + q2 s^2 + q3 s + q4) /. s -> (4 P/q1 - q2/(3 q1))];
cubN = Expand[cub (q1/4)^2];                   (* = 4P^3 - g2 P - g3 *)
Print["leading coefficient 4: ", N[Coefficient[cubN, P, 3], 20], "  P^2 term: ", N[Coefficient[cubN, P, 2], 3]];
g2 = -Coefficient[cubN, P, 1]; g3 = -Coefficient[cubN, P, 0];
inv = {g2, g3};

zOf[x_] := InverseWeierstrassP[{AA/x + BB, -AA Sqrt[Q4 /. r -> x]/x^2}, inv];
vs = Table[InverseWeierstrassP[{AA/c + BB, WeierstrassPPrime[
      InverseWeierstrassP[AA/c + BB, inv], inv]}, inv], {c, cs}];
lam = Table[-AA al[[k]]/(cs[[k]]^2 WeierstrassPPrime[vs[[k]], inv]), {k, 3}];
Lam0 = Lpoly - Total[al/cs] + Total[2 WeierstrassZeta[#, inv] & /@ vs lam];
phiC[z_] := Lam0 z + Sum[lam[[k]] Log[WeierstrassSigma[z - vs[[k]], inv]/
      WeierstrassSigma[z + vs[[k]], inv]], {k, 3}];

(* uniformisation check *)
Print["P(z(r)) - (A/r+B) at r=5: ", N[Abs[WeierstrassP[zOf[5], inv] - (AA/5 + BB)], 3]];

arcs = {{6, 60}, {36/10, 9}, {4, 20}, {355/100, 4}, {352/100, 36/10}};
fails = 0;
Do[
  {r1, r2} = arc;
  quadv = NIntegrate[G/Sqrt[Q4], {r, r1, r2}, WorkingPrecision -> 45, PrecisionGoal -> 35,
    MaxRecursion -> 40];
  dz = NIntegrate[1/Sqrt[Q4], {r, r1, r2}, WorkingPrecision -> 45, PrecisionGoal -> 30];
  z1 = zOf[r1]; z2 = zOf[r2];
  (* lattice representative: z2 - z1 must equal the integral of dr/sqrt(Q4) *)
  If[Abs[(z2 - z1) - dz] > 10^-20, Print["  lattice shift needed on ", arc]; fails++];
  (* the logarithms are summed along a path, so take the principal value of each
     increment and check it is continuous by subdividing *)
  zs = Table[zOf[r1 + (r2 - r1) t], {t, 0, 1, 1/400}];
  inc = Sum[Lam0 (zs[[i + 1]] - zs[[i]]) + Sum[lam[[k]] (Log[WeierstrassSigma[zs[[i + 1]] - vs[[k]], inv]/
          WeierstrassSigma[zs[[i + 1]] + vs[[k]], inv]] - Log[WeierstrassSigma[zs[[i]] - vs[[k]], inv]/
          WeierstrassSigma[zs[[i]] + vs[[k]], inv]]), {k, 3}], {i, Length[zs] - 1}];
  err = Abs[inc - quadv]/Abs[quadv];
  Print["  arc [", N[r1, 4], ", ", N[r2, 4], "]: closed form ", N[Re[inc], 16],
    "  quadrature ", N[quadv, 16], "  rel ", N[err, 3], "  Im ", N[Abs[Im[inc]], 2]];
  If[err > 10^-25, fails++],
  {arc, arcs}];

(* falsifiability: a wrong residue must break the agreement *)
lamBad = lam; lamBad[[3]] = 1.001 lam[[3]];
{r1, r2} = {4, 20};
zs = Table[zOf[r1 + (r2 - r1) t], {t, 0, 1, 1/400}];
incBad = Sum[Lam0 (zs[[i + 1]] - zs[[i]]) + Sum[lamBad[[k]] (Log[WeierstrassSigma[zs[[i + 1]] - vs[[k]], inv]/
        WeierstrassSigma[zs[[i + 1]] + vs[[k]], inv]] - Log[WeierstrassSigma[zs[[i]] - vs[[k]], inv]/
        WeierstrassSigma[zs[[i]] + vs[[k]], inv]]), {k, 3}], {i, Length[zs] - 1}];
quadv = NIntegrate[G/Sqrt[Q4], {r, r1, r2}, WorkingPrecision -> 45, PrecisionGoal -> 35];
If[Abs[incBad - quadv]/Abs[quadv] < 10^-6, Print["falsifiability FAILED"]; fails++,
  Print["falsifiability: lambda_3 off by 0.1% gives rel ", N[Abs[incBad - quadv]/Abs[quadv], 3]]];

Print[If[fails == 0, "ALL CHECKS PASSED", ToString[fails] <> " FAILED"]];
Exit[If[fails == 0, 0, 1]];
