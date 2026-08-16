#!/usr/bin/env wolframscript
(* ======================================================================
   Symbolic verification of the formal core of Paper I (JMP variant).

   Covers every algebraic step that paper1_JMP.tex states as proved:

     Lemma I.A   rail indicatrix: tricotomy, Vaidya ellipse, support function
     Theorem I.1 (ii)  H_v is the support function of the indicatrix
     Lemma I.4   det M = -32 (Ehat^2 - 1) disc_r S ; kernel at degeneration
     Lemma I.6   Euler homogeneity of S ; N_tot(r_d) = 0
     Theorem I.5 Randers-Euler, self-similarity, Theta H_v, boundary form

   Every residual printed below must be 0.  Run:  wolframscript -file <this>
   Cross-checked independently in SymPy by verify_paper1_core.py.
   ====================================================================== *)

ok[label_, val_] := Print[
   If[PossibleZeroQ[val] || val === 0, "  [OK]   ", "  [FAIL] "], label,
   "  ->  ", val];

Print["=================================================================="];
Print[" LEMMA I.A -- rail indicatrix"];
Print["=================================================================="];

(* (a) tricotomy in the rest frame of W: gamma is fixed by the rail alone *)
Print["\n(a) rest-frame decomposition u = gamma n + wbar,  n = W/|W|"];
Print["    rail  -u.W = |W| gamma = Ehat   =>  gamma = Ehat/|W|"];
Print["    shell -gamma^2 + |wbar|^2 = -1  =>  |wbar|^2 = Ehat^2/|W|^2 - 1"];
Print["    sign of |wbar|^2 for Ehat > |W| > 0 : ",
  Reduce[Ehat^2/Wn^2 - 1 > 0 && Ehat > 0 && Wn > 0, {Ehat, Wn}]];
Print["    => empty / single point / sphere  according as Ehat </=/> |W|."];
Print["    The single-point case Ehat = |W| is the freezing surface."];

(* (b) Vaidya: eliminate u^v directly from the metric *)
Print["\n(b) ingoing Vaidya, W = K = d_v, |W| = Sqrt[f]"];
f = 1 - 2 m/r;
w = Ehat^2 - f;
shell = -f uv^2 + 2 uv (uv Vr) + r^2 (uv Vphi/r)^2 + 1;  (* g(u,u) = -1   *)
rail = f uv - uv Vr - Ehat;                              (* -u_v = Ehat   *)
uvsol = uv /. First@Solve[rail == 0, uv];
ind = Simplify[Numerator[Together[shell /. uv -> uvsol]]];
ok["indicatrix  ==  r^2 [ (Vr+(E^2-f))^2 + E^2 Vphi^2 - E^2 w ]",
   FullSimplify[Expand[ind - r^2 ((Vr + (Ehat^2 - f))^2
       + Ehat^2 Vphi^2 - Ehat^2 w)]]];
ok["paper's eq:vaidya-ind parametrization satisfies it",
   FullSimplify[ind /. {Vr -> (f - Ehat^2) + Ehat Sqrt[w] Cos[th],
                        Vphi -> Sqrt[w] Sin[th]},
     Assumptions -> {r > 2 m, m > 0, Ehat > 1}]];
Print["    centre = ", f - Ehat^2, " (ingoing), semi-axes = ",
      Ehat Sqrt[w], " and ", Sqrt[w]];

(* (c) Theorem I.1(ii): H_v is the support function of that ellipse *)
Print["\n(c) support function of the ellipse  =  branch Hamiltonian eq:Hv"];
Hv = pr (f - Ehat^2) - 1 + Sqrt[w] Sqrt[Ehat^2 pr^2 + J^2/r^2];
sup = pr (f - Ehat^2) + Sqrt[(pr Ehat Sqrt[w])^2 + (J Sqrt[w]/r)^2];
ok["max_theta (pr r' + J phi') - 1 - H_v",
   FullSimplify[sup - 1 - Hv, Assumptions -> {r > 2 m, m > 0, Ehat > 1, J > 0}]];

Print["\n=================================================================="];
Print[" LEMMA I.4 -- the 11x11 reduction system"];
Print["=================================================================="];

DE = (E2 - 1) rr + 2 mm;
S = Expand[rr (rr - 2 mm) DE (rr^2 (rr - 2 mm) - JJ^2 DE)];
A = Sum[Subscript[aa, i] rr^i, {i, 0, 5}];
C5 = Sum[Subscript[cc, k] rr^k, {k, 0, 4}];
P = Expand[2 D[A, rr] S - A D[S, rr] + 2 S C5];
vars = Join[Table[Subscript[aa, i], {i, 0, 5}], Table[Subscript[cc, k], {k, 0, 4}]];
M = Table[Coefficient[Coefficient[P, rr, n], v], {n, 0, 10}, {v, vars}];
Print["\ndeg_rr S = ", Exponent[S, rr], " , deg_rr P = ", Exponent[P, rr],
      " , matrix ", Dimensions[M]];
ok["det M + 32 (Ehat^2-1) disc_r S",
   Factor[Det[M] + 32 (E2 - 1) Discriminant[S, rr]]];
Print["    => nonsingular  <=>  S squarefree AND Ehat != 1"];

(* kernel at the degeneration locus, E = 7/5, m = 1 *)
sub0 = {mm -> 1, E2 -> 49/25};
rdList = rr /. Solve[((E2 - 1) rr^2 - mm (E2 - 4) rr - 4 mm^2 /. sub0) == 0, rr];
J2of = Simplify[(3 rr^2 - 4 mm rr)/(E2 - 1) /. sub0];
rd = SelectFirst[rdList, N[J2of /. rr -> #] > 0 &];
Jd = Sqrt[Simplify[J2of /. rr -> rd]];
Print["\nr_d = ", N[rd, 16], "   J_deg = ", N[Jd, 16]];
ok["J_deg - paper's closed form 5 Sqrt[3011/3072 + 581 Sqrt[249]/9216]",
   RootReduce[Jd - 5 Sqrt[3011/3072 + 581 Sqrt[249]/9216]]];
Sn = Expand[S /. sub0 /. JJ -> Jd];
ok["S(r_d)", RootReduce[Sn /. rr -> rd]];
ok["S'(r_d)", RootReduce[D[Sn, rr] /. rr -> rd]];
ok["det M at J = J_deg", RootReduce[Det[M /. sub0 /. JJ -> Jd]]];
Print["  rank M at J = J_deg = ", MatrixRank[RootReduce[M /. sub0 /. JJ -> Jd]],
      "  (of 11: kernel is exactly one-dimensional)"];
Q4 = Expand[PolynomialQuotient[Sn, (rr - rd)^2, rr]];
Aker = Expand[(rr - rd) Q4];
ok["S - (r-r_d)^2 Q4", RootReduce[Expand[Sn - (rr - rd)^2 Q4]]];
Print["  deg Q4 = ", Exponent[Q4, rr], " , deg A_kernel = ", Exponent[Aker, rr],
      "  (must be 4 and 5)"];
{quot, rem} = PolynomialQuotientRemainder[
   Expand[2 D[Aker, rr] Sn - Aker D[Sn, rr]], Sn, rr];
ok["remainder of 2A'S - A S' modulo S", RootReduce[Expand[rem]]];
ok["C - (A' - Q4)", RootReduce[Expand[quot - (D[Aker, rr] - Q4)]]];
Print["  deg C = ", Exponent[quot, rr], "  (must be <= 4)"];

Print["\n=================================================================="];
Print[" LEMMA I.6 -- degeneration-family derivative"];
Print["=================================================================="];
ok["Euler:  rr S_rr + mm S_mm + JJ S_JJ - 6 S",
   Simplify[rr D[S, rr] + mm D[S, mm] + JJ D[S, JJ] - 6 S]];
ok["explicit scaling  S(a r, a m, a J) - a^6 S",
   Simplify[(S /. {rr -> al rr, mm -> al mm, JJ -> al JJ}) - al^6 S]];
Print["  => at a double root (S = S' = 0):  m S_m + J S_J = 6S - r S_r = 0"];

KK = JJ DE;
NN = S D[KK, mm] - KK D[S, mm]/2;
NJ = S D[KK, JJ] - KK D[S, JJ]/2;
Ntot = Expand[NN + (JJ/mm) NJ];
sub1 = {mm -> 1, E2 -> 49/25, JJ -> Jd};
ok["dS/dm + (J/m) dS/dJ  at r_d",
   RootReduce[(D[S, mm] + (JJ/mm) D[S, JJ]) /. sub1 /. rr -> rd]];
Print["  N(r_d)   = ", N[RootReduce[NN /. sub1 /. rr -> rd], 12],
      "   (nonzero)"];
Print["  N_J(r_d) = ", N[RootReduce[NJ /. sub1 /. rr -> rd], 12],
      "   (nonzero)"];
ok["N_tot(r_d)", RootReduce[Ntot /. sub1 /. rr -> rd]];

Print["\n=================================================================="];
Print[" THEOREM I.5 -- first-order response and the boundary form"];
Print["=================================================================="];
ok["Randers-Euler   (J d_J + pr d_pr) H - (H + 1)",
   FullSimplify[(J D[Hv, J] + pr D[Hv, pr]) - (Hv + 1)]];
ok["self-similarity  r H_r + m H_m + J H_J   (degree zero, off shell too)",
   FullSimplify[m D[Hv, m] + r D[Hv, r] + J D[Hv, J]]];
ok["explicit scaling H(a r, a m, a J) - H",
   FullSimplify[(Hv /. {r -> a r, m -> a m, J -> a J}) - Hv,
     Assumptions -> {a > 0, r > 2 m, m > 0, Ehat > 1, J > 0}]];
ok["Theta H_v - eq:ThetaHv",
   FullSimplify[m D[Hv, m]
     - m (-2 pr/r + Sqrt[Ehat^2 pr^2 + J^2/r^2]/(r Sqrt[w]))]];
ok["d(r pr)/dlambda - (1 + Theta H) - H   [= 0 identically, so = 0 on shell]",
   FullSimplify[(pr D[Hv, pr] - r D[Hv, r]) - (1 + m D[Hv, m]) - Hv]];
Print["  => on the shell H = 0:  Theta H = d(r pr)/dlambda - 1,"];
Print["     hence  S_D(lambda) = [r pr]_0^lambda - lambda.   (eq:SD-vaidya)"];

Print["\nDONE -- every [OK] above is an exact symbolic zero."];
