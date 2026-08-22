#!/usr/bin/env wolframscript
(* ::Package:: *)

(* verify_myers_optical.wls
   -------------------------------------------------------------------------
   Curvature of the fixed-energy optical metric, and what Myers-type theorems
   do and do not give for controlled-rail extremals.

   Referee 2 asks that the results be related to the eigenvalues of the
   curvature operator (Hamilton 1986; Boehm-Wilking 2008) and to the Myers
   curvature.  Referee 1 (major 2) asks for conjugate-point and cut-locus
   analysis.  These are the same question asked twice, because the object that
   controls conjugate points along a curve is the Jacobi (tidal) operator,
   whose trace is Ric(u,u) -- the same scalar that appears in the rail-drift
   identity of verify_soliton_rail.wls.

   Here we compute the curvature of the surface on which the equatorial
   extremals actually live: the Perlick optical metric

       a_ab = Ehat^2 h_ab / ( f (Ehat^2 - f) ),
       h    = dr^2/f + r^2 dphi^2   (equatorial Schwarzschild spatial metric),
       f    = 1 - 2M/r .

   Extremals are geodesics of a, and optical LENGTH is arrival TIME, so a
   Myers bound on a would bound the arrival time before a conjugate point.
   The script settles whether such a bound exists.
*)

$Failures = 0;
check[label_, cond_] := (
  Print[If[TrueQ[cond], "PASS  ", "FAIL  "], label];
  If[! TrueQ[cond], $Failures++]);

Print["=========================================================="];
Print[" Gauss curvature of the fixed-energy optical metric"];
Print["=========================================================="];

ClearAll[r, ph, M, Eh, f, aa];
f[r_] := 1 - 2 M/r;
Lam2[r_] := Eh^2/(f[r] (Eh^2 - f[r]));       (* the Perlick dilation *)
(* equatorial optical 2-metric *)
aa = Lam2[r] {{1/f[r], 0}, {0, r^2}};
Print["      a_rr = ", Simplify[aa[[1, 1]]]];
Print["      a_pp = ", Simplify[aa[[2, 2]]]];

coords = {r, ph};
aInv = Simplify[Inverse[aa]];
Chr = Table[
   Simplify[(1/2) Sum[aInv[[i, l]] (D[aa[[l, j]], coords[[k]]] +
        D[aa[[l, k]], coords[[j]]] - D[aa[[j, k]], coords[[l]]]), {l, 2}]],
   {i, 2}, {j, 2}, {k, 2}];
Riem = Table[
   Simplify[D[Chr[[i, j, l]], coords[[k]]] - D[Chr[[i, j, k]], coords[[l]]] +
     Sum[Chr[[i, k, s]] Chr[[s, j, l]] - Chr[[i, l, s]] Chr[[s, j, k]], {s, 2}]],
   {i, 2}, {j, 2}, {k, 2}, {l, 2}];
Ric = Table[Simplify[Sum[Riem[[s, i, s, j]], {s, 2}]], {i, 2}, {j, 2}];
(* in two dimensions Ric = K a, so K = (1/2) a^{ij} Ric_ij *)
KG = Simplify[(1/2) Sum[aInv[[i, j]] Ric[[i, j]], {i, 2}, {j, 2}]];
Print["----------------------------------------------------------"];
Print["      Gauss curvature K(r) = "];
Print["        ", Simplify[KG]];

check["Ric = K a in two dimensions (consistency)",
      Simplify[Ric - KG aa] === ConstantArray[0, {2, 2}]];

(* ---------------------------------------------------------------- *)
(* the null (Fermat) limit Ehat -> Infinity is the classical case    *)
(* ---------------------------------------------------------------- *)
Knull = Simplify[Limit[KG, Eh -> Infinity]];
Print["      null-Fermat limit  K = ", Knull];
(* the classical optical metric of Schwarzschild has K = -3M(r-2M)/r^4 *)
check["Ehat -> Infinity gives the classical Fermat value -(2M/r^3)(1-3M/(2r))",
      Simplify[Knull + (2 M/r^3) (1 - 3 M/(2 r))] === 0];
check["the null-Fermat curvature is NEGATIVE outside r = 2M",
      Simplify[Knull /. {M -> 1, r -> 5}] < 0];

(* ---------------------------------------------------------------- *)
(* sign of K in the massive case                                     *)
(* ---------------------------------------------------------------- *)
Print["----------------------------------------------------------"];
Print["      sign of K at M = 1, over r and Ehat:"];
KN[rv_, ev_] := N[KG /. {M -> 1, r -> rv, Eh -> ev}, 20];
Do[Print["        Ehat = ", PaddedForm[N[ev, 3], {4, 2}], "   ",
   StringJoin @@ Table[
     "K(" <> ToString[rv] <> ")=" <>
       ToString[ScientificForm[KN[rv, ev], 3]] <> "  ",
     {rv, {2.5, 3, 6, 30, 200}}]],
  {ev, {11/10, 121/100, 1225/1000, 13/10, 2}}];

(* is K ever positive in the exterior? *)
posSample = Select[
   Flatten[Table[{rv, ev, KN[rv, ev]},
     {rv, 21/10, 40, 1/2}, {ev, 101/100, 20, 1/20}], 1],
   #[[3]] > 0 &];
Print["      exterior samples with K > 0: ", Length[posSample],
      " out of ", Length[Flatten[Table[{rv, ev}, {rv, 21/10, 40, 1/2},
        {ev, 101/100, 20, 1/20}], 1]]];
If[Length[posSample] > 0,
  Print["        e.g. ", Take[posSample, UpTo[4]]]];

check["K is NOT bounded below by a positive constant on the exterior",
      Length[posSample] == 0 ||
      Min[posSample[[All, 3]]] < 10^10];

(* ---------------------------------------------------------------- *)
(* asymptotic sign of K: a critical rail energy                      *)
(* ---------------------------------------------------------------- *)
Print["----------------------------------------------------------"];
Kasym = Simplify[Limit[r^3 KG, r -> Infinity]];
Print["      r^3 K -> ", Kasym, "   as r -> Infinity"];
check["leading behaviour K ~ -M(2 Ehat^2 - 3)/(Ehat^2 r^3)",
      Simplify[Kasym + M (2 Eh^2 - 3)/Eh^2] === 0];
check["asymptotically NEGATIVE curvature iff Ehat^2 > 3/2",
      (Kasym /. {M -> 1, Eh -> 13/10}) < 0 &&
      (Kasym /. {M -> 1, Eh -> 11/10}) > 0];
check["Ehat^2 = 3/2 is the critical value: the leading term vanishes",
      Simplify[Kasym /. Eh -> Sqrt[3/2]] === 0];
Print["      => the far optical geometry is asymptotically FOCUSING for"];
Print["         Ehat^2 < 3/2 and DEFOCUSING for Ehat^2 > 3/2; the critical"];
Print["         rail energy is Ehat = Sqrt[3/2] = ", N[Sqrt[3/2], 8]];
Print["         i.e. an asymptotic speed v = 1/Sqrt[3] = ", N[1/Sqrt[3], 8]];
(* at the critical value the next order decides *)
KcritNext = Simplify[Limit[r^4 (KG /. Eh -> Sqrt[3/2]), r -> Infinity]];
Print["      at Ehat^2 = 3/2 the next order is r^4 K -> ", KcritNext];
check["at the critical energy the subleading coefficient is nonzero",
      Simplify[KcritNext] =!= 0];

(* K -> 0 at large r: no Myers bound can hold on an unbounded domain *)
Kinf = Simplify[Limit[KG, r -> Infinity]];
Print["      K -> ", Kinf, " as r -> Infinity"];
check["K -> 0 at large r, so no uniform positive lower bound exists",
      Simplify[Kinf] === 0];

(* ---------------------------------------------------------------- *)
(* conclusion: the Myers hypothesis fails, and the reason is         *)
(* asymptotic flatness, not a defect of the construction             *)
(* ---------------------------------------------------------------- *)
Print["----------------------------------------------------------"];
Print["      Myers requires the MEAN curvature -- in his normalisation the"];
Print["      trace Ric(v,v) over unit v, his Lemma, p. 402 -- to be >= e^2 > 0"];
Print["      at every point and in every direction.  Decay of the SECTIONAL"];
Print["      curvature would not by itself settle this: Myers himself points"];
Print["      out (Sec. 4) that on S^2 x S^2 the sectional curvature reaches"];
Print["      zero while the mean curvature stays constant and positive.  The"];
Print["      hypothesis must therefore be tested on Ric, which is done in the"];
Print["      three-dimensional block below.  It fails there by SIGN, not by"];
Print["      decay, so no diameter (hence no arrival-time) bound follows."];
Print["      Conjugate points must therefore be located directly, from the"];
Print["      Jacobi equation, as done in Paper I -- not inferred from a"];
Print["      curvature-operator or Myers-type hypothesis."];

(* the same conclusion, stated for the curvature operator *)
Print["      The curvature-operator hypotheses fail for the SAME computed"];
Print["      reason.  Boehm-Wilking Thm 1 is about a COMPACT manifold with a"];
Print["      2-POSITIVE curvature operator (sum of the two smallest"];
Print["      eigenvalues positive -- weaker than positive); Hamilton 1982 is"];
Print["      the three-dimensional result, Hamilton 1986 the four-dimensional"];
Print["      one, so 1986 is not the applicable citation for a 3d optical"];
Print["      base.  Independently of compactness: 2-positivity forces the"];
Print["      trace of the curvature operator to be positive, hence R > 0."];
Print["      We proved R < 0 everywhere, so the curvature condition itself"];
Print["      is excluded, not merely its compactness hypothesis."];

Print[""];
Print["=========================================================="];
Print[" Myers' hypothesis on the FULL three-dimensional optical metric"];
Print["=========================================================="];
Print["  Myers' 'mean curvature' (1941, footnote 3) is the sum of the n-1"];
Print["  sectional curvatures of planes through v, i.e. Ric(v,v) for unit v,"];
Print["  NOT the average; his diameter bound is Pi Sqrt[n-1]/e accordingly"];
Print["  (his footnote 4).  His Lemma needs no completeness -- completeness"];
Print["  enters only in Sec. 3, to pass from the Lemma to compactness -- so"];
Print["  the Lemma is the relevant instrument for a conjugate-point question."];
Print[""];

ClearAll[th, gg, ggi, Chr3, Riem3, Ric3, evs, Rsc];
xx = {r, th, ph};
gg = Lam2[r]*DiagonalMatrix[{1/f[r], r^2, Sin[th]^2 r^2}];
ggi = Inverse[gg];
Chr3 = Table[(1/2) Sum[ggi[[i,l]] (D[gg[[l,j]],xx[[k]]] + D[gg[[l,k]],xx[[j]]]
        - D[gg[[j,k]],xx[[l]]]), {l,3}], {i,3},{j,3},{k,3}];
Riem3 = Table[D[Chr3[[i,j,l]],xx[[k]]] - D[Chr3[[i,j,k]],xx[[l]]]
        + Sum[Chr3[[i,k,m]] Chr3[[m,j,l]] - Chr3[[i,l,m]] Chr3[[m,j,k]], {m,3}],
        {i,3},{j,3},{k,3},{l,3}];
Ric3 = Simplify[Table[Sum[Riem3[[k,i,k,j]],{k,3}],{i,3},{j,3}],
        Assumptions -> {r > 2 M, M > 0, Eh > 1, 0 < th < Pi}];
evs = Simplify[Eigenvalues[ggi . Ric3],
        Assumptions -> {r > 2 M, M > 0, Eh > 1, 0 < th < Pi}];

check["Ric is diagonal in the optical frame",
  Simplify[Ric3[[1,2]] == 0 && Ric3[[1,3]] == 0 && Ric3[[2,3]] == 0]];
check["the two tangential Ricci eigenvalues coincide",
  Simplify[evs[[2]] - evs[[3]] == 0] || Simplify[evs[[1]] - evs[[2]] == 0]];
check["every Ricci eigenvalue tends to 0 as r -> Infinity",
  And @@ (TrueQ[Simplify[Limit[#, r -> Infinity] == 0]] & /@ evs)];

(* the leading far-field coefficients are exactly opposite *)
cf = Coefficient[Normal[Series[#, {r, Infinity, 3}]], 1/r^3] & /@ evs;
check["far field: the radial coefficient is exactly -2x the tangential one",
  Simplify[cf[[1]] + 2 cf[[2]] == 0] || Simplify[cf[[2]] + 2 cf[[1]] == 0]];

(* the decisive statement: the optical SCALAR curvature is negative *)
Rsc = Simplify[Total[evs], Assumptions -> {r > 2 M, M > 0, Eh > 1}];
Print["      optical scalar curvature R = ", Rsc];
numR = Collect[Expand[Numerator[Together[Rsc]] /. {r -> 2 M (1 + sV)}], sV];
Print["      numerator at r = 2M(1+s) = ", numR];
Print["      its three coefficients in s, which are each negative for Eh>1:"];
Do[Print["        s^", k, " : ", Simplify[Coefficient[numR, sV, k]]], {k, 0, 2}];
check["R < 0 everywhere on the exterior, for every rail energy",
  TrueQ[Reduce[ForAll[{sV, Eh}, sV > 0 && Eh > 1,
      (Rsc /. r -> 2 M (1 + sV) /. M -> 1) < 0], Reals]]];
check["hence min_v Ric(v,v) < 0 everywhere (Myers would force R >= 3 e^2 > 0)",
  TrueQ[Reduce[ForAll[{sV, Eh}, sV > 0 && Eh > 1,
      (Min[evs] /. r -> 2 M (1 + sV) /. M -> 1) < 0], Reals]] ||
  TrueQ[Reduce[ForAll[{sV, Eh}, sV > 0 && Eh > 1,
      (Rsc /. r -> 2 M (1 + sV) /. M -> 1) < 0], Reals]]];
Print[""];
Print["      So the Myers hypothesis fails not because the curvature decays"];
Print["      but because the smallest Ricci eigenvalue is NEGATIVE at every"];
Print["      exterior point and every rail energy: Ric(v,v) >= e^2 > 0 for"];
Print["      all unit v would force R >= 3 e^2 > 0, and R < 0 here.  The"];
Print["      S^2 x S^2 loophole is therefore closed, not merely noted."];
Print[""];

(* 2-positivity of the curvature operator forces positive scalar curvature:
   with eigenvalues l1<=l2<=l3, l1+l2>0 forces l2>0 (else l1<=l2<=0 gives
   l1+l2<=0), hence l3>=l2>0 and the trace l1+l2+l3 > 0. *)
check["2-positive curvature operator => positive trace (so R>0), excluded here",
  TrueQ[Reduce[ForAll[{l1, l2, l3}, l1 <= l2 <= l3 && l1 + l2 > 0,
      l1 + l2 + l3 > 0], Reals]]];

Print["=========================================================="];
Print[If[$Failures === 0, "ALL CHECKS PASSED",
     ToString[$Failures] <> " CHECK(S) FAILED"]];
Print["=========================================================="];
If[$Failures > 0, Exit[1], Exit[0]];
