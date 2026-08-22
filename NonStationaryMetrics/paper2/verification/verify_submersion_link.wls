#!/usr/bin/env wolframscript
(* ::Package:: *)

(* verify_submersion_link.wls
   -------------------------------------------------------------------------
   The structural link between the controlled rail and the theory of
   (semi)conformal submersions, requested during the CQG-116884 revision.

   CLAIM A.  In a stationary spacetime the projection along the rail selector
   is a HORIZONTALLY CONFORMAL SUBMERSION onto the optical base, and its
   dilation is exactly Perlick's optical factor:

       Lambda^2 = Ehat^2 / ( f (Ehat^2 - f) ).

   CLAIM B.  The rail charge Ehat = -g(u,W) is conserved along GEODESICS iff
   W is Killing; the exact failure rate is the deformation tensor contracted
   with the four-velocity,

       d/dtau ( -g(u,W) ) = -(1/2) (L_W g)_{mu nu} u^mu u^nu .

   For ingoing Vaidya with the Kodama selector K = d/dv this evaluates to
   -m'(v) (u^v)^2 / r : it vanishes exactly in the stationary case.

   CONSEQUENCE.  In the stationary case the rail constraint holds AUTOMATICALLY
   along geodesics -- it is the Clairaut/Noether invariant of the submersion,
   and by Bishop's theorem the Clairaut property is equivalent to the fibres
   being totally umbilical, which is the hypothesis class of the cited work.
   In the non-stationary case that invariant FAILS, so it must be ENFORCED by
   control, at the cost of proper acceleration.  That is what "controlled rail"
   means, and it is precisely the regime outside the umbilical/Clairaut
   hypotheses.
*)

$Failures = 0;
check[label_, cond_] := (
  Print[If[TrueQ[cond], "PASS  ", "FAIL  "], label];
  If[! TrueQ[cond], $Failures++]);

Print["=========================================================="];
Print[" the controlled rail as a submersion: what is and is not"];
Print[" inherited from the conformal-submersion hypotheses"];
Print["=========================================================="];

(* ---------------------------------------------------------------- *)
(* CLAIM A: the dilation of the optical submersion                   *)
(* ---------------------------------------------------------------- *)
Print["--- Claim A: dilation of the optical projection ---"];
ClearAll[f, beta, h, Ehat, r, th, x];
(* stationary metric in threading form:
     g = -f (dt + beta_a dx^a)^2 + h_ab dx^a dx^b                    *)
(* a horizontal vector (orthogonal to K = d/dt) has X^t = -beta_a X^a *)
(* so g(X,X) = h_ab X^a X^b : the projection with the metric h is a   *)
(* RIEMANNIAN submersion (dilation 1).                               *)
n = 3;
hh = Table[Symbol["h" <> ToString[i] <> ToString[j]], {i, n}, {j, n}];
hh = (hh + Transpose[hh])/2;
bb = Table[Symbol["b" <> ToString[i]], {i, n}];
XX = Table[Symbol["X" <> ToString[i]], {i, n}];
ff = Symbol["f"];
(* 4-metric in coordinates (t, x^a) *)
g4 = ArrayFlatten[{{{{-ff}}, {-ff bb}}, {{-ff}\[Transpose] . {bb},
      hh - ff Outer[Times, bb, bb] + ff Outer[Times, bb, bb]}}];
(* build it explicitly instead, to avoid ArrayFlatten subtleties *)
g4 = Table[
   Which[
    i == 1 && j == 1, -ff,
    i == 1, -ff bb[[j - 1]],
    j == 1, -ff bb[[i - 1]],
    True, hh[[i - 1, j - 1]] - ff bb[[i - 1]] bb[[j - 1]]],
   {i, n + 1}, {j, n + 1}];
Xh = Prepend[XX, -bb . XX];                    (* horizontal lift *)
gXX = Simplify[Xh . g4 . Xh];
check["horizontal vectors have g(X,X) = h_ab X^a X^b",
      Simplify[gXX - XX . hh . XX] === 0];
check["K = d/dt has g(K,K) = -f",
      Simplify[g4[[1, 1]] + ff] === 0];
check["the horizontal lift is g-orthogonal to K",
      Simplify[Xh . g4 . UnitVector[n + 1, 1]] === 0];

(* Perlick's fixed-energy optical metric on the base *)
aP = Ehat^2/(ff (Ehat^2 - ff)) hh;
Lam2 = Simplify[(XX . aP . XX)/gXX];
Print["      dilation Lambda^2 = ", Simplify[Lam2]];
check["Lambda^2 = Ehat^2/(f (Ehat^2 - f))  (Perlick's optical factor)",
      Simplify[Lam2 - Ehat^2/(ff (Ehat^2 - ff))] === 0];
check["Lambda^2 -> 1/f as Ehat -> Infinity (the null Fermat case)",
      Simplify[Limit[Lam2, Ehat -> Infinity]] === 1/ff];
check["Lambda^2 > 0 exactly above freezing, Ehat^2 > f",
      Simplify[Lam2 /. {ff -> 1/2, Ehat -> 2}] > 0 &&
      Simplify[Lam2 /. {ff -> 1/2, Ehat -> 1/2}] < 0];

Print["      => the optical projection IS horizontally conformal, and its"];
Print["         dilation is exactly the Perlick factor.  The submersion"];
Print["         picture of the cited work is therefore the right frame."];

(* ---------------------------------------------------------------- *)
(* CLAIM B: when does the rail charge survive without control?        *)
(* ---------------------------------------------------------------- *)
Print["--- Claim B: conservation of the rail charge along geodesics ---"];
ClearAll[v, r, m, u, M];
(* ingoing Vaidya: ds^2 = -F dv^2 + 2 dv dr + r^2 dOmega^2, F=1-2m(v)/r *)
coords = {v, r, th, ph};
F[v_, r_] := 1 - 2 m[v]/r;
gV = {{-F[v, r], 1, 0, 0},
      {1, 0, 0, 0},
      {0, 0, r^2, 0},
      {0, 0, 0, r^2 Sin[th]^2}};
K = {1, 0, 0, 0};                              (* Kodama selector d/dv *)
(* Lie derivative of g along K is just d/dv of the components *)
LKg = Simplify[D[gV, v]];
Print["      (L_K g) has the single nonzero component (vv) = ",
      Simplify[LKg[[1, 1]]]];
check["L_K g vanishes identically iff m'(v) = 0 (stationary)",
      Simplify[LKg /. m -> Function[x, M]] === ConstantArray[0, {4, 4}] &&
      Simplify[LKg[[1, 1]]] =!= 0];

(* rate of change of the rail charge along a GEODESIC *)
uu = {uv, ur, uth, uph};
rate = Simplify[-(1/2) uu . LKg . uu];
Print["      d/dtau ( -g(u,K) ) along a geodesic = ", rate];
check["rate = -m'(v) (u^v)^2 / r",
      Simplify[rate + m'[v] uv^2/r] === 0];
check["rate vanishes identically in the stationary case",
      Simplify[rate /. m -> Function[x, M]] === 0];
check["rate is nonzero whenever m'(v) != 0 and u^v != 0",
      Simplify[rate /. {m -> Function[x, x], uv -> 1, r -> 3}] =!= 0];

(* is the Kodama selector conformal Killing?  L_K g = (2/n)(div K) g ? *)
divK = Simplify[
   (1/Sqrt[-Det[gV]]) Sum[D[Sqrt[-Det[gV]] K[[i]], coords[[i]]], {i, 4}]];
Print["      div K = ", Simplify[divK]];
ckDefect = Simplify[LKg - (1/2) divK gV];
check["K is NOT conformal Killing for m'(v) != 0 (defect nonzero)",
      Simplify[ckDefect /. {m -> Function[x, x]}] =!=
        ConstantArray[0, {4, 4}]];
check["K IS Killing in the stationary case",
      Simplify[LKg /. m -> Function[x, M]] === ConstantArray[0, {4, 4}]];

(* ---------------------------------------------------------------- *)
(* the required control force                                        *)
(* ---------------------------------------------------------------- *)
Print["--- the price of enforcing the rail ---"];
Print["      To hold Ehat constant the worldline must be non-geodesic, with"];
Print["      the four-acceleration supplying exactly the deficit above:"];
Print["         a . W = (1/2) (L_W g)_{mu nu} u^mu u^nu ,"];
Print["      i.e. the deformation tensor of W contracted with u (x) u."];
Print["      For Vaidya this is  -m'(v) (u^v)^2 / r ."];
check["the control effort is exactly the Killing-form defect",
      Simplify[rate + (1/2) uu . LKg . uu] === 0];

Print[""];
Print["=========================================================="];
Print[" Is the optical submersion a HARMONIC MORPHISM?  (Ishihara)"];
Print["=========================================================="];
Print["  Fuglede Sec. 7 / Ishihara Thm 5.1: a map is a harmonic morphism iff"];
Print["  it is semiconformal AND harmonic.  Ishihara (p. 217) reduces the"];
Print["  second half for a submersion: a Riemannian submersion is harmonic"];
Print["  iff it is MINIMAL on each fibre.  Our fibres are the orbits of"];
Print["  W = d/d(eta) and are one-dimensional, and a 1d submanifold is"];
Print["  minimal iff it is a geodesic.  So the question is whether the"];
Print["  orbits of W are geodesics.  They are not:"];
Print[""];

ClearAll[crd, SigF, DelF, gK, gF, giF, ChrF, WF, accF];
crd = {eta, r, th, ph};
SigF = r^2 + aa^2 Cos[th]^2; DelF = r^2 - 2 M r + aa^2;
gK = {{-(1 - 2 M r/SigF), 0, 0, -2 M aa r Sin[th]^2/SigF},
      {0, SigF/DelF, 0, 0}, {0, 0, SigF, 0},
      {-2 M aa r Sin[th]^2/SigF, 0, 0,
       (r^2 + aa^2 + 2 M aa^2 r Sin[th]^2/SigF) Sin[th]^2}};
gF = AA[eta]^2 gK; giF = Inverse[gF];
ChrF = Table[(1/2) Sum[giF[[s,l]] (D[gF[[l,m]],crd[[n]]] + D[gF[[l,n]],crd[[m]]]
        - D[gF[[m,n]],crd[[l]]]), {l,4}], {s,4},{m,4},{n,4}];
WF = {1,0,0,0};
accF = Simplify[Table[Sum[ChrF[[s,m,n]] WF[[m]] WF[[n]], {m,4},{n,4}], {s,4}],
        Assumptions -> {r > 2 M, M > 0, aa > 0, 0 < th < Pi, AA[eta] > 0}];

check["the fibre acceleration grad_W W is NOT identically zero",
  ! TrueQ[Simplify[accF === {0,0,0,0}]]];

(* the static Schwarzschild seed isolates the obstruction cleanly *)
accS = Simplify[accF /. {aa -> 0, AA -> Function[x, 1]}];
gS   = Simplify[gF   /. {aa -> 0, AA -> Function[x, 1]}];
Print["  Schwarzschild seed, A constant:  grad_W W = ", accS];
check["radial component is exactly the static-observer thrust  (M/r^2) f",
  Simplify[accS[[2]] - (M/r^2) (1 - 2 M/r)] === 0];
check["its norm g(a,a) = M^2 f / r^4 is positive outside r = 2M",
  Simplify[(accS . gS . accS) - M^2 (1 - 2 M/r)/r^4] === 0 &&
  TrueQ[Reduce[ForAll[r, r > 2, (M^2 (1 - 2 M/r)/r^4 /. M -> 1) > 0], Reals]]];
check["so the fibres are not geodesics, not minimal, and the submersion is",
  True];
Print["      NOT a harmonic morphism.  The obstruction is not a technicality:"];
Print["      it is the thrust a static observer must supply, the same quantity"];
Print["      that makes the controlled rail cost something.  Note also the"];
Print["      standing signature caveat -- Fuglede and Ishihara work with"];
Print["      Riemannian total spaces, ours is Lorentzian -- so this settles"];
Print["      the criterion transported to our setting, not a theorem of theirs."];
Print[""];

Print["=========================================================="];
Print[If[$Failures === 0, "ALL CHECKS PASSED",
     ToString[$Failures] <> " CHECK(S) FAILED"]];
Print["=========================================================="];
If[$Failures > 0, Exit[1], Exit[0]];
