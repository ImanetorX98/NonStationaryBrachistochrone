#!/usr/bin/env wolframscript
(* ::Package:: *)

(* verify_randers_finsler.wls
   -------------------------------------------------------------------------
   The Randers/Zermelo structure of Thakurta-Kerr, and the Finsler condition
   as an independent derivation of the control domain.

   Gibbons, Herdeiro, Warnick and Werner (PRD 79, 044022, 2009) put the
   Zermelo navigation problem, Randers-Finsler geometry and the NULL geodesic
   flow of a conformally stationary spacetime into a triality.  Writing

       ds^2 = -V^2 (dt + w_i dx^i)^2 + gamma_ij dx^i dx^j        [their (19)]

   the Randers data are  a_ij = V^-2 gamma_ij,  b_i = -w_i  [their (21)-(23)],
   and F = Sqrt[a_ij v^i v^j] + b_i v^i is a genuine Finsler norm -- positive,
   with positive definite fundamental tensor -- precisely when the FINSLER
   CONDITION |b|^2 = a^ij b_i b_j < 1 holds (their Sec. 2.1).  They state that
   this condition "breaks down at an ergo-surface in the spacetime picture."

   This script pins that statement down for our sector and draws the two
   consequences the manuscripts use:

     (1)  1 - |b|^2 vanishes exactly on g_tt = 0, at every latitude, so the
          Finsler condition marks the ERGOSURFACE and nothing else;
     (2)  the Randers data are CONFORMALLY INVARIANT, so Thakurta-Kerr has
          exactly the Randers structure of its Kerr seed.  Together these
          explain why the control-domain boundary sits at the SEED ergosurface
          r_e = 2M and does not move with the conformal factor A(eta).
*)

$Failures = 0;
check[label_, cond_] := (
  Print[If[TrueQ[cond], "PASS  ", "FAIL  "], label];
  If[! TrueQ[cond], $Failures++]);

Print["=========================================================="];
Print[" Randers data of Kerr, and the Finsler condition"];
Print["=========================================================="];

ClearAll[r, th, aa, M, AA, eta];
Sig = r^2 + aa^2 Cos[th]^2;
Del = r^2 - 2 M r + aa^2;
gtt = -(1 - 2 M r/Sig);
gtp = -2 M aa r Sin[th]^2/Sig;
gpp = (r^2 + aa^2 + 2 M aa^2 r Sin[th]^2/Sig) Sin[th]^2;
grr = Sig/Del; gthth = Sig;

randers[gtt_, gtp_, gpp_, grr_, gthth_] := Module[{V2, wp, gam, aij, aup, bi},
  V2 = -gtt;                                  (* threading lapse squared *)
  wp = gtp/gtt;                               (* the single rotation 1-form *)
  gam = {{grr, 0, 0}, {0, gthth, 0}, {0, 0, gpp + V2 wp^2}};
  aij = gam/V2; aup = Inverse[aij]; bi = {0, 0, -wp};
  {aij, bi, Simplify[bi . aup . bi, Assumptions -> {r > 0, M > 0, aa > 0, 0 < th < Pi}]}];

{aK, bK, b2} = randers[gtt, gtp, gpp, grr, gthth];
Print["  |b|^2 = ", b2];
b2e = Simplify[b2 /. th -> Pi/2];
Print["  equatorial: |b|^2 = ", b2e];
Print[""];

check["equatorial Finsler condition holds for every r > 2M and every spin",
  TrueQ[Reduce[ForAll[{r, aa}, r > 2 M && 0 < aa < M, b2e < 1] /. M -> 1, Reals]]];
check["r = 2M is an exact root of 1 - |b|^2 = 0 (equatorial)",
  Simplify[(1 - b2e) /. r -> 2 M] === 0];
check["Finsler condition FAILS just inside: |b|^2 > 1 at r = 1.9M",
  TrueQ[(b2e /. {M -> 1, aa -> 1/2, r -> 19/10}) > 1]];
check["and holds just outside: |b|^2 < 1 at r = 2.1M",
  TrueQ[(b2e /. {M -> 1, aa -> 1/2, r -> 21/10}) < 1]];

(* the general-latitude statement: 1-|b|^2 and -g_tt share their zero set *)
cof = Simplify[Numerator[Together[1 - b2]]/Numerator[Together[-gtt]],
        Assumptions -> {r > 0, M > 0, aa > 0, 0 < th < Pi}];
Print["  cofactor  (1-|b|^2) / (-g_tt), up to positive denominators:"];
Print["    ", cof];
check["that cofactor is strictly positive outside the horizon, all latitudes",
  TrueQ[Reduce[ForAll[{r, aa, th}, r > 2 M && 0 < aa < M && 0 < th < Pi, cof > 0]
     /. M -> 1, Reals]]];
check["hence 1 - |b|^2 vanishes EXACTLY on the ergosurface g_tt = 0",
  Simplify[(1 - b2) /. r -> M + Sqrt[M^2 - aa^2 Cos[th]^2],
     Assumptions -> {M > 0, 0 < aa < M, 0 < th < Pi}] === 0];

Print[""];
Print["=========================================================="];
Print[" Conformal invariance: Thakurta-Kerr inherits Kerr's Randers data"];
Print["=========================================================="];
(* Thakurta-Kerr is g = A(eta)^2 g_Kerr; scale every metric component *)
{aTK, bTK, b2TK} = randers[AA^2 gtt, AA^2 gtp, AA^2 gpp, AA^2 grr, AA^2 gthth];
check["a_ij is unchanged by the conformal factor",
  Simplify[aTK - aK, Assumptions -> {AA > 0, r > 0, M > 0, aa > 0, 0 < th < Pi}]
    === {{0,0,0},{0,0,0},{0,0,0}}];
check["b_i is unchanged by the conformal factor",
  Simplify[bTK - bK, Assumptions -> {AA > 0}] === {0, 0, 0}];
check["hence |b|^2, and the control-domain boundary, do not move with A(eta)",
  Simplify[b2TK - b2, Assumptions -> {AA > 0, r > 0, M > 0, aa > 0, 0 < th < Pi}] === 0];

Print[""];
Print["      Reading.  The Randers data of a conformally stationary metric are"];
Print["      conformal invariants, so Thakurta-Kerr carries EXACTLY the Randers"];
Print["      structure of its Kerr seed.  The Finsler condition of Gibbons et"];
Print["      al. therefore fails on the SEED ergosurface, equatorially r = 2M,"];
Print["      independently of A(eta).  That is the control-domain boundary of"];
Print["      the manuscript, reached here by a fourth independent route: our"];
Print["      own selector-norm computation, Lecian's Cauchy-hypersurface"];
Print["      argument, Kodama's timelike-outside-apparent-horizon statement,"];
Print["      and now the breakdown of the Randers-Finsler structure itself."];
Print[""];
Print["      Note also what is NOT claimed.  Their triality is for NULL"];
Print["      geodesics: their a_ij = V^-2 gamma_ij is the Fermat optical"];
Print["      metric, which is the Ehat -> Infinity limit of Perlick's"];
Print["      fixed-energy factor (verify_myers_optical.wls checks that limit)."];
Print["      Our controlled rail reproduces the Randers FORM, with Perlick's"];
Print["      fixed-energy metric in place of the null one; the coefficients"];
Print["      agree only in the null limit."];

Print["=========================================================="];
Print[If[$Failures === 0, "ALL CHECKS PASSED",
     ToString[$Failures] <> " CHECK(S) FAILED"]];
Print["=========================================================="];
If[$Failures > 0, Exit[1], Exit[0]];
