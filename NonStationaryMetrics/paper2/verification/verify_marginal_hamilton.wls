#!/usr/bin/env wolframscript
(* ::Package:: *)

(* verify_marginal_hamilton.wls -- referee major comment 4 (CQG-116884), part 2.
   The referee asks that the case J = -Jc be rederived "from the original
   Hamilton equations, including one-sided expansions, signs, branch
   orientation, and numerical convergence tests".  verify_cusp_corner.wls
   settles the ORBIT SHAPE dphi/dr.  This script settles the DYNAMICS: whether
   r_e is reached in finite parameter time, and why the two marginal signs
   behave differently even though the shape radical S(r) is even in J.

   Frozen conformal factor (A = 1, Ehat = E), equatorial, tau-branch:

     H_tau = ptil phi0' + R Sqrt[(Delta/r^2) pr^2 + ptil^2/Pbar] - f/E,
     ptil  = J - 2 M a/(r E),        phi0' = (2 M a/r) vbar^2/Pbar,
     R^2   = Delta vbar^2/Pbar,      vbar^2 = 1 - f/E^2,
     Pbar  = P + (2 M a/r)^2/E^2,    P = r^2 + a^2 + 2 M a^2/r,
     f     = 1 - 2M/r,               Delta = r^2 - 2 M r + a^2.

   Key structural fact established below:  at r_e = 2M one has
   ptil(r_e) = J - a/E = J - Jc.  So the PROGRADE marginal J = +Jc kills the
   conformal gravitomagnetic momentum exactly at the stationary limit, while
   the RETROGRADE marginal J = -Jc leaves ptil(r_e) = -2 Jc != 0.  The shape
   radical cannot see this (it depends on J^2); the Hamilton equations do.
*)

$Failures = 0;
check[label_, cond_] := (
  Print[If[TrueQ[cond], "PASS  ", "FAIL  "], label];
  If[! TrueQ[cond], $Failures++]);

Print["=========================================================="];
Print[" marginal |J| = Jc from the Hamilton equations (tau-branch)"];
Print["=========================================================="];

ClearAll[r, M, a, EE, J, pr];
f[r_]     := 1 - 2 M/r;
Dl[r_]    := r^2 - 2 M r + a^2;
P[r_]     := r^2 + a^2 + 2 M a^2/r;
Pbar[r_]  := P[r] + (2 M a/r)^2/EE^2;
vbar2[r_] := 1 - f[r]/EE^2;
phi0p[r_] := (2 M a/r) vbar2[r]/Pbar[r];
R2[r_]    := Dl[r] vbar2[r]/Pbar[r];
ptil[r_]  := J - 2 M a/(r EE);
Jc        = a/EE;

(* ---------------------------------------------------------------- *)
(* 1. ptil(r_e) = J - Jc   -- the sign-breaking quantity             *)
(* ---------------------------------------------------------------- *)
check["ptil(2M) = J - Jc", Simplify[ptil[2 M] - (J - Jc)] === 0];
Print["      ptil(2M)          = ", Simplify[ptil[2 M]]];
Print["      ptil(2M)|_{J=+Jc} = ", Simplify[ptil[2 M] /. J -> Jc], "   (vanishes)"];
Print["      ptil(2M)|_{J=-Jc} = ", Simplify[ptil[2 M] /. J -> -Jc]];

(* Pbar at the stationary limit, as printed in the paper *)
check["Pbar(2M) = 4M^2 + 2a^2 + a^2/E^2",
      Simplify[Pbar[2 M] - (4 M^2 + 2 a^2 + a^2/EE^2)] === 0];

(* ---------------------------------------------------------------- *)
(* 2. solve H_tau = 0 for pr^2 and identify the radial potential     *)
(* ---------------------------------------------------------------- *)
(* H = 0  =>  R Sqrt[(Dl/r^2) pr^2 + ptil^2/Pbar] = f/EE - ptil phi0'   *)
rhs[r_]  := f[r]/EE - ptil[r] phi0p[r];
pr2[r_]  := Simplify[(r^2/Dl[r]) (rhs[r]^2/R2[r] - ptil[r]^2/Pbar[r])];

pr2s = Together[pr2[r]];
Print["----------------------------------------------------------"];
Print["      pr^2(r) computed from H_tau = 0."];

(* the shape radical of the paper *)
w[r_] := (EE^2 - 1) r + 2 M;
Ssex[r_] := r (r - 2 M) w[r] (r Dl[r] - J^2 w[r]);

(* pr^2 must be proportional to Ssex up to a nonvanishing rational factor;
   verify that the ratio has no zero or pole at r = 2M for generic J.       *)
ratio = Simplify[Together[pr2s/Ssex[r]]];
Print["      pr^2 / S  = ", Simplify[ratio]];
ratAt2M = Simplify[ratio /. r -> 2 M];
check["pr^2 is S(r) times a factor regular and nonzero at r_e",
      FreeQ[ratAt2M, DirectedInfinity] && Simplify[ratAt2M] =!= 0];
Print["      (pr^2/S)|_{r=2M} = ", Simplify[ratAt2M]];

(* multiplicity of the root r = 2M of pr^2 *)
mult[expr_, pt_] := Module[{k = 0, e = Together[expr]},
   While[Simplify[Numerator[Together[e]] /. r -> pt] === 0 && k < 6,
     e = D[e, r]; k++]; k];
check["generic J: pr^2 has a SIMPLE zero at r_e", mult[pr2s, 2 M] === 1];
pr2marg = Simplify[pr2s /. J^2 -> Jc^2];
check["|J| = Jc: pr^2 has a DOUBLE zero at r_e", mult[pr2marg, 2 M] === 2];

(* ---------------------------------------------------------------- *)
(* 3. the radial Hamilton equation and the time to reach r_e         *)
(*      dr/ds = dH/dpr = R^2 (Dl/r^2) pr / rhs                       *)
(* ---------------------------------------------------------------- *)
drds[r_, prv_] := R2[r] (Dl[r]/r^2) prv/rhs[r];

(* denominator at the stationary limit: rhs(2M) = -Jc? no -- compute it *)
rhs2M = Simplify[rhs[2 M]];
Print["----------------------------------------------------------"];
Print["      rhs(2M) = f/E - ptil phi0'  at r_e = ", Simplify[rhs2M]];
check["rhs(2M) = -a (J - Jc)/Pbar(2M)",
      Simplify[rhs2M + a (J - Jc)/Pbar[2 M]] === 0];
check["rhs(2M) vanishes exactly at the PROGRADE marginal J = +Jc",
      Simplify[rhs2M /. J -> Jc] === 0];
check["rhs(2M) is nonzero at the RETROGRADE marginal J = -Jc",
      Simplify[rhs2M /. J -> -Jc] =!= 0];
Print["      rhs(2M)|_{J=-Jc} = ", Simplify[rhs2M /. J -> -Jc]];

(* ---------------------------------------------------------------- *)
(* 4. numerics at the paper's representative point M=1, a=0.9, E=1.2 *)
(* ---------------------------------------------------------------- *)
Print["----------------------------------------------------------"];
Print[" numerics at M = 1, a = 0.9, E = 1.2   (Jc = 0.75)"];
Print["----------------------------------------------------------"];
num = {M -> 1, a -> 9/10, EE -> 6/5};
JcN = 3/4;

pr2N[rv_, Jv_] := N[pr2s /. num /. J -> Jv /. r -> rv, 40];
drdsN[rv_, Jv_] := Module[{p2 = pr2N[rv, Jv], pv},
   pv = -Sqrt[p2];                       (* ingoing branch: pr < 0 *)
   N[drds[r, pv] /. num /. J -> Jv /. r -> rv, 40]];

(* (a) retrograde marginal: dr/ds vanishes LINEARLY in (r - 2M)      *)
Print["      J = -Jc :   eps        pr             dr/ds        (dr/ds)/eps"];
Do[Module[{ep = e, p, v},
    p = -Sqrt[pr2N[2 + ep, -JcN]];
    v = drdsN[2 + ep, -JcN];
    Print["                 ", ScientificForm[N[ep, 3]], "   ",
      N[p, 10], "   ", N[v, 10], "   ", N[v/ep, 10]]],
  {e, {1/10^2, 1/10^4, 1/10^6, 1/10^8}}];

lin1 = drdsN[2 + 1/10^6, -JcN]/(1/10^6);
lin2 = drdsN[2 + 1/10^8, -JcN]/(1/10^8);
check["J = -Jc: dr/ds ~ C (r - r_e), C finite nonzero (linear vanishing)",
      Abs[N[lin1 - lin2, 30]] < 10^-4 && Abs[N[lin2, 30]] > 10^-3];
Print["      linear rate C = lim (dr/ds)/(r-r_e) = ", N[lin2, 12]];

(* the elapsed parameter to reach r_e diverges logarithmically *)
tof[eps_, Jv_] := NIntegrate[
   Evaluate[1/Abs[drds[r, -Sqrt[pr2s]] /. num /. J -> Jv]],
   {r, 2 + eps, 2 + 1/10}, WorkingPrecision -> 30, AccuracyGoal -> 20];
t1 = tof[1/10^4, -JcN]; t2 = tof[1/10^8, -JcN]; t3 = tof[1/10^12, -JcN];
Print["      elapsed parameter from r=2.1 down to r_e+eps:"];
Print["        eps=1e-4 : ", N[t1, 10]];
Print["        eps=1e-8 : ", N[t2, 10]];
Print["        eps=1e-12: ", N[t3, 10]];
(* logarithmic divergence: equal decades of eps add equal increments *)
d1 = t2 - t1; d2 = t3 - t2;
Print["      increments per 4 decades: ", N[d1, 10], " , ", N[d2, 10]];
check["J = -Jc: r_e is reached only asymptotically (log-divergent parameter)",
      N[d2, 30] > 1 && Abs[N[(d2 - d1)/d1, 30]] < 10^-2];
check["J = -Jc: increment matches -ln(10^4)/C",
      Abs[N[d2 - Log[10^4]/Abs[lin2], 30]] < 10^-2];

(* closed form of the linear rate C at the retrograde marginal.
   Assembling  C = [R^2 Delta/r^2]_{r_e} (dpr/deps) / rhs(2M)|_{J=-Jc}
   with  [R^2 Delta/r^2]_{r_e} = a^4/(4 M^2 Pbar_e),
         dpr/deps               = -Sqrt[4M^2 + Jc^2]/(E a^2),
         rhs(2M)|_{J=-Jc}       = 2 a^2/(E Pbar_e),
   every parameter cancels except Jc:                                     *)
Cclosed = -Sqrt[4 M^2 + Jc^2]/(8 M^2);
check["closed form of the retrograde linear rate: C = -Sqrt(4M^2+Jc^2)/(8M^2)",
      Abs[N[Cclosed /. num, 30] - N[lin2, 30]] < 10^-6];
Print["      closed form C = -Sqrt(4M^2 + a^2/E^2)/(8 M^2) = ",
      N[Cclosed /. num, 12]];
check["log-divergence increment equals ln(10^4)/|C|",
      Abs[N[d2 - Log[10^4]/Abs[Cclosed /. num], 30]] < 10^-2];

(* (b) sub-marginal |J| < Jc: r_e IS reached in finite parameter time.
   The integrand has an integrable 1/Sqrt[r-r_e] endpoint singularity, so
   substitute r = 2 + u^2 to remove it before integrating.                *)
Jsub = 6/10;
tofReg[epsLo_, Jv_] := Module[{integ},
   integ = Evaluate[2 u/Abs[drds[r, -Sqrt[pr2s]] /. num /. J -> Jv
        /. r -> 2 + u^2]];
   NIntegrate[integ, {u, Sqrt[epsLo], Sqrt[1/10]},
     WorkingPrecision -> 30, AccuracyGoal -> 20, PrecisionGoal -> 20]];
t1b = tofReg[1/10^4, -Jsub]; t2b = tofReg[1/10^8, -Jsub];
t3b = tofReg[1/10^16, -Jsub]; t4b = tofReg[0, -Jsub];
Print["      J = -0.6 (< Jc), same integral, singularity removed by r = 2+u^2:"];
Print["        eps=1e-4  : ", N[t1b, 12]];
Print["        eps=1e-8  : ", N[t2b, 12]];
Print["        eps=1e-16 : ", N[t3b, 12]];
Print["        eps=0     : ", N[t4b, 12]];
(* the truncated tail scales like Sqrt[eps] -- the signature of an integrable
   inverse-square-root endpoint, i.e. finite arrival parameter.            *)
tail1 = Abs[N[t4b - t2b, 30]];        (* eps = 1e-8  -> tail ~ 1e-4 *)
tail2 = Abs[N[t4b - t3b, 30]];        (* eps = 1e-16 -> tail ~ 1e-8 *)
Print["      truncated tails: ", N[tail1, 8], " , ", N[tail2, 8],
      "   ratio = ", N[tail1/tail2, 8], " (expect 1e4)"];
check["|J| < Jc: r_e reached in FINITE parameter time (integral converges)",
      N[t4b, 30] < 10 && tail2 < 10^-6 && Abs[N[tail1/tail2, 30] - 10^4] < 10^2];

(* (c) prograde marginal: rhs -> 0 too, so dr/ds is a 0/0 limit       *)
Print["      J = +Jc :   eps        pr             rhs          dr/ds"];
Do[Module[{ep = e, p, rh, v},
    p = -Sqrt[pr2N[2 + ep, JcN]];
    rh = N[rhs[r] /. num /. J -> JcN /. r -> 2 + ep, 40];
    v = drdsN[2 + ep, JcN];
    Print["                 ", ScientificForm[N[ep, 3]], "   ",
      N[p, 10], "   ", N[rh, 10], "   ", N[v, 10]]],
  {e, {1/10^2, 1/10^4, 1/10^6, 1/10^8}}];
vP = drdsN[2 + 1/10^8, JcN];
check["J = +Jc: dr/ds tends to a FINITE nonzero limit (0/0 resolves)",
      Abs[N[vP, 30]] > 10^-3 &&
      Abs[N[vP - drdsN[2 + 1/10^6, JcN], 30]] < 10^-3];
Print["      dr/ds -> ", N[vP, 12], "  : the prograde marginal reaches and",
      " crosses r_e in finite parameter"];

(* both numerator and denominator vanish linearly; the ratio of their slopes
   is the crossing rate.  Get it symbolically.                             *)
rhsSlope = Simplify[D[rhs[r], r] /. J -> Jc /. r -> 2 M];
Print["      rhs'(2M)|_{J=+Jc} = ", rhsSlope];
prSlope = -Sqrt[4 M^2 + Jc^2]/(EE a^2);
check["dpr/deps at r_e (both marginals) = -Sqrt(4M^2+Jc^2)/(E a^2)",
      Abs[N[prSlope /. num, 30]
          - N[-Sqrt[pr2N[2 + 1/10^8, JcN]]/(1/10^8), 30]] < 10^-6];
prefac = a^4/(4 M^2 Pbar[2 M]);
vPclosed = Simplify[prefac prSlope/rhsSlope];
Print["      closed form of the prograde crossing rate = ", vPclosed];
check["closed form matches the numerical prograde crossing rate",
      Abs[N[vPclosed /. num, 30] - N[vP, 30]] < 10^-6];
check["prograde crossing rate is strictly negative (genuinely ingoing)",
      N[vPclosed /. num, 30] < 0];

(* the two marginal signs therefore differ dynamically *)
check["the two marginal signs are dynamically INEQUIVALENT",
      Abs[N[vP, 30]] > 10^-3 && Abs[N[lin2 (1/10^8), 30]] < 10^-6];

Print["=========================================================="];
Print[If[$Failures === 0, "ALL CHECKS PASSED",
     ToString[$Failures] <> " CHECK(S) FAILED"]];
Print["=========================================================="];
If[$Failures > 0, Exit[1], Exit[0]];
