#!/usr/bin/env wolframscript
(* ::Package:: *)

(* verify_exterior_separatrix.wls -- referee major comment 5 (CQG-116884).

   The referee proposes that the EXTERIOR retrograde separatrix replace the
   interior classification as the flagship exact result, quoting

       J = -8.053516,   r_d = 3.513905 M

   and warns that the prograde degeneration of the same classification,

       J =  2.93635,    r_d = 1.51229 M

   lies inside r = 2M and must NOT be called a physical control separatrix for
   the selector W = d/d(eta).

   This script (a) reproduces both numbers to high precision from the printed
   t-branch turning function, (b) confirms which side of r = 2M each lies on,
   (c) checks the closed form of Jc^- printed in the manuscript, and
   (d) certifies that the exterior one satisfies the two control-domain
   inequalities, so that it is a bona fide type-(E) object.

   t-branch turning function (frozen A = 1, Ehat = E, equatorial):
       N(r) = (Pbar - J b vbar^2)^2 - J^2 Delta vbar^2
   with b = 2 M a/r, vbar^2 = 1 - f/E^2, Pbar = P + b^2/E^2,
        P = r^2 + a^2 + 2 M a^2/r, f = 1 - 2M/r, Delta = r^2 - 2 M r + a^2.
   A saddle-node (marginal capture) is a DOUBLE root: N = N' = 0.
*)

$Failures = 0;
check[label_, cond_] := (
  Print[If[TrueQ[cond], "PASS  ", "FAIL  "], label];
  If[! TrueQ[cond], $Failures++]);

Print["=========================================================="];
Print[" exterior retrograde separatrix (referee major 5)"];
Print["=========================================================="];

ClearAll[r, J, M, a, EE];
f[r_]     := 1 - 2 M/r;
Dl[r_]    := r^2 - 2 M r + a^2;
P[r_]     := r^2 + a^2 + 2 M a^2/r;
b[r_]     := 2 M a/r;
vbar2[r_] := 1 - f[r]/EE^2;
Pbar[r_]  := P[r] + b[r]^2/EE^2;
N6[r_, J_] := (Pbar[r] - J b[r] vbar2[r])^2 - J^2 Dl[r] vbar2[r];

num = {M -> 1, a -> 9/10, EE -> 6/5};

(* ---------------------------------------------------------------- *)
(* 1. solve N = N' = 0 at the representative parameters              *)
(* ---------------------------------------------------------------- *)
eqs = {N6[r, J] == 0, D[N6[r, J], r] == 0} /. num;

solRetro = FindRoot[eqs, {{r, 7/2}, {J, -8}}, WorkingPrecision -> 40];
rdR = r /. solRetro; JR = J /. solRetro;
Print["   retrograde saddle-node:"];
Print["     r_d = ", N[rdR, 20]];
Print["     J   = ", N[JR, 20]];
check["retrograde r_d reproduces the referee's 3.513905 M",
      Abs[N[rdR - 3513905/1000000, 30]] < 10^-6];
check["retrograde J reproduces the referee's -8.053516",
      Abs[N[JR + 8053516/1000000, 30]] < 10^-6];
check["retrograde degeneration is EXTERIOR: r_d > 2M", N[rdR, 30] > 2];

solPro = FindRoot[eqs, {{r, 3/2}, {J, 3}}, WorkingPrecision -> 40];
rdP = r /. solPro; JP = J /. solPro;
Print["   prograde degeneration:"];
Print["     r_d = ", N[rdP, 20]];
Print["     J   = ", N[JP, 20]];
check["prograde r_d reproduces the referee's 1.51229 M",
      Abs[N[rdP - 151229/100000, 30]] < 10^-5];
check["prograde J reproduces the referee's 2.93635",
      Abs[N[JP - 293635/100000, 30]] < 10^-5];
check["prograde degeneration is INTERIOR: r_d < 2M", N[rdP, 30] < 2];
rplus = N[1 + Sqrt[1 - (9/10)^2], 30];
Print["     seed Kerr null surface r_+ = ", rplus];
check["prograde degeneration lies outside the seed r_+ but inside 2M",
      N[rdP, 30] > rplus && N[rdP, 30] < 2];

(* ---------------------------------------------------------------- *)
(* 2. it really is a DOUBLE root, not two nearby simple roots        *)
(* ---------------------------------------------------------------- *)
mult[rv_, Jv_] := Module[{k = 0, e = N6[r, J] /. num /. J -> Jv},
   While[Abs[N[e /. r -> rv, 30]] < 10^-18 && k < 5, e = D[e, r]; k++]; k];
check["retrograde: r_d is a double root of N (multiplicity 2)",
      mult[rdR, JR] == 2];
check["prograde:   r_d is a double root of N (multiplicity 2)",
      mult[rdP, JP] == 2];

(* ---------------------------------------------------------------- *)
(* 3. closed form of Jc^- printed in the manuscript                  *)
(*      Jc^- = -Pbar/(vbar (Sqrt[Delta] - b vbar))  evaluated at r_* *)
(* ---------------------------------------------------------------- *)
JcMinusClosed = -Pbar[r]/(Sqrt[vbar2[r]] (Sqrt[Dl[r]] - b[r] Sqrt[vbar2[r]]));
val = N[JcMinusClosed /. num /. r -> rdR, 30];
Print["   printed closed form at r_*: ", val];
check["printed closed form of Jc^- agrees with the saddle-node solution",
      Abs[val - N[JR, 30]] < 10^-12];
check["Jc^- is negative (retrograde)", N[val, 30] < 0];

(* ---------------------------------------------------------------- *)
(* 4. the exterior separatrix satisfies the control-domain           *)
(*    inequalities of Lemma (compactness): r > 2M and vbar^2 > 0     *)
(* ---------------------------------------------------------------- *)
vb = N[vbar2[r] /. num /. r -> rdR, 30];
Print["   at r_d:  vbar^2 = ", vb, ",  g(W,W) = ", N[-f[r] /. num /. r -> rdR, 30]];
check["exterior separatrix: W is timelike there (g(W,W) < 0)",
      N[-f[r] /. num /. r -> rdR, 30] < 0];
check["exterior separatrix: Ehat > |W|, i.e. vbar^2 > 0", vb > 0];
check["exterior separatrix is a bona fide type-(E) object",
      N[rdR, 30] > 2 && vb > 0];

(* the whole approach arc stays exterior *)
check["the entire arc r in [r_d, 6M] is inside the control domain",
      AllTrue[Table[N[vbar2[r] /. num /. r -> rv, 30], {rv, rdR, 6, 1/10}],
              # > 0 &] && N[rdR, 30] > 2];

(* ---------------------------------------------------------------- *)
(* 5. spin dependence: does the retrograde one stay exterior?        *)
(* ---------------------------------------------------------------- *)
Print["----------------------------------------------------------"];
Print["   spin scan at E = 1.2 (r_d of the retrograde saddle-node)"];
scan = Table[Module[{s, rr},
    s = FindRoot[{N6[r, J] == 0, D[N6[r, J], r] == 0}
        /. {M -> 1, a -> av, EE -> 6/5}, {{r, 7/2}, {J, -8}},
        WorkingPrecision -> 30];
    rr = r /. s; {av, N[rr, 12], N[J /. s, 12]}], {av, {1/10, 1/2, 9/10, 99/100}}];
Do[Print["     a = ", N[s[[1]], 3], "   r_d = ", s[[2]], "   J = ", s[[3]]],
   {s, scan}];
check["the retrograde saddle-node stays exterior (r_d > 2M) for all scanned spins",
      AllTrue[scan, #[[2]] > 2 &]];

Print["=========================================================="];
Print[If[$Failures === 0, "ALL CHECKS PASSED",
     ToString[$Failures] <> " CHECK(S) FAILED"]];
Print["=========================================================="];
If[$Failures > 0, Exit[1], Exit[0]];
