#!/usr/bin/env wolframscript
(* ::Package:: *)

(* verify_thrust_bound.wls
   -------------------------------------------------------------------------
   The proper acceleration required to hold the rail, and where it diverges.

   Referee 1, major comment 9:
     "The rail minimises a selected elapsed clock; it does not minimise thrust,
      fuel, impulse, or integrated acceleration.  Paper I notes that the
      required proper acceleration can diverge near freezing or horizon limits.
      These distinctions should remain visible in the abstract, figures, and
      conclusion."

   That was a caveat in prose.  It can be made a sharp theorem.

   Differentiating the rail charge along an arbitrary (non-geodesic) worldline,
        d/dtau (-g(u,W)) = -g(a,W) - (1/2)(L_W g)_{mu nu} u^mu u^nu,
   with a = nabla_u u the four-acceleration.  Holding the rail therefore forces
        g(a,W) = -(1/2)(L_W g)_{mu nu} u^mu u^nu,
   and for Thakurta-Kerr, where L_W g = 2 (A'/A) g and g(u,u) = -1,
        g(a,W) = A'/A = epsilon.                                          [1]

   Split W into its parts along and orthogonal to u:  since g(a,u) = 0,
   [1] involves only the orthogonal part, and Cauchy-Schwarz in the (positive
   definite) rest space of u gives the SHARP bound
        |a| >= epsilon / |W_perp|,      |W_perp|^2 = g(W,W) + (g(W,u))^2 .
   With g(W,W) = -A^2 f and g(W,u) = -Ehat this is |W_perp|^2 = Ehat^2 - A^2 f
   = Ehat^2 vbar^2, hence
        |a| >= epsilon / (Ehat vbar),   vbar^2 = 1 - A^2 f/Ehat^2 .        [2]

   The bound is attained by a parallel to W_perp -- the minimum-thrust control.
   It DIVERGES at the freezing surface vbar -> 0 and stays FINITE at the
   conformal stationary limit f -> 0.  That distinction is sharp and is the
   content of the referee's caveat.
*)

$Failures = 0;
check[label_, cond_] := (
  Print[If[TrueQ[cond], "PASS  ", "FAIL  "], label];
  If[! TrueQ[cond], $Failures++]);

Print["=========================================================="];
Print[" minimum proper acceleration needed to hold the rail"];
Print["=========================================================="];

(* ---------------------------------------------------------------- *)
(* 1. the identity g(a,W) = epsilon, done abstractly                 *)
(* ---------------------------------------------------------------- *)
n = 4;
gg = Table[Symbol["g" <> ToString[i] <> ToString[j]], {i, n}, {j, n}];
gg = (gg + Transpose[gg])/2;
uu = Table[Symbol["u" <> ToString[i]], {i, n}];
aa = Table[Symbol["a" <> ToString[i]], {i, n}];
eps = Symbol["epsilon"];

(* conformal Killing selector: L_W g = 2 eps g  (Thakurta-Kerr) *)
LWg = 2 eps gg;
(* rail held: 0 = -g(a,W) - (1/2)(L_W g)(u,u)  =>  g(a,W) = -(1/2)(L_W g)(u,u) *)
gaW = Simplify[-(1/2) uu . LWg . uu];
gaWonShell = Simplify[gaW /. Solve[uu . gg . uu == -1, Symbol["g11"]][[1]]];
Print["      g(a,W) = ", gaWonShell];
check["for a conformal Killing selector, g(a,W) = epsilon exactly",
      Simplify[gaWonShell - eps] === 0];

(* ---------------------------------------------------------------- *)
(* 2. the orthogonal decomposition and the sharp bound               *)
(* ---------------------------------------------------------------- *)
Print["----------------------------------------------------------"];
ClearAll[Eh, Af, ff, vb, epsv, w0, w1];
(* work in the local rest frame of u: u = (1,0,0,0), metric = diag(-1,1,1,1) *)
eta4 = DiagonalMatrix[{-1, 1, 1, 1}];
uRest = {1, 0, 0, 0};
(* W has time component Ehat (since -g(u,W) = Ehat) and a spatial part *)
Wrest = {Eh, w1, 0, 0};
check["-g(u,W) = Ehat fixes the time component of W in the rest frame",
      Simplify[-uRest . eta4 . Wrest - Eh] === 0];
(* g(W,W) = -A^2 f  fixes the magnitude of the spatial part *)
sol = Solve[Wrest . eta4 . Wrest == -Af^2 ff, w1];
w1val = w1 /. sol[[2]];
Print["      |W_perp| = ", Simplify[w1val]];
check["|W_perp|^2 = Ehat^2 - A^2 f = Ehat^2 vbar^2",
      Simplify[w1val^2 - (Eh^2 - Af^2 ff)] === 0];

(* the acceleration is orthogonal to u, hence purely spatial *)
aRest = {0, a1, a2, a3};
check["a is orthogonal to u (four-acceleration is spatial)",
      Simplify[aRest . eta4 . uRest] === 0];
check["g(a,W) reduces to a1 |W_perp| (the time part drops, since g(a,u)=0)",
      Simplify[(aRest . eta4 . Wrest /. w1 -> w1val) - a1 w1val] === 0];
(* imposing g(a,W) = eps: a1 = eps/|W_perp|, and |a|^2 = a1^2 + a2^2 + a3^2 *)
a1val = epsv/w1val;
normSq = Simplify[a1val^2 + a2^2 + a3^2];
Print["      |a|^2 = ", Simplify[normSq]];
check["|a| >= epsilon/|W_perp|, with equality iff a is parallel to W_perp",
      Simplify[(normSq /. {a2 -> 0, a3 -> 0}) - epsv^2/w1val^2] === 0 &&
      Simplify[normSq - epsv^2/w1val^2] === a2^2 + a3^2];

aMin = Simplify[epsv/w1val];
Print["      minimum thrust |a|_min = ", aMin];
check["|a|_min = epsilon/(Ehat vbar) with vbar^2 = 1 - A^2 f/Ehat^2",
      Simplify[aMin - epsv/Sqrt[Eh^2 - Af^2 ff]] === 0];

(* ---------------------------------------------------------------- *)
(* 3. where it diverges, and where it does not                       *)
(* ---------------------------------------------------------------- *)
Print["----------------------------------------------------------"];
(* freezing surface: A^2 f -> Ehat^2, i.e. vbar -> 0 *)
(* approach freezing along A^2 f = Ehat^2 (1 - delta), delta -> 0+ *)
aNearFreeze = Simplify[aMin /. Af^2 -> Eh^2 (1 - dlt)/ff,
   Assumptions -> {Eh > 0, epsv > 0, 0 < dlt < 1}];
Print["      near freezing, A^2 f = Ehat^2 (1-delta):  |a|_min = ",
      Simplify[aNearFreeze, Assumptions -> {Eh > 0, epsv > 0, 0 < dlt < 1}]];
check["thrust DIVERGES at the freezing surface (delta -> 0+)",
      Limit[aNearFreeze, dlt -> 0, Direction -> "FromAbove",
        Assumptions -> {Eh > 0, epsv > 0}] === DirectedInfinity[1]];
check["the divergence is the inverse square root of the margin",
      Simplify[aNearFreeze Sqrt[dlt] - epsv/Eh,
        Assumptions -> {Eh > 0, epsv > 0, 0 < dlt < 1}] === 0];
(* conformal stationary limit: f -> 0, so vbar -> 1 *)
aStat = Simplify[aMin /. ff -> 0];
Print["      at the conformal stationary limit f -> 0:  |a|_min = ", aStat];
check["thrust stays FINITE at the conformal stationary limit f -> 0",
      Simplify[aStat - epsv/Eh, Assumptions -> {Eh > 0, epsv > 0}] === 0];
check["and vanishes at frozen A (epsilon = 0): no control needed",
      Simplify[aMin /. epsv -> 0] === 0];

(* numbers at the paper's representative point *)
Print["----------------------------------------------------------"];
Print["      at M = 1, a = 0.9, Ehat = 1.4, A = 1, epsilon = 0.04:"];
num = {epsv -> 4/100, Eh -> 14/10, Af -> 1};
Do[Module[{fv = 1 - 2/rv},
   Print["        r = ", ToString[N[rv, 4]],
     "   vbar = ", ToString[N[Sqrt[1 - fv/(14/10)^2], 6]],
     "   |a|_min = ", ToString[N[aMin /. num /. ff -> fv, 6]]]],
  {rv, {2.05, 2.5, 4, 10, 100}}];
Print["      (finite throughout the exterior: for Ehat > A the freezing"];
Print["       surface is not reached, since f < 1 everywhere outside)"];
check["for Ehat > A the exterior thrust is bounded",
      N[aMin /. num /. ff -> 999/1000, 20] < 10];

Print["=========================================================="];
Print[If[$Failures === 0, "ALL CHECKS PASSED",
     ToString[$Failures] <> " CHECK(S) FAILED"]];
Print["=========================================================="];
If[$Failures > 0, Exit[1], Exit[0]];
