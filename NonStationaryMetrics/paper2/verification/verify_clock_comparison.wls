#!/usr/bin/env wolframscript
(* ::Package:: *)

(* verify_clock_comparison.wls
   -------------------------------------------------------------------------
   The scalar comparison lemma behind the saturation step of Paper I's
   existence theorem (Lemma "scalar comparison for the clock").

   The relaxed problem replaces the clock equation by the inclusion
   chi' >= F(x, xdot, chi).  Recovering an admissible arc requires that the
   saturated solution, chi' = F, has the SMALLEST terminal clock among all
   such trajectories.  That does not follow from monotonicity in t alone,
   because F depends on chi as well; it needs the comparison argument, whose
   only hypothesis is the Lipschitz bound in the slow argument already
   established by the Gronwall step.

   Here the lemma is checked numerically on families of F that satisfy the
   hypothesis, including cases where F DECREASES in chi -- the case in which
   a naive monotonicity argument would fail.
*)

$Failures = 0;
check[l_, c_] := (Print[If[TrueQ[c], "PASS  ", "FAIL  "], l]; If[! TrueQ[c], $Failures++]);

Print["=========================================================="];
Print[" comparison lemma: saturated clock is the smallest"];
Print["=========================================================="];

(* a family of admissible F: positive, Lipschitz in the slow argument chi *)
Ftest = {
  {"F increasing in chi",  Function[{s, c}, 1 + 0.3 c]},
  {"F decreasing in chi",  Function[{s, c}, 2 - 0.3 c]},
  {"F oscillating in s",   Function[{s, c}, 1.5 + 0.4 Sin[6 s] + 0.2 c]},
  {"F nearly chi-free",    Function[{s, c}, 1.2 + 0.01 c]}
};

Do[
  lbl = f[[1]]; FF = f[[2]];
  (* saturated solution *)
  sat = NDSolveValue[{x'[s] == FF[s, x[s]], x[0] == 0}, x, {s, 0, 1}];
  (* several relaxed trajectories: chi' = F + g(s) with g >= 0 *)
  worst = Infinity; ok = True;
  Do[
    g = gg;
    rel = NDSolveValue[{y'[s] == FF[s, y[s]] + g[s], y[0] == 0}, y, {s, 0, 1}];
    d = rel[1] - sat[1];
    If[d < -1.*^-9, ok = False];
    worst = Min[worst, d],
    {gg, {Function[s, 0.], Function[s, 0.5], Function[s, 2. s],
          Function[s, 3. s^2 (1 - s)], Function[s, 0.8 (1 + Sin[10 s])^2]}}];
  check[lbl <> ": every relaxed clock >= saturated clock (min gap " <>
        ToString[NumberForm[worst, 4]] <> ")", ok],
  {f, Ftest}];

Print[""];
Print["      The decreasing case is the point: there the terminal clock is NOT"];
Print["      monotone in the pointwise value of chi', so the naive argument"];
Print["      fails while the comparison lemma still gives the conclusion."];

Print["=========================================================="];
Print[If[$Failures === 0, "ALL CHECKS PASSED",
     ToString[$Failures] <> " CHECK(S) FAILED"]];
If[$Failures > 0, Exit[1], Exit[0]];
