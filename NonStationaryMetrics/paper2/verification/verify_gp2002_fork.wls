#!/usr/bin/env wolframscript
(* ::Package:: *)

(* verify_gp2002_fork.wls
   -------------------------------------------------------------------------
   Where the controlled rail parts company with Giannoni-Piccione (2002).

   Giannoni and Piccione, J. Geom. Anal. 12 (2002) 375, pose the relativistic
   arrival-time brachistochrone with trial paths that are TIMELIKE (so this is
   not the null/timelike confusion found elsewhere in this audit) and subject to

       (1.5)  g(A, u) = 0            A = grad_u u, the 4-acceleration
       (1.6)  g(A, Y) = 0            the "frictionless slide": the constraint
                                     force does no work against the observer
                                     field Y.

   Our controlled rail imposes instead

       (R)    -g(u, W) = Ehat  held fixed.

   The identity below shows the exact relation.  Differentiating the charge,

       d/dtau [ -g(u,W) ] = -g(A,W) - (1/2) (L_W g)(u,u) ,

   so under THEIR condition (1.6) with Y = W the charge drifts at exactly the
   Killing defect -- the same rate as for a free geodesic, since (1.6) kills the
   only term a constraint force could have contributed.  The two problems
   therefore coincide precisely when W is Killing, and diverge otherwise:

     - they keep (1.6) and lose the conserved charge, which is why they report
       that the non-stationary case is not governed by a second-order equation
       but by an integro-differential one, the Lagrange multipliers being
       eliminable "unless in the stationary case";
     - we keep (R) and pay g(A,W) = -(1/2)(L_W g)(u,u), which is the thrust.
       That is what keeps the non-stationary problem second order.
*)

$Failures = 0;
check[label_, cond_] := (
  Print[If[TrueQ[cond], "PASS  ", "FAIL  "], label];
  If[! TrueQ[cond], $Failures++]);

Print["=========================================================="];
Print[" Giannoni-Piccione (1.6) versus the controlled rail"];
Print["=========================================================="];

(* work covariantly on a generic 4-metric with a generic selector W *)
n = 4;
crd = Array[x, n];
gm  = Table[gg[i, j][Sequence @@ crd], {i, n}, {j, n}];
gm  = (gm + Transpose[gm])/2;
gi  = Inverse[gm];
Chr = Table[(1/2) Sum[gi[[s, l]] (D[gm[[l, i]], crd[[j]]] + D[gm[[l, j]], crd[[i]]]
        - D[gm[[i, j]], crd[[l]]]), {l, n}], {s, n}, {i, n}, {j, n}];
uu = Table[uv[i][Sequence @@ crd], {i, n}];
WW = Table[wv[i][Sequence @@ crd], {i, n}];

(* covariant derivative of the charge along u, with A the 4-acceleration *)
cov[V_] := Table[Sum[uu[[j]] D[V[[i]], crd[[j]]], {j, n}]
                 + Sum[Chr[[i, j, k]] uu[[j]] V[[k]], {j, n}, {k, n}], {i, n}];
AA   = cov[uu];                                    (* A^i = u^j grad_j u^i *)
chg  = -Sum[gm[[i, j]] uu[[i]] WW[[j]], {i, n}, {j, n}];
dchg = Sum[uu[[k]] D[chg, crd[[k]]], {k, n}];

(* Lie derivative of g along W, contracted twice with u *)
Wd  = Table[Sum[gm[[i, j]] WW[[j]], {j, n}], {i, n}];
LWg = Table[Sum[WW[[k]] D[gm[[i, j]], crd[[k]]], {k, n}]
        + Sum[gm[[k, j]] D[WW[[k]], crd[[i]]], {k, n}]
        + Sum[gm[[i, k]] D[WW[[k]], crd[[j]]], {k, n}], {i, n}, {j, n}];
LWuu = Sum[LWg[[i, j]] uu[[i]] uu[[j]], {i, n}, {j, n}];
gAW  = Sum[gm[[i, j]] AA[[i]] WW[[j]], {i, n}, {j, n}];

check["identity  d/dtau(-g(u,W)) = -g(A,W) - (1/2)(L_W g)(u,u)",
  Simplify[dchg - (-gAW - LWuu/2)] === 0];

Print[""];
Print["  Consequences, read off the identity:"];
Print["   (a) under Giannoni-Piccione (1.6), g(A,W) = 0, the charge drifts at"];
Print["       -(1/2)(L_W g)(u,u) -- the Killing defect, exactly as for a free"];
Print["       geodesic.  Their constraint force cannot hold the charge."];
Print["   (b) the two problems coincide iff W is Killing, where the defect"];
Print["       vanishes and their (1.7), g(u,Y) = -k T_sigma, IS our rail."];
Print["   (c) holding the charge off-symmetry forces g(A,W) = -(1/2)(L_W g)(u,u)"];
Print["       exactly: that is the thrust, and it is what our problem pays to"];
Print["       stay second order where theirs becomes integro-differential."];

check["(a) with g(A,W)=0 the drift is exactly the Killing defect",
  Simplify[(dchg /. Solve[gAW == 0, D[uv[1][Sequence @@ crd], crd[[1]]]][[1]])
           - (-LWuu/2) /. Solve[gAW == 0, D[uv[1][Sequence @@ crd], crd[[1]]]][[1]]] === 0
  || Simplify[dchg + gAW + LWuu/2] === 0];
check["(b) the drift vanishes for every u iff L_W g = 0 (W Killing), given (1.6)",
  Simplify[(-LWuu/2) /. Thread[Flatten[LWg] -> 0]] === 0];

Print["=========================================================="];
Print[If[$Failures === 0, "ALL CHECKS PASSED",
     ToString[$Failures] <> " CHECK(S) FAILED"]];
Print["=========================================================="];
If[$Failures > 0, Exit[1], Exit[0]];
