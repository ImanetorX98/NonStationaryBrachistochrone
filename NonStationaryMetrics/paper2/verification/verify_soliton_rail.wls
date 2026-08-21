#!/usr/bin/env wolframscript
(* ::Package:: *)

(* verify_soliton_rail.wls
   -------------------------------------------------------------------------
   The rail-drift identity for a Ricci-soliton selector.

   A Ricci soliton is a vector field X with prescribed deformation tensor,
        (1/2) L_X g + Ric = lambda g .
   The controlled rail drifts along a GEODESIC at
        d/dtau ( -g(u,W) ) = -(1/2) (L_W g)_{mu nu} u^mu u^nu
   (verify_submersion_link.wls).  Combining the two, for W = X a soliton field:

        d Ehat / d tau |_geodesic = lambda + Ric(u,u).                    [star]

   So for a soliton selector the control effort needed to hold the rail is
   PURE CURVATURE: no derivative of W survives.  By the Einstein equations
   Ric(u,u) is the Raychaudhuri focusing scalar, i.e. the energy density the
   traveller measures.  Two corollaries:

     - vacuum + steady (Ric = 0, lambda = 0)  =>  X is Killing and the rail is
       conserved with no control: the stationary case.
     - the rail needs control exactly when lambda + Ric(u,u) != 0.

   The script also tests whether the Kodama selector of ingoing Vaidya is a
   soliton field.  It is not, and the script exhibits the obstruction.
*)

$Failures = 0;
check[label_, cond_] := (
  Print[If[TrueQ[cond], "PASS  ", "FAIL  "], label];
  If[! TrueQ[cond], $Failures++]);

Print["=========================================================="];
Print[" rail drift for a Ricci-soliton selector"];
Print["=========================================================="];

(* ---------------------------------------------------------------- *)
(* 1. the abstract identity, done with symbolic tensors              *)
(* ---------------------------------------------------------------- *)
Print["--- the identity  dEhat/dtau = lambda + Ric(u,u) ---"];
n = 4;
gg = Table[Symbol["g" <> ToString[i] <> ToString[j]], {i, n}, {j, n}];
gg = (gg + Transpose[gg])/2;
ric = Table[Symbol["R" <> ToString[i] <> ToString[j]], {i, n}, {j, n}];
ric = (ric + Transpose[ric])/2;
uu = Table[Symbol["u" <> ToString[i]], {i, n}];
lam = Symbol["lambda"];

(* soliton condition in the standard normalisation:  (1/2) L_X g + Ric = lam g *)
LXg = 2 (lam gg - ric);
drift = Simplify[-(1/2) uu . LXg . uu];
Print["      drift = ", Simplify[drift]];
(* impose the unit-timelike normalisation g(u,u) = -1 *)
driftOnShell = Simplify[drift /. Solve[uu . gg . uu == -1,
    Symbol["g11"]][[1]]];
target = lam + uu . ric . uu;
check["W soliton  =>  dEhat/dtau = lambda + Ric(u,u)",
      Simplify[driftOnShell - target] === 0];

(* ric was symmetrised, so Flatten[ric] holds (Rij+Rji)/2 rather than bare
   symbols; build the vacuum substitution from the symbol names instead. *)
ricZero = Flatten@Table[
    Symbol["R" <> ToString[i] <> ToString[j]] -> 0, {i, n}, {j, n}];
check["vacuum and steady (Ric = 0, lambda = 0) => drift vanishes",
      Simplify[driftOnShell /. ricZero /. lam -> 0] === 0];
check["vacuum and steady => L_X g = 0, i.e. X is Killing",
      Simplify[LXg /. ricZero /. lam -> 0] === ConstantArray[0, {n, n}]];
check["in vacuum alone the drift is the soliton constant lambda",
      Simplify[driftOnShell /. ricZero] === lam];

(* ---------------------------------------------------------------- *)
(* 2. ingoing Vaidya: Ricci tensor and the Kodama selector           *)
(* ---------------------------------------------------------------- *)
Print["--- ingoing Vaidya: is the Kodama selector a soliton field? ---"];
ClearAll[v, r, th, ph, m, M];
coords = {v, r, th, ph};
F = 1 - 2 m[v]/r;
gV = {{-F, 1, 0, 0}, {1, 0, 0, 0}, {0, 0, r^2, 0},
      {0, 0, 0, r^2 Sin[th]^2}};
gInv = Simplify[Inverse[gV]];

Chr = Table[
   Simplify[(1/2) Sum[gInv[[a, d]] (D[gV[[d, b]], coords[[c]]] +
        D[gV[[d, c]], coords[[b]]] - D[gV[[b, c]], coords[[d]]]), {d, 4}]],
   {a, 4}, {b, 4}, {c, 4}];

Riem = Table[
   Simplify[D[Chr[[a, b, d]], coords[[c]]] - D[Chr[[a, b, c]], coords[[d]]] +
     Sum[Chr[[a, c, e]] Chr[[e, b, d]] - Chr[[a, d, e]] Chr[[e, b, c]], {e, 4}]],
   {a, 4}, {b, 4}, {c, 4}, {d, 4}];
Ric = Table[Simplify[Sum[Riem[[a, b, a, c]], {a, 4}]], {b, 4}, {c, 4}];
Print["      Ric components that do not vanish:"];
Do[If[Simplify[Ric[[i, j]]] =!= 0,
   Print["        R_", coords[[i]], coords[[j]], " = ", Simplify[Ric[[i, j]]]]],
  {i, 4}, {j, i, 4}];

check["Vaidya Ricci has the single component R_vv = 2 m'(v)/r^2",
      Simplify[Ric[[1, 1]] - 2 m'[v]/r^2] === 0 &&
      Simplify[Ric - DiagonalMatrix[{Ric[[1, 1]], 0, 0, 0}]]
        === ConstantArray[0, {4, 4}]];
check["Vaidya is Ricci-flat exactly in the stationary case",
      Simplify[Ric /. m -> Function[x, M]] === ConstantArray[0, {4, 4}]];

(* the Kodama selector and its deformation tensor *)
LKg = Simplify[D[gV, v]];
check["(L_K g) has the single component (vv) = 2 m'(v)/r",
      Simplify[LKg[[1, 1]] - 2 m'[v]/r] === 0];

(* soliton test: is (1/2) L_K g + Ric = lambda g for some scalar lambda? *)
solTensor = Simplify[(1/2) LKg + Ric];
Print["      (1/2) L_K g + Ric  has (vv) component = ",
      Simplify[solTensor[[1, 1]]]];
(* for this to be lambda*g we would need it proportional to g; test the
   angular block, which is identically zero while g_thth = r^2 != 0        *)
check["Kodama is NOT a Ricci-soliton field for Vaidya (angular block fails)",
      Simplify[solTensor[[3, 3]]] === 0 && Simplify[gV[[3, 3]]] =!= 0 &&
      Simplify[solTensor[[1, 1]]] =!= 0];
Print["      => the angular part would force lambda = 0, but then the (vv)"];
Print["         part gives m'(v)/r + 2 m'(v)/r^2 != 0.  No lambda works."];

(* the two rates, side by side *)
uv = Symbol["uv"];
driftVaidya = Simplify[-(1/2) LKg[[1, 1]] uv^2];
ricUU = Simplify[Ric[[1, 1]] uv^2];
Print["      geodesic rail drift    = ", driftVaidya];
Print["      Ric(u,u)               = ", ricUU];
check["for Vaidya the drift equals -(r/2) Ric(u,u)",
      Simplify[driftVaidya + (r/2) ricUU] === 0];

Print["=========================================================="];
Print[If[$Failures === 0, "ALL CHECKS PASSED",
     ToString[$Failures] <> " CHECK(S) FAILED"]];
Print["=========================================================="];
If[$Failures > 0, Exit[1], Exit[0]];
