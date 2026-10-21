#!/usr/bin/env wolframscript
(* verify_tk_komar.wls
   Komar angular momentum of Thakurta-Kerr, g = A(eta)^2 g_Kerr, for the exact
   axial Killing field xi = d_phi, on the Boyer-Lindquist spheres {eta, r} = const.
   Only sqrt(-g) F^{eta r} (F = d xi_flat) enters the pullback of the dual 2-form;
   the A' terms of F sit in F_{eta phi}. Claim: J = A(eta)^2 M a for every r.
   Run: wolframscript -file verify_tk_komar.wls *)
$Failures = 0;
check[label_, cond_] := (If[! TrueQ[cond], $Failures++]; Print[If[TrueQ[cond], "PASS  ", "FAIL  "], label]);
xs = {t, r, th, ph};
SS = r^2 + a^2 Cos[th]^2; DD = r^2 - 2 M r + a^2;
gK = {{-(1 - 2 M r/SS), 0, 0, -2 M a r Sin[th]^2/SS}, {0, SS/DD, 0, 0}, {0, 0, SS, 0},
      {-2 M a r Sin[th]^2/SS, 0, 0, (r^2 + a^2 + 2 M a^2 r Sin[th]^2/SS) Sin[th]^2}};
dens[g_] := Module[{gi = Simplify[Inverse[g]], xi, F, Fup},
  xi = g[[All, 4]];
  F = Table[D[xi[[n]], xs[[m]]] - D[xi[[m]], xs[[n]]], {m, 4}, {n, 4}];
  Fup = gi.F.gi;
  Simplify[Sqrt[-Simplify[Det[g]]] Fup[[1, 2]], {r > 2 M, M > 0, 0 < th < Pi, Aa[t] > 0}]];
dK = dens[gK]; dT = dens[Aa[t]^2 gK];
check["density on {eta,r} spheres: TK = A(eta)^2 x Kerr, with no A' terms",
  Simplify[dT - Aa[t]^2 dK] === 0 && FreeQ[Simplify[dT], Derivative[1][Aa]]];
JK = Simplify[2 Pi Integrate[dK, {th, 0, Pi}]/(16 Pi), {r > 2 M, M > 0, a > 0}];
check["Kerr: J = +M a at every radius (orientation of the sphere fixed so that J_Kerr = +M a)", Simplify[JK - M a] === 0];
Print[If[$Failures === 0, "ALL CHECKS PASSED", ToString[$Failures] <> " CHECK(S) FAILED"]];
If[$Failures =!= 0, Exit[1]];
