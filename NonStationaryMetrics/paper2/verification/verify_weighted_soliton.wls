#!/usr/bin/env wolframscript
(* ::Package:: *)

(* verify_weighted_soliton.wls
   -------------------------------------------------------------------------
   The conformal-Killing selector of Thakurta-Kerr as a WEIGHTED soliton.

   g = A(eta)^2 gbar, gbar = Kerr (Boyer-Lindquist, eta = t), U = ln A.
   Unweighted, W = d_eta cannot be a Ricci-soliton field: that would need
   g Einstein, and Ric_{eta r} != 0 where A' != 0 (verify_tk_ckv_rail.wls).
   With the weight phi = 2U and the Bakry-Emery tensor

        Ric_phi^m = Ric + Hess phi - (1/m) dphi (x) dphi ,   m = 2 - n = -2,

   the conformal transformation of the Ricci tensor gives

        Ric_phi^m = mu g ,   mu = -(Box_gbar U)/A^2 = -gbar^{eta eta} U''/A^2    [1]

   (Kerr is Ricci-flat and its Boyer-Lindquist time is harmonic).  Since
   L_W g = 2 psi g with psi = A'/A,

        (1/2) L_W g + Ric_phi^m = lambda g ,   lambda = psi + mu .                [2]

   W is therefore a weighted soliton field with a FUNCTION lambda (an "almost"
   soliton), and a weighted soliton with constant lambda exactly when
   U'' = 0, i.e. A = A0 exp(c eta); then lambda = c.  Along a geodesic the rail
   drifts at  dEhat/dtau = lambda + Ric_phi^m(u,u) = psi,  the weighted form of
   eq. (soliton-drift).

   The tensor is the generalised m-Bakry-Emery tensor, defined algebraically
   and taken at the negative parameter m = -2; it is not the conformal
   rescaling of the metric used in the cited soliton papers.  For M > 0,
   lambda constant on the spacetime forces U'' = 0 because gbar^{eta eta}
   depends on r (checked on the equator below); the flat seed is excluded.

   Checks: [1] and [2] from the full Ricci tensor of g = A(eta)^2 gbar, for
   Thakurta-Kerr (gbar = Kerr) and for its a = 0 case (gbar = Schwarzschild); lambda constant iff A is exponential; the drift; and
   a control that must fail (no weight: Ric is not pure trace).
*)

$Failures = 0;
check[label_, cond_] := (
  Print[If[TrueQ[cond], "PASS  ", "FAIL  "], label];
  If[! TrueQ[cond], $Failures++]);

Print["=========================================================="];
Print[" Thakurta-Kerr selector as a weighted (Bakry-Emery) soliton"];
Print["=========================================================="];

ClearAll[eta, r, th, ph, M, a, A, c, A0];
coords = {eta, r, th, ph};
n = 4;

ricci[g_] := Module[{gi, chr, riem},
  gi = Simplify[Inverse[g]];
  chr = Table[Simplify[(1/2) Sum[gi[[i, d]] (D[g[[d, j]], coords[[k]]] +
        D[g[[d, k]], coords[[j]]] - D[g[[j, k]], coords[[d]]]), {d, n}]],
     {i, n}, {j, n}, {k, n}];
  {chr, Table[Simplify[Sum[D[chr[[i, j, k]], coords[[i]]] -
        D[chr[[i, j, i]], coords[[k]]] +
        Sum[chr[[i, i, e]] chr[[e, j, k]] - chr[[i, k, e]] chr[[e, j, i]], {e, n}],
       {i, n}]], {j, n}, {k, n}]}];

hess[chr_, s_] := Table[D[s, coords[[i]], coords[[j]]] -
    Sum[chr[[k, i, j]] D[s, coords[[k]]], {k, n}], {i, n}, {j, n}];

runTK[name_, gbar_] := Module[{g, chr, ric, U, phi, m, dphi, ricW, mu, psi,
    LWg, lam, gbarInv, rest, ricChk},
  Print["--- ", name, ": g = A(eta)^2 gbar ---"];
  g = A[eta]^2 gbar;
  {chr, ric} = ricci[g];
  U = Log[A[eta]]; phi = 2 U; m = 2 - n;
  dphi = D[phi, #] & /@ coords;
  ricW = Simplify[ric + hess[chr, phi] - (1/m) Outer[Times, dphi, dphi]];
  gbarInv = Simplify[Inverse[gbar]];
  mu = Simplify[-gbarInv[[1, 1]] D[U, {eta, 2}]/A[eta]^2];
  check[name <> ": its seed gbar is Ricci-flat",
    Simplify[ricci[gbar][[2]]] === ConstantArray[0, {n, n}]];
  check[name <> ": Box_gbar(U) = gbar^{eta eta} U''  (harmonic time)",
    Simplify[1/Sqrt[-Det[gbar]] Sum[D[Sqrt[-Det[gbar]] gbarInv[[i, j]] D[U, coords[[j]]],
        coords[[i]]], {i, n}, {j, n}] - gbarInv[[1, 1]] D[U, {eta, 2}]] === 0];
  check[name <> ": Ric_phi^m = mu g with phi = 2 ln A, m = -2   [1]",
    Simplify[ricW - mu g] === ConstantArray[0, {n, n}]];
  (* control that must fail: without the weight Ric is not pure trace *)
  ricChk = Simplify[ric - (ric[[3, 3]]/g[[3, 3]]) g];
  check[name <> ": CONTROL - without the weight Ric is not pure trace",
    ricChk =!= ConstantArray[0, {n, n}] && Simplify[ric[[1, 2]]] =!= 0];
  psi = D[A[eta], eta]/A[eta];
  LWg = D[g, eta];
  check[name <> ": L_W g = 2 psi g, psi = A'/A (conformal Killing)",
    Simplify[LWg - 2 psi g] === ConstantArray[0, {n, n}]];
  lam = Simplify[psi + mu];
  check[name <> ": (1/2) L_W g + Ric_phi^m = lambda g, lambda = psi + mu   [2]",
    Simplify[(1/2) LWg + ricW - lam g] === ConstantArray[0, {n, n}]];
  check[name <> ": lambda depends on r for generic A (A = 1 + eta^2 sampled)",
    Simplify[D[lam /. A -> Function[x, 1 + x^2], r]] =!= 0];
  check[name <> ": A = A0 exp(c eta)  =>  lambda = c, constant",
    Simplify[(lam /. A -> Function[x, A0 Exp[c x]]) - c] === 0];
  (* geodesic drift along a unit timelike u:  -(1/2) L_W g(u,u) = lambda + Ric_phi(u,u) *)
  check[name <> ": geodesic drift lambda + Ric_phi^m(u,u) = psi for g(u,u) = -1",
    Simplify[lam + mu (-1) - psi] === 0];
  mu];

gS = DiagonalMatrix[{-(1 - 2 M/r), 1/(1 - 2 M/r), r^2, r^2 Sin[th]^2}];
muS = runTK["Thakurta (a = 0)", gS];
check["Thakurta (a = 0): mu = (ln A)''/(f A^2)",
  Simplify[muS - D[Log[A[eta]], {eta, 2}]/((1 - 2 M/r) A[eta]^2)] === 0];

Sig = r^2 + a^2 Cos[th]^2; Del = r^2 - 2 M r + a^2;
gK = {{-(1 - 2 M r/Sig), 0, 0, -2 M a r Sin[th]^2/Sig},
      {0, Sig/Del, 0, 0},
      {0, 0, Sig, 0},
      {-2 M a r Sin[th]^2/Sig, 0, 0,
       (r^2 + a^2 + 2 M a^2 r Sin[th]^2/Sig) Sin[th]^2}};
runTK["Thakurta-Kerr", gK];

(* necessity of U'' = 0 for constant lambda, M > 0: on the equator gbar^{eta eta}
   is not constant in r, so lambda at two radii differs unless U'' = 0 *)
gtt = Simplify[Inverse[gK][[1, 1]] /. th -> Pi/2];
check["Thakurta-Kerr: equatorial gbar^{eta eta} = -1 - 2M(r^2+a^2)/(r Delta)",
  Simplify[gtt - (-1 - 2 M (r^2 + a^2)/(r (r^2 - 2 M r + a^2)))] === 0];
check["Thakurta-Kerr: d/dr gbar^{eta eta} != 0 on the equator for M > 0",
  Simplify[D[gtt, r] /. {M -> 1, a -> 1/2, r -> 3}] != 0];
check["flat seed (M = 0): gbar^{eta eta} = -1, the necessity argument does not apply",
  Simplify[gtt /. M -> 0] === -1];

Print[];
If[$Failures == 0, Print["ALL CHECKS PASSED"], Print[$Failures, " CHECK(S) FAILED"]];
Exit[If[$Failures == 0, 0, 1]];
