#!/usr/bin/env wolframscript
(* ::Package:: *)

(* verify_optical_planes_umbilicity.wls
   -------------------------------------------------------------------------
   Symbolic check, for general r, M, iota and s, of the Riemannian statements
   of the appendix paragraph "Central planes of the optical data"
   (app:rail-surfaces) and of eq:inclined-plane-force.

   (A) Conformal reduction: alpha_K^2 = Psi^2 hbar with
       Psi^2 = Sigma^2/(Delta Xi), Xi = Sigma - 2Mr, and
       hbar = dr^2 + Delta dtheta^2 + Delta^2 sin^2(theta) dphi^2/Xi, exactly in a.
       h_t = n^2 alpha_K^2 and h_tau = (A^2 N^2/Ehat)^2 h_t are then conformal
       to hbar, and ordinary umbilicity is decided for hbar.
   (B) Equatorial and meridional coordinate planes: the second fundamental form
       of hbar vanishes on their tangent vectors, exactly in a.
   (C) Central plane of inclination 0 < iota < pi/2, as the level set
       S = sin(iota) sin(theta) cos(phi) + cos(iota) cos(theta) = 0, with
       II = Hess_hbar(S)/|dS|_hbar on the tangents d_r and d_s, where the point
       of the plane is r (cos s e1 + sin s e2), e1 = (cos iota, 0, -sin iota),
       e2 = (0, 1, 0). The coefficients of a^0, a^1 and a^2 are taken by exact
       differentiation in a:
         a^0 and a^1: II = 0 (totally geodesic through first order);
         induced metric at a = 0: diag(1, B), B = r(r - 2M);
         a^2: the two minors of eq:inclined-plane-umbilicity,
           -(r - M) cos(iota) sin(iota) sin(s)/B^(3/2),  -sin(2 iota) cos(s)/B^(1/2),
         which never vanish together (Reduce).
       Mutation: (r - 2M) in place of (r - M) is rejected.
   (D) eq:inclined-plane-force: d beta_K(T, m) at first order in a, from
       beta_K = -2 M a sin^2(theta) dphi/(r - 2M), in the Cartesian chart built
       on Boyer-Lindquist coordinates.

   verify_optical_planes.py checks the same minors at sampled points (40 digits)
   and by an exact-metric Cartesian route; this script proves them symbolically.
   Run: wolframscript -file verify_optical_planes_umbilicity.wls
   ------------------------------------------------------------------------- *)

$Failures = 0;
check[label_, cond_] := (If[! TrueQ[cond], $Failures++];
  Print[If[TrueQ[cond], "PASS  ", "FAIL  "], label]);

xs = {r, th, ph};
SS = r^2 + a^2 Cos[th]^2; DD = r^2 - 2 M r + a^2; Xi = SS - 2 M r;
asmR = r > 2 M && M > 0 && 0 < th < Pi;

Print["--- (A) conformal reduction of the optical metric"];
gtt = -(1 - 2 M r/SS); gtp = -2 M a r Sin[th]^2/SS;
gpp = (r^2 + a^2 + 2 M a^2 r Sin[th]^2/SS) Sin[th]^2;
alpha2 = DiagonalMatrix[{SS/DD, SS, gpp - gtp^2/gtt}]/(-gtt);
hb = DiagonalMatrix[{1, DD, DD^2 Sin[th]^2/Xi}];
check["alpha_K^2 = Psi^2 hbar, Psi^2 = Sigma^2/(Delta Xi), exactly in a",
  Simplify[alpha2 - SS^2/(DD Xi) hb, asmR] === ConstantArray[0, {3, 3}]];

hbi = Inverse[hb];
Gam = Table[Sum[hbi[[k, l]] (D[hb[[l, j]], xs[[i]]] + D[hb[[l, i]], xs[[j]]] - D[hb[[i, j]], xs[[l]]])/2, {l, 3}],
   {k, 3}, {i, 3}, {j, 3}];
hess[Sf_] := Table[D[Sf, xs[[i]], xs[[j]]] - Sum[Gam[[k, i, j]] D[Sf, xs[[k]]], {k, 3}], {i, 3}, {j, 3}];

Print["--- (B) equatorial and meridional planes, exactly in a"];
Heq = hess[Cos[th]];
check["equator (theta = pi/2): II on d_r, d_phi vanishes",
  Simplify[Table[u.Heq.w, {u, {{1, 0, 0}, {0, 0, 1}}}, {w, {{1, 0, 0}, {0, 0, 1}}}] /. th -> Pi/2] === {{0, 0}, {0, 0}}];
Hmer = hess[Sin[ph]];
check["meridional plane (phi = 0): II on d_r, d_theta vanishes",
  Simplify[Table[u.Hmer.w, {u, {{1, 0, 0}, {0, 1, 0}}}, {w, {{1, 0, 0}, {0, 1, 0}}}] /. ph -> 0] === {{0, 0}, {0, 0}}];

Print["--- (C) inclined central plane, orders a^0, a^1, a^2"];
Sf = Sin[io] Sin[th] Cos[ph] + Cos[io] Cos[th];
dS = D[Sf, {xs}]; H = hess[Sf]; nrm = Sqrt[dS.hbi.dS];
st = Sqrt[1 - Sin[io]^2 Cos[s]^2];
onPlane = {th -> ArcCos[-Sin[io] Cos[s]], ph -> ArcTan[Cos[io] Cos[s], Sin[s]]};
tR = {1, 0, 0}; tS = {0, -Sin[io] Sin[s]/st, Cos[io]/st^2};   (* d_r and d_s *)
asm = r > 2 M && M > 0 && 0 < io < Pi/2 && Element[s, Reals];
ord[e_, k_] := FullSimplify[((D[e, {a, k}] /. a -> 0)/k!) /. onPlane, asm];
check["the parametrised point lies on the plane, and d_r, d_s are tangent",
  FullSimplify[{Sf, dS.tR, dS.tS} /. onPlane, asm] === {0, 0, 0}];
II[k_] := Table[ord[(u.H.w)/nrm, k], {u, {tR, tS}}, {w, {tR, tS}}];
II0 = II[0]; II1 = II[1]; II2 = II[2];
g0 = Table[ord[u.hb.w, 0], {u, {tR, tS}}, {w, {tR, tS}}];
B = r (r - 2 M);
check["order a^0: II = 0", II0 === {{0, 0}, {0, 0}}];
check["order a^1: II = 0 (totally geodesic through first order in a)", II1 === {{0, 0}, {0, 0}}];
check["induced metric at a = 0: diag(1, B), B = r(r - 2M)", FullSimplify[g0 - {{1, 0}, {0, B}}, asm] === {{0, 0}, {0, 0}}];
m1 = FullSimplify[II2[[1, 2]] g0[[1, 1]] - II2[[1, 1]] g0[[1, 2]], asm];
m2 = FullSimplify[II2[[1, 1]] g0[[2, 2]] - II2[[2, 2]] g0[[1, 1]], asm];
m3 = FullSimplify[II2[[1, 2]] g0[[2, 2]] - II2[[2, 2]] g0[[1, 2]], asm];
t1 = -(r - M) Cos[io] Sin[io] Sin[s]/B^(3/2); t2 = -Sin[2 io] Cos[s]/Sqrt[B];
check["order a^2, first minor = -(r - M) cos(iota) sin(iota) sin(s)/B^(3/2)", FullSimplify[m1 - t1, asm] === 0];
check["order a^2, second minor = -sin(2 iota) cos(s)/B^(1/2)", FullSimplify[m2 - t2, asm] === 0];
check["order a^2, third minor = B x first minor (no further condition)", FullSimplify[m3 - B m1, asm] === 0];
check["the two minors never vanish together for r > 2M, 0 < iota < pi/2",
  Reduce[t1 == 0 && t2 == 0 && asm, {r, M, io, s}, Reals] === False];
check["mutation: (r - 2M) in place of (r - M) is rejected",
  FullSimplify[m1 + (r - 2 M) Cos[io] Sin[io] Sin[s]/B^(3/2), asm] =!= 0];

Print["--- (D) normal magnetic force on the inclined plane, first order in a"];
b1 = Simplify[Normal[Series[-gtp/gtt, {a, 0, 1}]]];
check["beta_K = -2 M a sin^2(theta) dphi/(r - 2M) + O(a^2)", Simplify[b1 + 2 M a Sin[th]^2/(r - 2 M)] === 0];
cs = {X, Y, Z};
cc = Simplify[b1/(r^2 Sin[th]^2)] /. r -> Sqrt[X^2 + Y^2 + Z^2];   (* beta = cc (X dY - Y dX) *)
bet = {-cc Y, cc X, 0};
Fm = Table[D[bet[[j]], cs[[i]]] - D[bet[[i]], cs[[j]]], {i, 3}, {j, 3}];
e1 = {Cos[io], 0, -Sin[io]}; e2 = {0, 1, 0}; m = Cross[e1, e2];
pt = rho (Cos[s] e1 + Sin[s] e2); T = -Sin[s] e1 + Cos[s] e2;
force = Simplify[T.(Fm /. Thread[cs -> pt]).m, rho > 2 M && M > 0];
check["d beta_K(T, m) = 4 M a cos(s) sin(iota)/((rho - 2M) rho^2)",
  Simplify[force - 4 M a Cos[s] Sin[io]/((rho - 2 M) rho^2), rho > 2 M] === 0];

Print[If[$Failures === 0, "ALL CHECKS PASSED", ToString[$Failures] <> " CHECK(S) FAILED"]];
If[$Failures =!= 0, Exit[1]];
