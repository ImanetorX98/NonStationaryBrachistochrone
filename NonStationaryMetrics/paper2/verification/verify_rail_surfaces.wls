#!/usr/bin/env wolframscript
(* ::Package:: *)

(* verify_rail_surfaces.wls
   -------------------------------------------------------------------------
   Extrinsic geometry of the timelike surfaces swept by three equatorial rail
   families in Thakurta-Kerr, g = A(eta)^2 g_Kerr, with A an ARBITRARY positive
   function of the conformal time eta (Lemma lem:rail-surfaces, appendix
   app:rail-surfaces).

   Everything is computed in four dimensions from the Christoffel symbols of the
   full metric; no intermediate formula of the text is assumed.

   (A) Kerr seed: surface Sigma(r,psi) = (T(r), r, pi/2, psi), q = T'(r) generic.
       Induced metric h0, det h0 = -r^2 S, polar normal component of II zero,
       b0_{r psi}, b0_{psi psi}, and the umbilicity minor.
   (B) Conformal change on the same immersion: II_g = II_0 - h0 (grad0 ln A)^perp,
       h_g = A^2 h0; minor with the physical unit normal = A^3 x seed minor.
   (C) a = 0, outgoing (and ingoing) radial rails of charge Ehat with A(eta)
       varying along the rail: h, b in the physical unit normal; umbilicity iff
       d(A r)/d eta = 0; b_{psi psi} != 0 in that case; the A r = c profile.
   (D) Common launch: leading algebra in an orthonormal frame at the launch.
   (E) The tangential component of the magnetic variation field (introduction):
       eta' = -B xi, and the normal component obeys xi'' + (K + B^2 + dB(N)) xi = 0.

   An independent SymPy route is verify_rail_surfaces_sympy.py.
   Run: wolframscript -file verify_rail_surfaces.wls
   ------------------------------------------------------------------------- *)

$Failures = 0;
check[label_, cond_] := (If[! TrueQ[cond], $Failures++];
  Print[If[TrueQ[cond], "PASS  ", "FAIL  "], label]);
zr[e_] := TrueQ[Simplify[Together[e]] === 0];

xs = {t, r, th, ph};
SS = r^2 + a^2 Cos[th]^2; DD = r^2 - 2 M r + a^2;
gK = {{-(1 - 2 M r/SS), 0, 0, -2 M a r Sin[th]^2/SS}, {0, SS/DD, 0, 0}, {0, 0, SS, 0},
      {-2 M a r Sin[th]^2/SS, 0, 0, (r^2 + a^2 + 2 M a^2 r Sin[th]^2/SS) Sin[th]^2}};
christ[g_] := Module[{gi = Simplify[Inverse[g]]},
  Table[Simplify[Sum[gi[[i, l]] (D[g[[l, j]], xs[[k]]] + D[g[[l, k]], xs[[j]]] - D[g[[j, k]], xs[[l]]])/2, {l, 4}]],
    {i, 4}, {j, 4}, {k, 4}]];
GK = christ[gK]; gT = Aa[t]^2 gK; GT = christ[gT];
eq = {th -> Pi/2};
X = {q[r], 1, 0, 0}; Y = {0, 0, 0, 1}; nc = {-1, q[r], 0, 0};
cov[G_, U_, V_, dV_] := Table[dV[[i]] + Sum[G[[i, j, k]] U[[j]] V[[k]], {j, 4}, {k, 4}], {i, 4}];
(* raw coefficients c_ab = N_mu (nabla_a e_b)^mu with the unnormalised covector N = (-1,q,0,0) *)
surf[g_, G_] := Module[{ge = g /. eq, Ge = G /. eq, vecs},
  vecs = {cov[Ge, X, X, {q'[r], 0, 0, 0}], cov[Ge, X, Y, {0, 0, 0, 0}], cov[Ge, Y, Y, {0, 0, 0, 0}]};
  {Simplify[nc.Inverse[ge].nc], Simplify[{{X.ge.X, X.ge.Y}, {Y.ge.X, Y.ge.Y}}],
   Simplify[nc.# & /@ vecs], Simplify[{0, 0, 1, 0}.# & /@ vecs]}];
sub = {Aa[t] -> AA, Aa'[t] -> A1, t -> T[r]};   (* A1 = dA/d eta, evaluated on the surface *)

f = 1 - 2 M/r; k = -2 M a/r; p = r^2 + a^2 + 2 M a^2/r; Del = r^2 - 2 M r + a^2;
Sv = Del q[r]^2/r^2 - p/Del;

Print["--- (A) Kerr seed"];
{nnK, hK, cK, polK} = surf[gK, GK];
check["g0^-1(N,N) = S = Delta q^2/r^2 - p/Delta", zr[nnK - Sv]];
check["h0 = [[r^2/Delta - f q^2, k q], [k q, p]], det h0 = -r^2 S",
  (And @@ (zr /@ Flatten[hK - {{r^2/Del - f q[r]^2, k q[r]}, {k q[r], p}}])) && zr[Det[hK] + r^2 Sv]];
check["polar normal component of II vanishes (reflection theta -> pi - theta)", polK === {0, 0, 0}];
check["sqrt(S) b0_{r psi} = (p k' - k p')/(2 Delta) - q^2 Delta k'/(2 r^2)",
  zr[cK[[2]] - ((p D[k, r] - k D[p, r])/(2 Del) - q[r]^2 Del D[k, r]/(2 r^2))]];
check["sqrt(S) b0_{psi psi} = -q Delta p'/(2 r^2)", zr[cK[[3]] + q[r] Del D[p, r]/(2 r^2)]];
check["p k' - k p' = 2 M a (3 r^2 + a^2)/r^2 and f p + k^2 = Delta",
  zr[p D[k, r] - k D[p, r] - 2 M a (3 r^2 + a^2)/r^2] && zr[f p + k^2 - Del]];
check["raw seed minor c_{r psi} h_{psi psi} - c_{psi psi} h_{r psi} = -M a (3 r^2 + a^2) S/r^2",
  zr[cK[[2]] hK[[2, 2]] - cK[[3]] hK[[1, 2]] + M a (3 r^2 + a^2) Sv/r^2]];

Print["--- (B) Thakurta-Kerr, same immersion, arbitrary A(eta)"];
{nnT, hT, cT, polT} = surf[gT, GT] /. sub;
check["polar normal component of II vanishes for every A", Simplify[polT] === {0, 0, 0}];
check["h_g = A^2 h0 and g^-1(N,N) = S/A^2", (And @@ (zr /@ Flatten[hT - AA^2 hK])) && zr[nnT - Sv/AA^2]];
nupK = Inverse[gK /. eq].nc;
check["raw: c_g = c0 - (A'/A) N^eta h0  (II_g = II_0 - h0 (grad0 ln A)^perp)",
  And @@ (zr /@ (cT - cK + (A1/AA) nupK[[1]] {hK[[1, 1]], hK[[1, 2]], hK[[2, 2]]}))];
minors[c_, h_] := {c[[2]] h[[2, 2]] - c[[3]] h[[1, 2]], c[[1]] h[[2, 2]] - c[[3]] h[[1, 1]], c[[1]] h[[1, 2]] - c[[2]] h[[1, 1]]};
check["trace-free part unchanged: all raw minors scale by A^2", And @@ (zr /@ (minors[cT, hT] - AA^2 minors[cK, hK]))];
(* physical unit normal covector A N/sqrt(S): b = (A/sqrt S) c. Multiply the claim by sqrt(S) to stay rational. *)
check["unit normal: b_{r psi} h_{psi psi} - b_{psi psi} h_{r psi} = -A^3 M a (3 r^2 + a^2) sqrt(S)/r^2",
  zr[AA (cT[[2]] hT[[2, 2]] - cT[[3]] hT[[1, 2]]) + AA^3 M a (3 r^2 + a^2) Sv/r^2]];

Print["--- (C) a = 0, radial rails, A varying along the rail"];
Dv = EE^2 - AA^2 f; asmR = M > 0 && r > 2 M && AA > 0 && EE > 0 && Dv > 0;
(* charge Ehat = A^2 f etadot and mass shell: etadot = Ehat/(A^2 f), rdot^2 = D/A^4 *)
etad = EE/(AA^2 f);
check["rail: rdot^2 = (Ehat^2 - A^2 f)/A^4 from -A^2 f etadot^2 + A^2 rdot^2/f = -1",
  zr[(AA^2 f etad^2 - 1) f/AA^2 - Dv/AA^4]];
qr = EE/(f Sqrt[Dv]);                       (* q = d eta/dr on the outgoing branch *)
qpr = D[qr, r] + D[qr, AA] A1 qr;           (* total derivative: dA/dr = A_eta q *)
rad = Simplify[(Sqrt[Dv] cT /. a -> 0) /. {q[r] -> qr, q'[r] -> qpr}, asmR];   (* unit covector sqrt(D) N *)
hR = Simplify[(hT /. a -> 0) /. q[r] -> qr, asmR];
check["unit normal is sqrt(D)(-1,q,0,0): S/A^2 = 1/D", zr[Simplify[(Sv/AA^2 /. a -> 0) /. q[r] -> qr, asmR] - 1/Dv]];
check["h_rr = -A^4/D, h_{r psi} = 0, h_{psi psi} = A^2 r^2", zr[hR[[1, 1]] + AA^4/Dv] && zr[hR[[1, 2]]] && zr[hR[[2, 2]] - AA^2 r^2]];
check["b_rr = -A^3 A_eta/D^(3/2), b_{r psi} = 0", zr[Simplify[rad[[1]] + AA^3 A1/Dv^(3/2), asmR]] && zr[rad[[2]]]];
check["b_{psi psi} = -Ehat r - A_eta r^2 sqrt(D)/(A f)", zr[Simplify[rad[[3]] + EE r + A1 r^2 Sqrt[Dv]/(AA f), asmR]]];
minR = Simplify[rad[[1]] hR[[2, 2]] - rad[[3]] hR[[1, 1]], asmR];
check["diagonal minor = -(A^3 Ehat r/D^(3/2)) (A sqrt(D) + A_eta Ehat r/f)",
  zr[Simplify[minR + AA^3 EE r/Dv^(3/2) (AA Sqrt[Dv] + A1 EE r/f), asmR]]];
A1u = -AA f Sqrt[Dv]/(EE r);
check["minor vanishes iff A_eta = -(A/r) dr/deta, i.e. d(A r)/d eta = 0 (coefficient of A_eta nonzero)",
  zr[Simplify[minR - D[minR, A1] (A1 - A1u), asmR]] && zr[Simplify[D[minR, A1] + AA^3 EE^2 r^2/(Dv^(3/2) f), asmR]] &&
  zr[A1u + AA/(r qr)]];
check["umbilic case: b = H_n h with H_n = -f/(Ehat r), and b_{psi psi} = -A^2 f r/Ehat != 0",
  zr[Simplify[(rad[[1]] /. A1 -> A1u)/hR[[1, 1]] + f/(EE r), asmR]] &&
  zr[Simplify[(rad[[3]] /. A1 -> A1u)/hR[[2, 2]] + f/(EE r), asmR]] &&
  zr[Simplify[rad[[3]] /. A1 -> A1u, asmR] + AA^2 f r/EE]];
check["A_eta = 0 (frozen A): minor nonzero, b_rr = 0", zr[rad[[1]] /. A1 -> 0] &&
  Resolve[ForAll[{r, M, EE, AA}, asmR, (minR /. A1 -> 0) < 0], Reals]];
qin = -qr; qpin = D[qin, r] + D[qin, AA] A1 qin;
radin = Simplify[(Sqrt[Dv] cT /. a -> 0) /. {q[r] -> qin, q'[r] -> qpin}, asmR];
check["ingoing branch: minor vanishes iff d(A r)/d eta = 0 (A_eta = +A f sqrt(D)/(Ehat r))",
  zr[Simplify[(radin[[1]] hR[[2, 2]] - radin[[3]] hR[[1, 1]]) /. A1 -> -A1u, asmR]]];
check["A r = c profile: D = Ehat^2 - c^2 f/r^2; max_{r>2M} f/r^2 = 1/(27 M^2)",
  zr[(Dv /. AA -> c/r) - (EE^2 - c^2 f/r^2)] && zr[(f/r^2 /. r -> 3 M) - 1/(27 M^2)] &&
  Resolve[ForAll[{r, M}, M > 0 && r > 2 M, f/r^2 <= 1/(27 M^2)], Reals]];
check["A r is the areal radius (h_{psi psi} = (A r)^2); for A_eta/A = eps constant the condition on A r = c is r^5 [f^2 D - eps^2 Ehat^2 r^2] = 0, a polynomial of degree 7 in r with leading coefficient -eps^2 Ehat^2 and r^5 coefficient Ehat^2: isolated instants only",
  zr[(hR[[2, 2]] /. a -> 0) - (AA r)^2] &&
  Module[{P = Expand[r^5 (f^2 (Dv /. AA -> c/r) - eps^2 EE^2 r^2)]},
    PolynomialQ[P, r] && Exponent[P, r] == 7 && zr[Coefficient[P, r, 7] + eps^2 EE^2] && zr[Coefficient[P, r, 5] - EE^2]]];

Print["--- (D) common launch: orthonormal frame at the launch"];
etaM = DiagonalMatrix[{-1, 1, 1}];
u0 = {gam, sig Cos[chi[lam]], sig Sin[chi[lam]]};
n0 = {sig, gam Cos[chi[lam]], gam Sin[chi[lam]]};
ul = D[u0, lam]; ull = D[u0, {lam, 2}];
rep = sig -> Sqrt[gam^2 - 1];
check["u0 unit timelike; n0 unit, orthogonal to u0 and u0_lambda",
  zr[(u0.etaM.u0 + 1) /. rep] && zr[(n0.etaM.n0 - 1) /. rep] && zr[n0.etaM.u0] && zr[n0.etaM.ul]];
check["g(u0_lambda, u0_lambda) = sigma^2 chi_lambda^2, g(n0, u0_lambdalambda) = -gamma sigma chi_lambda^2",
  zr[ul.etaM.ul - sig^2 chi'[lam]^2] && zr[n0.etaM.ull + gam sig chi'[lam]^2]];

Print["--- (E) tangential component of the magnetic variation field"];
(* frame (T,N): nabla_T T = -B N, N = J T, nabla_T N = B T; V = eta T + xi N; vectors as {T, N} components *)
cd[{x_, y_}] := {D[x, s] + B[s] y, D[y, s] - B[s] x};
Jr[{x_, y_}] := {-y, x};
V = {et[s], xi[s]}; dV = cd[V];
jac = cd[dV] + {0, K[s] xi[s]} + {0, et[s] B'[s] + xi[s] BN[s]} + B[s] Jr[dV];
unit = {et'[s] -> -B[s] xi[s], et''[s] -> -D[B[s] xi[s], s]};
check["unit alpha-speed: g(nabla_T V, T) = eta' + B xi, so eta' = -B xi (eta not zero in general)",
  zr[dV[[1]] - (et'[s] + B[s] xi[s])]];
check["after projection: xi'' + (K + B^2 + dB(N)) xi = 0, the -B' eta terms cancel; tangential part consistent",
  zr[(jac /. unit)[[2]] - (xi''[s] + (K[s] + B[s]^2 + BN[s]) xi[s])] && zr[(jac /. unit)[[1]]]];

Print[If[$Failures === 0, "ALL CHECKS PASSED", ToString[$Failures] <> " CHECK(S) FAILED"]];
If[$Failures =!= 0, Exit[1]];
