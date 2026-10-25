(* ::Package:: *)
(* =====================================================================
   SYMBOLIC FAY / beta in closed form - Vaidya separatrix (tau branch)
   ---------------------------------------------------------------------
   Goal: close in SYMBOLIC FORM the coefficients beta of the 5 independent
   elliptic dilogarithms of  delta phi|_sep, which in Python stay numerical
   (~1e-10) because of ill-conditioning. Mathematica gets round the wall with:
     (A) arbitrary precision (>=60 digits) -> clean rank, beta to ~50 digits;
     (B) symbolic Weierstrass addition theorems -> EXACT Fay relations.
   Conventions IDENTICAL to the Python script vaidya_sep_5term.py:
     curve E: w^2 = Q4(r),  sqrt(S)=(r-r_d)sqrt(Q4),  z = Int dr/sqrt(Q4),
     base point z0 = z(r0=12),  period lattice (2 om1, 2 I w_im).
   ===================================================================== *)

$MinPrecision = 80; $MaxPrecision = 80; prec = 70;
Off[NIntegrate::ncvb]; Off[NIntegrate::slwcon]; Off[NIntegrate::nlim];
Off[Power::infy]; Off[Infinity::indet]; Off[General::stop]; Off[NIntegrate::inumr];

(* ---------- 0. physical data ---------- *)
Em = 7/5; mm = 1; r0 = 12;                       (* E, M, base point *)
DEpol[r_] := (Em^2 - 1) r + 2 mm;
Spoly[r_, J_] := r (r - 2 mm) DEpol[r] (r^2 (r - 2 mm) - J^2 DEpol[r]);

(* Jc = separatrix value: double root of S  <=>  Resultant(S,S',r)=0 *)
JcExact = Jc /. Solve[Resultant[Spoly[r, Jc], D[Spoly[r, Jc], r], r] == 0 &&
                       Jc > 1, Jc, Reals][[1]];
Print["Jc (exact) = ", JcExact, "  ~ ", N[JcExact, 20]];

Jcn = N[JcExact, prec];                          (* high-precision NUMERICAL Jc *)
Sc = Expand[Spoly[r, Jcn]];                      (* polynomial with numerical coefficients *)
a4 = N[Coefficient[Sc, r, 6], prec];             (* leading coefficient of S = leading coefficient of Q4 (degree 6 -> a4 r^6) *)
rootsS = Sort[N[r /. Solve[Sc == 0, r], prec], Re[#1] < Re[#2] &];
(* double root = closest pair *)
rd = First@SelectFirst[
   Table[{rootsS[[i]], rootsS[[j]]}, {i, Length[rootsS]}, {j, i + 1, Length[rootsS]}] // Flatten[#, 1] &,
   Abs[#[[1]] - #[[2]]] < 10^-6 &];
rd = Re[rd]; Print["r_d = ", rd];
(* Q4 = product of the 4 NON-double factors (roots != rd) * a4 *)
others = DeleteCases[rootsS, x_ /; Abs[x - rd] < 10^-6, 1, 2];  (* remove the two copies of rd *)
erts = Sort[Re /@ others];
{e1, e2, e3, e4} = erts;
Q4[x_] := a4 (x - e1)(x - e2)(x - e3)(x - e4);
Print["roots of Q4 = ", erts, "   a4 = ", a4];

(* ---------- 1. lattice and Weierstrass invariants ---------- *)
(* half-periods from the Legendre moduli (as in the Python script) *)
k2 = ((e3 - e2)(e4 - e1))/((e4 - e2)(e3 - e1));
pref = 2/Sqrt[(e4 - e2)(e3 - e1)]/Sqrt[a4];
om1 = N[pref EllipticK[k2], prec];
wim = N[pref EllipticK[1 - k2], prec];
Print["om1 = ", om1, "   w_im = ", wim, "   tau = ", N[I wim/om1, 12]];
(* invariants g2, g3 from the half-periods (om1, I w_im) *)
{g2, g3} = WeierstrassInvariants[{om1, I wim}];
WP[z_] := WeierstrassP[z, {g2, g3}];
WZ[z_] := WeierstrassZeta[z, {g2, g3}];
WS[z_] := WeierstrassSigma[z, {g2, g3}];
WPp[z_] := WeierstrassPPrime[z, {g2, g3}];

(* z(r) = Int_{e4}^r dr/sqrt(Q4)  (Abel coordinate) *)
zr[rv_] := NIntegrate[1/Sqrt[Q4[x]], {x, e4, rv}, WorkingPrecision -> prec];
z0 = zr[r0];
zinf = NIntegrate[1/Sqrt[Q4[x]], {x, e4, Infinity}, WorkingPrecision -> prec];
zd = zinf + NIntegrate[1/Sqrt[Q4[x]], {x, -Infinity, rd}, WorkingPrecision -> prec];
iw = I wim;
Print["z_inf = ", zinf, "   z_d = ", zd];

(* check: explicit r(z)  = c_r - (1/Sqrt a4)[zeta(z-zinf)-zeta(z+zinf)] *)
sa = Sqrt[a4]; cr = e4 - (2/sa) WZ[zinf];
rOfz[z_] := cr - (1/sa)(WZ[z - zinf] - WZ[z + zinf]);
Print["check r(z): r(z(10)) - 10 = ", N[rOfz[zr[10]] - 10, 12]];

(* ---------- 2. residues (coefficients of R and eta') in CLOSED form ---------- *)
(* N_m = S dm K - 1/2 K dm S  with M symbolic, then M->1 *)
DEm[r_, M_] := (Em^2 - 1) r + 2 M;
Sm[r_, M_, J_] := r (r - 2 M) DEm[r, M] (r^2 (r - 2 M) - J^2 DEm[r, M]);
Km[r_, M_, J_] := J DEm[r, M];
NmSym = Expand[ (Sm[r, M, Jcn] D[Km[r, M, Jcn], M]
                 - (1/2) Km[r, M, Jcn] D[Sm[r, M, Jcn], M]) /. M -> mm ];
Nmf[x_] := NmSym /. r -> x;
Qd[x_] := (D[Q4[y], y] /. y -> x);  Qdd[x_] := (D[Q4[y], {y, 2}] /. y -> x);
Frat[x_] := Nmf[x]/Q4[x];
(* triple pole z_d :  s=sqrt(Q4(rd)), a1,a2, h0,h1,h2  (Laurent identity) *)
sP = Sqrt[Q4[rd]]; a1c = Qd[rd]/(4 sP); a2c = Qdd[rd]/12;
h0 = Frat[rd];
h1 = (D[Frat[x], x] /. x -> rd) sP;
h2 = (1/2)((D[Frat[x], {x, 2}] /. x -> rd) sP^2 + (D[Frat[x], x] /. x -> rd)(Qd[rd]/2));
b1zd = (h2 - 3 a1c h1 + (6 a1c^2 - 3 a2c) h0)/sP^3;   (* Res_{z_d}(R) *)
b2zd = (h1 - 3 a1c h0)/sP^3;
b3zd = h0/sP^3;
b2h[ei_] := Nmf[ei]/(ei - rd)^3 (4/Qd[ei]^2);          (* order-2 residue at the half-periods *)
e1zd = (rd^3 - 2 rd^2)/sP;                              (* Res_{z_d}(eta') *)
Bc = cr + (1/sa) WZ[2 zinf]; Aq = -1/sa;
e1zi = Aq (2 Bc + rd - 2); e2zi = Aq^2;                 (* Res_{z_inf}(eta') *)
Print["residues: b1zd=", N[b1zd,10], " b2zd=", N[b2zd,10], " b3zd=", N[b3zd,10]];
Print["          e1zd=", N[e1zd,10], " e1zi=", N[e1zi,10], " e2zi=1/a4=", N[e2zi,10]];

(* which e_i maps to the half-period I w_im: the one with N_m(e_i)!=0 (e4->0) *)
eIw = SelectFirst[{e1, e2, e3}, Abs[N[Nmf[#]]] > 10^-3 Abs[N[Nmf[e4]]] &];
Print["e -> z=I w_im : ", eIw, "   e -> z=0 : ", e4];

(* ---------- 3. the 8 EVEN functions e_i and their primitives Pe_i ---------- *)
(* e = {1, Z_zd, P_zd, Pp_zd, wp0, wpiw, P_zi, Z_zi} ;  Z_a=z(-a)-z(+a) etc. *)
ef[1][z_] := 1;
ef[2][z_] := WZ[z - zd] - WZ[z + zd];
ef[3][z_] := WP[z - zd] + WP[z + zd];
ef[4][z_] := WPp[z - zd] - WPp[z + zd];
ef[5][z_] := WP[z];
ef[6][z_] := WP[z - iw];
ef[7][z_] := WP[z - zinf] + WP[z + zinf];
ef[8][z_] := WZ[z - zinf] - WZ[z + zinf];
(* primitives Pe_i = Int_{z0}^z e_i  (CLOSED form: log-sigma / zeta / z) *)
Pe[1][z_] := z - z0;
Pe[2][z_] := (Log[WS[z - zd]] - Log[WS[z + zd]]) - (Log[WS[z0 - zd]] - Log[WS[z0 + zd]]);
Pe[3][z_] := (-(WZ[z - zd] + WZ[z + zd])) - (-(WZ[z0 - zd] + WZ[z0 + zd]));
Pe[4][z_] := (WP[z - zd] - WP[z + zd]) - (WP[z0 - zd] - WP[z0 + zd]);
Pe[5][z_] := (-WZ[z]) - (-WZ[z0]);
Pe[6][z_] := (-WZ[z - iw]) - (-WZ[z0 - iw]);
Pe[7][z_] := (-(WZ[z - zinf] + WZ[z + zinf])) - (-(WZ[z0 - zinf] + WZ[z0 + zinf]));
Pe[8][z_] := (Log[WS[z - zinf]] - Log[WS[z + zinf]]) - (Log[WS[z0 - zinf]] - Log[WS[z0 + zinf]]);

(* coefficients of R (c) and eta' (d) on the 8 functions ; C0, Ce fixed by matching *)
cv = {0, b1zd, b2zd, -b3zd/2, b2h[e4], b2h[eIw], 0, 0};
dv = {0, e1zd, 0, 0, 0, 0, e2zi, e1zi};
zt = 191/1000;
Rexact = Nmf[rOfz[zt]]/((rOfz[zt] - rd)^3 Q4[rOfz[zt]]);
Epexact = (rOfz[zt]^3 - 2 rOfz[zt]^2)/(rOfz[zt] - rd);
cv[[1]] = Rexact - Sum[cv[[i]] ef[i][zt], {i, 2, 8}];
dv[[1]] = Epexact - Sum[dv[[i]] ef[i][zt], {i, 2, 8}];
Print["C0 = ", N[cv[[1]], 12], "   Ce = ", N[dv[[1]], 12]];

(* ---------- 4. antisymmetric dilogarithms and Fay reduction ---------- *)
(* A[i,j] = Int_{z0}^z (e_i Pe_j - e_j Pe_i) dz  (antisymmetric = genuine dilogarithm) *)
Rall = {1, 2, 3, 4, 5, 6}; Eall = {1, 2, 7, 8};   (* non-zero (1-based) indices in R, eta' *)
pairs = Select[Subsets[Range[8], {2}],
               (MemberQ[Rall, #[[1]]] && MemberQ[Eall, #[[2]]]) ||
               (MemberQ[Eall, #[[1]]] && MemberQ[Rall, #[[2]]]) &];
Print["A pairs = ", Length[pairs]];
wk = Table[(cv[[p[[1]]]] dv[[p[[2]]]] - cv[[p[[2]]]] dv[[p[[1]]]])/2, {p, pairs}]; (* coefficients of A_k *)

(* --- 4a. EXACT REDUCTION through Weierstrass addition (the symbolic core) ---
   Idea: the integrand A'_k = e_i Pe_j - e_j Pe_i has weight 1. The Fay relations
   are the combinations Sum lambda_k A'_k that are total derivatives of weight-1 functions.
   Represent every e_i(z), Pe_i(z) as a function of {P=WP[z], Pp=WPp[z], Z=WZ[z], z}
   through the addition theorems, with the values at the poles as numerical CONSTANTS
   (WP[zd], WPp[zd], WZ[zd], ...). The relations then become rational identities. *)

(* addition theorems (u+-v):  constants at the poles *)
PZval[a_] := {WP[a], WPp[a], WZ[a]};        (* numerical {℘(a),℘'(a),ζ(a)} *)
addP[Pv_, Ppv_, pa_, ppa_] :=              (* ℘(z - a) from ℘(z),℘'(z) and the constants at a *)
   (1/4)((Ppv + ppa)/(Pv - pa))^2 - Pv - pa;
addZ[Pv_, Ppv_, Zv_, pa_, ppa_, za_] :=     (* ζ(z - a) *)
   Zv - za + (1/2)(Ppv + ppa)/(Pv - pa);

Print["--- setup complete. Ready for the symbolic reduction / null space. ---"];
Print["NB: the exact algebraic reduction of part (4a) is still to be developed here: ",
      "rewrite e_i, Pe_i in {P,Pp,Z,z}, form Sum lambda_k A'_k, impose that it is a ",
      "total derivative, exact Solve[] -> Fay relations -> symbolic beta."];

(* ---------- 5. high-precision NUMERICAL extraction of beta (robust fallback) ---------- *)
(* A'_k(z) = e_i Pe_j - e_j Pe_i  (closed). Weight-1 basis W'. Solved on many points at 70 digits. *)
dApair[i_, j_, z_] := ef[i][z] Pe[j][z] - ef[j][z] Pe[i][z];
(* weight-1 basis (derivatives): d/dz(Pe_i Pe_j), d/dz(z Pe_i), d/dz z^2, d/dz z *)
Wb = Join[
   Flatten[Table[With[{ii = i, jj = j}, Function[z, ef[ii][z] Pe[jj][z] + ef[jj][z] Pe[ii][z]]], {i, 8}, {j, i, 8}], 1],
   Table[With[{ii = i}, Function[z, Pe[ii][z] + z ef[ii][z]]], {i, 8}],
   {Function[z, 2 z], Function[z, 1]}];
npts = 40;
zsamp = Table[zr[rr], {rr, Subdivide[N[119/10, prec], N[89/10, prec], npts - 1]}];
LHS = Table[N[Sum[wk[[k]] dApair[pairs[[k, 1]], pairs[[k, 2]], zz], {k, Length[pairs]}], prec], {zz, zsamp}];
(* choose 5 basis dilogarithms (as in Python): {wp0,P_zi},{wp0,Z_zi},{Z_zd,wp0},{1,wp0},{Pp_zd,P_zi} *)
baseNames = {{5, 7}, {5, 8}, {2, 5}, {1, 5}, {4, 7}};
baseIdx = Flatten[Position[pairs, Sort[#]] & /@ baseNames];
Amat = Table[N[Join[
    Table[dApair[pairs[[b, 1]], pairs[[b, 2]], zz], {b, baseIdx}],
    Table[w[zz], {w, Wb}]], prec], {zz, zsamp}];
sol = LeastSquares[Amat, LHS];
resid = Max[Abs[Amat.sol - LHS]];
beta = sol[[1 ;; 5]];
Print["residual (should be ~10^-60 if the rank is exactly 5): ", N[resid, 5]];
Print["beta (high precision):"]; Do[Print["  beta[", b, "] = ", N[beta[[b]], 40]], {b, 5}];

(* ---------- 6. closed-form recognition of the beta ---------- *)
(* try to recognise each beta as algebraic or as a combination of the residues/periods *)
consts = {b1zd, b2zd, b3zd, b2h[e4], b2h[eIw], e1zd, e1zi, e2zi, om1, wim, 1};
Print["--- recognition attempt (RootApproximant / FindIntegerNullVector) ---"];
Do[
  Module[{v = N[Append[consts, beta[[b]]], 55], ra, nv},
    ra = Quiet@RootApproximant[N[beta[[b]], 45]];
    nv = Quiet@Check[FindIntegerNullVector[v], $Failed];
    Print["  beta[", b, "]: RootApproximant -> ", ra];
    Print["            FindIntegerNullVector(consts,beta) -> ", nv]],
  {b, 5}];
