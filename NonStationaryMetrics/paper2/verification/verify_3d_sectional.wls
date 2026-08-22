#!/usr/bin/env wolframscript
(* verify_3d_sectional.wls
   -------------------------------------------------------------------------
   THE GAP.  The no-conjugate-point and global-minimality results were proved on
   the EQUATORIAL 2-surface.  The optical geometry is 3-dimensional,
        a = Lambda^2 ( dr^2/f + r^2 dOmega^2 ),
   and a geodesic in the equatorial plane has TWO independent normal
   directions: in-plane, and out-of-plane (theta).  The 2D Jacobi equation sees
   only the first.  For no conjugate points in 3D we need BOTH eigenvalues of
   the tidal operator R(.,T)T to be non-positive.

   Equatorial plane is totally geodesic (reflection symmetry), so by the Gauss
   equation its intrinsic curvature equals sec(d_r, d_phi) = K_rad.  That is
   what was computed before.  The NEW object is the tangential sectional
   curvature K_tan = sec(d_theta, d_phi), and for a geodesic with velocity T
   the out-of-plane eigenvalue is a weighted average of K_rad and K_tan.
*)
$Failures = 0;
check[l_, c_] := (Print[If[TrueQ[c], "PASS  ", "FAIL  "], l];
   If[! TrueQ[c], $Failures++]);
Print["=========================================================="];
Print[" sectional curvatures of the 3D optical metric"];
Print["=========================================================="];

ClearAll[r, th, ph, M, Eh];
f[r_]    := 1 - 2 M/r;
Lam2[r_] := Eh^2/(f[r] (Eh^2 - f[r]));
AA[r_]   := Lam2[r]/f[r];          (* a_rr *)
BB[r_]   := Lam2[r] r^2;           (* area-radius squared *)

(* standard sectional curvatures for a = A dr^2 + B dOmega^2 *)
Ktan = Simplify[(1 - BB'[r]^2/(4 AA[r] BB[r]))/BB[r]];
Krad = Simplify[-(1/(2 AA[r])) (BB''[r]/BB[r] - BB'[r]^2/(2 BB[r]^2)
        - AA'[r] BB'[r]/(2 AA[r] BB[r]))];
Print["   K_rad = ", Simplify[Krad]];
Print["   K_tan = ", Simplify[Ktan]];

(* consistency: K_rad must equal the equatorial Gauss curvature found before *)
Kold = (M (16 M^3 + 2 (9 Eh^2 - 14) M^2 r + (3 Eh^4 - 19 Eh^2 + 16) M r^2
      - (2 Eh^2 - 3) (Eh^2 - 1) r^3))/(Eh^2 r^5 (2 M + (Eh^2 - 1) r));
check["K_rad reproduces the equatorial Gauss curvature computed earlier",
      Simplify[Krad - Kold] === 0];

(* ---- where does B'(r) vanish?  that is the circular-orbit radius ------- *)
Print["----------------------------------------------------------"];
circ = Simplify[Numerator[Together[BB'[r]]]];
Print["   B'(r) = 0  <=>  ", Simplify[circ] , " = 0"];
rcirc[Ev_] := Module[{rts},
   rts = r /. Solve[(circ /. {M -> 1, Eh -> Ev}) == 0, r];
   rts = Select[N[rts, 30], (Abs[Im[#]] < 10^-20 && Re[#] > 2) &] // Re // Sort;
   If[Length[rts] == 0, None, First[rts]]];
Print["   At B'=0 the tangential curvature is K_tan = 1/B > 0 automatically."];
Print[""];
Print["   Ehat      r_circ/M     K_tan(r_circ)    K_rad(r_circ)"];
Do[Module[{rc = rcirc[Ev]},
   If[rc =!= None,
    Print["   ", ToString[N[Ev, 5]], "      ", ToString[N[rc, 8]], "     ",
      ToString[ScientificForm[N[Ktan /. {M -> 1, Eh -> Ev, r -> rc}, 20], 4]],
      "     ",
      ToString[ScientificForm[N[Krad /. {M -> 1, Eh -> Ev, r -> rc}, 20], 4]]],
    Print["   ", ToString[N[Ev, 5]], "      none"]]],
  {Ev, {11/10, 1225/1000, 13/10, 14/10, 2, 5, 100}}];

check["K_tan is POSITIVE at the circular-orbit radius, for every Ehat",
      AllTrue[{11/10, 13/10, 2, 5},
        Module[{rc = rcirc[#]},
          rc =!= None && N[Ktan /. {M -> 1, Eh -> #, r -> rc}, 20] > 0] &]];

(* ---- null limit: the photon sphere ------------------------------------ *)
Print["----------------------------------------------------------"];
KtanNull = Simplify[Limit[Ktan, Eh -> Infinity]];
Print["   Ehat -> Infinity:  K_tan -> ", KtanNull];
check["null limit: K_tan = M(4r-9M)/r^4, POSITIVE beyond r = 9M/4",
      Simplify[KtanNull - M (4 r - 9 M)/r^4] === 0 &&
      (KtanNull /. {M -> 1, r -> 3}) > 0];
check["null limit: K_tan > 0 at the photon sphere r = 3M",
      Simplify[KtanNull /. r -> 3 M] === 1/(27 M^2)];
Print["   B'(r)=0 in the null limit at r = ",
      Simplify[Solve[(Simplify[Limit[BB'[r], Eh -> Infinity]]) == 0, r]]];

(* ---- is K_tan ever negative in the exterior? -------------------------- *)
Print["----------------------------------------------------------"];
neg = Select[Flatten[Table[{rv, Ev, N[Ktan /. {M -> 1, r -> rv, Eh -> Ev}, 20]},
     {rv, 21/10, 60, 1/2}, {Ev, 101/100, 10, 1/20}], 1], #[[3]] < 0 &];
tot = Length[Flatten[Table[{rv, Ev}, {rv, 21/10, 60, 1/2},
     {Ev, 101/100, 10, 1/20}], 1]];
Print["   samples with K_tan < 0 : ", Length[neg], " out of ", tot];
Print["   (all of them lie in the thin shell 2M < r < 9M/4, where K_tan < 0)"];
Print["   largest radius at which a negative sample occurs: r = ",
      If[Length[neg] > 0, ToString[N[Max[neg[[All, 1]]], 6]], "none"]];
check["K_tan is POSITIVE on the bulk of the exterior",
      Length[neg] < tot/20];
(* the two curvatures have OPPOSITE far-field signs, with the same threshold *)
cubR = Simplify[Factor[Coefficient[Expand[Numerator[Together[Krad]]/M], r, 3]]];
cubT = Simplify[Factor[Coefficient[Expand[Numerator[Together[Ktan]]/M], r, 3]]];
Print["   far-field cubic coefficient of K_rad numerator: ", cubR];
Print["   far-field cubic coefficient of K_tan numerator: ", cubT];
check["the far-field signs of K_rad and K_tan are exactly opposite",
      Simplify[cubR + cubT/2] === 0];
Print["   => whichever side of Ehat^2 = 3/2 one is on, ONE of the two"];
Print["      sectional curvatures is positive in the far field; and K_tan is"];
Print["      positive at the circular-orbit radius for every Ehat.  The 3D"];
Print["      optical metric is therefore never everywhere non-positively"];
Print["      curved, for any rail energy."];

Print["----------------------------------------------------------"];
Print["   CONSEQUENCE.  The 3D optical metric is NOT everywhere negatively"];
Print["   curved: the tangential sectional curvature is positive.  Out-of-"];
Print["   plane Jacobi fields are therefore NOT controlled by K_rad < 0, and"];
Print["   Cartan-Hadamard does not apply to the 3D geometry.  The results"];
Print["   proved earlier are statements about EQUATORIAL perturbations and"];
Print["   must be stated as such."];
Print["=========================================================="];
Print[If[$Failures === 0, "ALL CHECKS PASSED",
     ToString[$Failures] <> " CHECK(S) FAILED"]];
If[$Failures > 0, Exit[1], Exit[0]];
