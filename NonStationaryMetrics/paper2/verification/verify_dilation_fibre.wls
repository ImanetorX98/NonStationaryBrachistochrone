#!/usr/bin/env wolframscript
(* verify_dilation_fibre.wls
   Meena-Zawadzki (Bull. Malays. Math. Sci. Soc. 47:101, 2024) prove that a
   Clairaut conformal submersion has CONSTANT DILATION ALONG ITS FIBRES, the
   fibres being totally umbilical with mean curvature the gradient of a
   function.  We test our optical submersion against that necessary condition. *)
$Failures = 0;
check[l_, c_] := (Print[If[TrueQ[c], "PASS  ", "FAIL  "], l];
   If[! TrueQ[c], $Failures++]);
ClearAll[eta, r, M, Eh, A];
f[r_] := 1 - 2 M/r;
(* Perlick dilation for the conformal metric: A^2 f replaces f *)
Lam2[r_, eta_] := Eh^2/(A[eta]^2 f[r] (Eh^2 - A[eta]^2 f[r]));
Print["=========================================================="];
Print[" is the optical submersion Clairaut?  (Meena-Zawadzki test)"];
Print["=========================================================="];
Print["      Lambda^2 = ", Lam2[r, eta]];
(* the fibres are the orbits of W = d/d(eta): constancy along a fibre means
   d(Lambda^2)/d(eta) = 0 *)
dfib = Simplify[D[Lam2[r, eta], eta]];
Print["      d(Lambda^2)/d(eta) = ", Simplify[dfib]];
check["the dilation is NOT constant along the fibres when A runs",
      Simplify[dfib] =!= 0];
check["it IS constant along the fibres at frozen A",
      Simplify[dfib /. A -> Function[x, 1]] === 0];
(* the obstruction is proportional to A' , i.e. to psi *)
psi = A'[eta]/A[eta];
ratio = Simplify[dfib/psi];
check["the obstruction is proportional to A'/A = psi",
      FreeQ[Simplify[ratio], Derivative[1][A]]];
Print["      d(Lambda^2)/d(eta) / psi = ", Simplify[ratio]];
Print["      => by the Meena-Zawadzki necessary condition the optical"];
Print["         submersion is Clairaut exactly at frozen A, and the failure"];
Print["         is again measured by psi."];
Print["=========================================================="];
Print[If[$Failures === 0, "ALL CHECKS PASSED",
     ToString[$Failures] <> " CHECK(S) FAILED"]];
If[$Failures > 0, Exit[1], Exit[0]];
