#!/usr/bin/env wolframscript
(* verify_which_clock.wls
   WHICH clock does the optical metric measure?  The global-minimality argument
   claims that a-length IS arrival time, so this had better be checked.

   Static case: g = -f dt^2 + h.  Rail: -g(u,d/dt) = f u^t = Ehat.
   Normalisation: -f (u^t)^2 + |v|_h^2 = -1  =>  |v|_h^2 = (Ehat^2 - f)/f.     *)
$Failures = 0;
check[l_, c_] := (Print[If[TrueQ[c], "PASS  ", "FAIL  "], l];
   If[! TrueQ[c], $Failures++]);
ClearAll[f, Eh];
Print["=========================================================="];
Print[" which clock is the optical length?"];
Print["=========================================================="];
ut = Eh/f;                                   (* dt/dtau *)
vh = Sqrt[(Eh^2 - f)/f];                     (* |dx/dtau|_h *)
LamT = Simplify[ut/vh];                      (* dt per unit h-length *)
LamTau = Simplify[1/vh];                     (* dtau per unit h-length *)
Print["      dt   per unit h-length : ", Simplify[LamT]];
Print["      dtau per unit h-length : ", Simplify[LamTau]];
check["coordinate-time factor is Perlick's: Lambda_t^2 = Ehat^2/(f(Ehat^2-f))",
      Simplify[LamT^2 - Eh^2/(f (Eh^2 - f))] === 0];
check["proper-time factor is DIFFERENT: Lambda_tau^2 = f/(Ehat^2-f)",
      Simplify[LamTau^2 - f/(Eh^2 - f)] === 0];
check["they differ by exactly Ehat^2/f^2",
      Simplify[LamT^2/LamTau^2 - Eh^2/f^2] === 0];
Print["      => Perlick's metric measures COORDINATE time, not proper time."];

(* ---- Gauss curvature of BOTH optical surfaces, Schwarzschild equatorial --- *)
Print["----------------------------------------------------------"];
ClearAll[r, M, Eh];
ff[r_] := 1 - 2 M/r;
gauss[Lam2_] := Module[{aa, ai, Ch, Rm, Rc, co},
   co = {r, ph};
   aa = Lam2 {{1/ff[r], 0}, {0, r^2}};
   ai = Simplify[Inverse[aa]];
   Ch = Table[Simplify[(1/2) Sum[ai[[i, l]] (D[aa[[l, j]], co[[k]]] +
        D[aa[[l, k]], co[[j]]] - D[aa[[j, k]], co[[l]]]), {l, 2}]],
      {i, 2}, {j, 2}, {k, 2}];
   Rm = Table[Simplify[D[Ch[[i, j, l]], co[[k]]] - D[Ch[[i, j, k]], co[[l]]] +
       Sum[Ch[[i, k, s]] Ch[[s, j, l]] - Ch[[i, l, s]] Ch[[s, j, k]], {s, 2}]],
      {i, 2}, {j, 2}, {k, 2}, {l, 2}];
   Rc = Table[Simplify[Sum[Rm[[s, i, s, j]], {s, 2}]], {i, 2}, {j, 2}];
   Simplify[(1/2) Sum[ai[[i, j]] Rc[[i, j]], {i, 2}, {j, 2}]]];

Kt   = gauss[Eh^2/(ff[r] (Eh^2 - ff[r]))];
Ktau = gauss[ff[r]/(Eh^2 - ff[r])];
Print["      K_t   (coordinate time) = ", Kt];
Print["      K_tau (proper time)     = ", Ktau];
check["the two curvatures are genuinely different functions",
      Simplify[Kt - Ktau] =!= 0];

(* far-field behaviour of each *)
at = Simplify[Limit[r^3 Kt, r -> Infinity]];
atau = Simplify[Limit[r^3 Ktau, r -> Infinity]];
Print["      r^3 K_t   -> ", at];
Print["      r^3 K_tau -> ", atau];
check["the t-branch threshold is Ehat^2 = 3/2, as previously found",
      Simplify[at + M (2 Eh^2 - 3)/Eh^2] === 0];

(* is K_tau negative throughout the exterior? scan *)
Print["----------------------------------------------------------"];
Print["      sign of K_tau at M = 1:"];
KtN[rv_, ev_] := N[Ktau /. {M -> 1, r -> rv, Eh -> ev}, 20];
Do[Print["        Ehat = ", ToString[N[ev, 4]], "   ",
   StringJoin @@ Table["K(" <> ToString[rv] <> ")=" <>
     ToString[ScientificForm[KtN[rv, ev], 3]] <> "  ",
     {rv, {2.5, 3, 6, 30, 200}}]],
  {ev, {11/10, 1225/1000, 13/10, 2, 5}}];

Print["=========================================================="];
Print[If[$Failures === 0, "ALL CHECKS PASSED",
     ToString[$Failures] <> " CHECK(S) FAILED"]];
If[$Failures > 0, Exit[1], Exit[0]];
