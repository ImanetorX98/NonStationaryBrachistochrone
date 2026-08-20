#!/usr/bin/env wolframscript
ok[l_, v_] := Print[If[PossibleZeroQ[Simplify[v]], "  [OK]   ", "  [FAIL] "], l,
   "  ->  ", Simplify[v]];
n = 2;
(* threading form of a GENERAL metric adapted to W = d_chi ; f, om, h may all
   depend on chi -- nothing below uses chi-independence.                      *)
om = Array[ww, n];
hab = Table[hh[Min[i, j], Max[i, j]], {i, n}, {j, n}];
gca = ff om;                                  (* g_{chi a} *)
gab = hab - ff Outer[Times, om, om];
beta = gca/ff;
aab = Simplify[Eh^2 (Outer[Times, gca, gca] + ff gab)/(ff^2 (Eh^2 - ff))];

Print["=== the optical data in threading form ==="];
ok["beta_a - omega_a", Total[Flatten[Simplify[beta - om]]]];
ok["a_ab - Ehat^2 h_ab/(f(Ehat^2-f))",
   Total[Flatten[Simplify[aab - Eh^2 hab/(ff (Eh^2 - ff))]]]];
Print["  |W|^2 = -g_chichi = f, so (H2) reads Ehat^2 > f."];

Print["\n=== (H2) is exactly positive-definiteness of a_ab ==="];
Print["  a_ab = [Ehat^2/(f(Ehat^2-f))] h_ab  with h positive definite,"];
Print["  so a_ab > 0  <=>  Ehat^2 - f > 0  <=>  Ehat > |W|.   (H2)"];
Print["  conformal factor = ", Simplify[Eh^2/(ff (Eh^2 - ff))],
      "  -> +infinity as Ehat -> |W| (freezing)"];

Print["\n=== Randers condition ||beta||_a < 1 ==="];
ainv = Simplify[Inverse[aab]];
nrm = Simplify[beta . ainv . beta];
Print["  a^{ab} beta_a beta_b = ", Simplify[nrm]];
target = ff (Eh^2 - ff) (om . Inverse[hab] . om)/Eh^2;
ok["equals f(Ehat^2-f)|omega|^2_h / Ehat^2", Simplify[nrm - target]];
Print["\n  monotonicity in Ehat:  d/dEhat [ f(Ehat^2-f)/Ehat^2 ] = ",
  Simplify[D[ff (Eh^2 - ff)/Eh^2, Eh]], "  > 0 for f>0, Ehat>0"];
Print["  limit Ehat -> infinity : ", Simplify[Limit[ff (Eh^2 - ff)/Eh^2, Eh -> Infinity]],
      "  = f  (the Fermat value)"];
Print[];
Print["  => ||beta||^2_a = [f(Ehat^2-f)/Ehat^2] |omega|^2_h increases with Ehat"];
Print["     and tends to f |omega|^2_h, the null/Fermat-Randers norm."];
Print["     Hence if the Fermat-Randers condition f|omega|^2_h < 1 holds"];
Print["     (the stationary domain outside the ergosphere), then the rail"];
Print["     Randers condition holds for EVERY Ehat: the rail Finsler"];
Print["     structure inherits, and never enlarges, the null validity domain."];
