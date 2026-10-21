(* Same-launch depth ordering, Thakurta-Kerr equatorial, prograde sector (Paper II, sec:inversion).
   Protocol P1: same launch event, same numerical axial costate J for the t- and tau-branches,
   each radial momentum on its own shell.  M = 1 (a/M <= 1 covers the general case by scaling).
   Checks:
   (1) turning conditions at p_r = 0 give J_t(r) = A Pb/Q  and  J_tau(r) = (A^2/E)(b + f Pb/Q);
   (2) J_t - J_tau = -(A/(E Q)) cond_t,  cond_t = A b Q - Pb (E - A f);
   (3) cond_t < 0 on the whole domain r > 2M, 0 < a <= M, E > 1, A > 0, A^2 f < E^2 (exact, Reduce/CAD);
   (4) contrast: the conformal-normalised condition A^2 b Q - Pb (E - A^2 f) does take positive values;
   (5) the analytic proof of the sign used in the paper: Delta = f P + b^2, P - Delta > 0, the quadratic
       equation of q = Q/Pb, and (J_tau/J_t)^2 = Delta/P - P (q - b/P)^2 < 1, with J_tau/J_t > 0. *)
m = 1; f = 1 - 2 m/r; Dl = r^2 - 2 m r + a^2; P = r^2 + a^2 + 2 m a^2/r; b = 2 m a/r;
vb2 = 1 - A^2 f/e^2; Pb = P + A^2 b^2/e^2; Q = b vb2 + Sqrt[vb2 Dl];
php0 = b vb2/Pb; S = Sqrt[vb2 Dl]/Pb;
Ht0 = j php0 + S j - A;                                  (* J > 0 *)
Htau0 = (j - A^2 b/e) php0 + S (j - A^2 b/e) - A^2 f/e;  (* J - A^2 b/E > 0 *)
Jt = A Pb/Q; Jtau = (A^2/e) (b + f Pb/Q);
fail = 0; ok[name_, c_] := (Print[If[TrueQ[c], "[PASS] ", fail++; "[FAIL] "], name]);
ok["(1) H_t(r,0,J_t(r)) = 0", Simplify[(Ht0 /. j -> Jt) == 0]];
ok["(1) H_tau(r,0,J_tau(r)) = 0", Simplify[(Htau0 /. j -> Jtau) == 0]];
condt = A b Q - Pb (e - A f);
ok["(2) J_t - J_tau = -(A/(E Q)) cond_t", Simplify[Jt - Jtau + (A/(e Q)) condt == 0]];
dom = r > 2 && 0 < a <= 1 && e > 1 && A > 0 && A^2 f < e^2;
ok["(3) Reduce[domain && cond_t >= 0] is False", Reduce[dom && condt >= 0, {r, a, e, A}, Reals] === False];
conde = A^2 b Q - Pb (e - A^2 f);
inst = FindInstance[dom && conde > 0, {r, a, e, A}, Reals];
ok["(4) conformal-normalised condition takes positive values (instance checked)",
   MatchQ[inst, {{__Rule}}] && TrueQ[(dom && conde > 0) /. First[inst]]];
(* (5) analytic proof of the sign, independent of Reduce (x = A/E, u = 1 - x^2 f, q = Q/Pb) *)
asm = {r > 2, 0 < a <= 1, x > 0, 1 - x^2 f > 0};
uu = 1 - x^2 f; Pbx = P + x^2 b^2; qx = (b uu + Sqrt[uu Dl])/Pbx;
ok["(5a) Delta = f P + b^2", Simplify[Dl - (f P + b^2)] === 0];
ok["(5b) P - Delta = 2M(r^2+a^2)/r > 0", Simplify[P - Dl - 2 m (r^2 + a^2)/r] === 0];
ok["(5c) Pb q^2 - 2 b u q - f u = 0", FullSimplify[Pbx qx^2 - 2 b uu qx - f uu, asm] === 0];
ok["(5d) (J_tau/J_t)^2 = x^2 (f + b q)^2 = Delta/P - P (q - b/P)^2",
   FullSimplify[x^2 (f + b qx)^2 - (Dl/P - P (qx - b/P)^2), asm] === 0 &&
   Simplify[(Jtau/Jt) - x (f + b Q/Pb) /. A -> x e] === 0];
ok["(5e) J_tau = (A^2/E) Sqrt[Delta/u], increasing in r",
   FullSimplify[(Jtau - (A^2/e) Sqrt[Dl/vb2]) /. A -> x e, Join[asm, {e > 1}]] === 0];
Print[If[fail == 0, "ALL CHECKS PASSED", ToString[fail] <> " FAILED"]];
If[fail > 0, Exit[1]];
