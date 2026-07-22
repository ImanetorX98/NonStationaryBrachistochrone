(* ===== Independent Mathematica verification of the frozen-Schwarzschild derivations ===== *)
Off[Solve::ratnz];
b_ := b; (* b=E^2-1>0, M=1 *)
Vf[r_]  := r (r^2 - 2 r)/(b r + 2);
Vpf[r_] := D[Vf[r], r];
Kf[r_]  := 1/Sqrt[r (r - 2)];
Wf[r_]  := Kf[r]/Vpf[r];

Print["=== (1) 2 V W  =?=  DE Sqrt[Delta]/N ==="];
lhs = Simplify[2 Vf[r] Wf[r], Assumptions -> {r > 2, b > 0}];
rhs = (b r + 2) Sqrt[r^2 - 2 r]/(b r^2 + (3 - b) r - 4);
Print["  Simplify[lhs-rhs] = ", Simplify[lhs - rhs, Assumptions -> {r > 2, b > 0}]];

Print["=== (2) large-r expansion of 2VW  (expect 1 - (1/b)/r + ...) ==="];
ser = Series[rhs /. r -> 1/e, {e, 0, 2}] // Normal // Simplify;
Print["  2VW (e=1/r) = ", ser];

Print["=== (2b) W(V) large-V:  1/(2V) - 1/(2 b^{3/2} V^{3/2}) + ... ==="];
(* parametrize: as r->inf, V ~ r^2/b. Get W as series in 1/r then substitute r->Sqrt[b V] leading *)
Wser = Series[Wf[r], {r, Infinity, 3}] // Normal;
Print["  W(r->inf) = ", Simplify[Wser]];

Print["=== (3) MODEL Abel integral: sqrt(x) Int_x^V0 (1/(2A))/Sqrt[A-x] dA = arctan Sqrt[V0/x-1] ==="];
Imodel = Integrate[(1/(2 A))/Sqrt[A - x], {A, x, V0}, Assumptions -> {0 < x < V0}];
Print["  sqrt(x)*I = ", Simplify[Sqrt[x] Imodel, Assumptions -> {0 < x < V0}]];
Print["  minus ArcTan[Sqrt[V0/x-1]] = ", Simplify[Sqrt[x] Imodel - ArcTan[Sqrt[V0/x - 1]], Assumptions -> {0 < x < V0}]];

Print["=== (4) CORRECTION Abel: -1/(2 b^{3/2}) sqrt(x) Int A^{-3/2}(A-x)^{-1/2} dA = -b^{-3/2} sqrt(1/x-1/V0) ==="];
Icorr = Integrate[A^(-3/2)/Sqrt[A - x], {A, x, V0}, Assumptions -> {0 < x < V0}];
lhsC = -(1/(2 b^(3/2))) Sqrt[x] Icorr;
rhsC = -b^(-3/2) Sqrt[1/x - 1/V0];
Print["  lhsC - rhsC = ", Simplify[lhsC - rhsC, Assumptions -> {0 < x < V0, b > 0}]];

Print["=== (5) Phi' from Phi=Phi_model+Phi_corr:  expect [1/(2 sqrt(x(V0-x)))](b^{-3/2} sqrt(V0)/x - 1) ==="];
Phi = ArcTan[Sqrt[V0/x - 1]] - b^(-3/2) Sqrt[1/x - 1/V0];
Phip = Simplify[D[Phi, x], Assumptions -> {0 < x < V0, b > 0}];
target = (1/(2 Sqrt[x (V0 - x)])) (b^(-3/2) Sqrt[V0]/x - 1);
Print["  Phi' - target = ", Simplify[Phip - target, Assumptions -> {0 < x < V0, b > 0}]];

Print["=== (6) x_pk (bracket=0) and Phi''(x_pk) leading coeff  (expect -1/2 b^{9/4} V0^{-5/4}) ==="];
xpk = b^(-3/2) Sqrt[V0];
Phipp = D[Phi, {x, 2}];
(* evaluate Phi'' at x=xpk, then leading order as V0->inf *)
val = Phipp /. x -> xpk;
lead = Series[val, {V0, Infinity, 1}] // Normal // Simplify;
Print["  Phi''(x_pk) leading = ", lead];
Print["  compare -1/2 b^{9/4} V0^{-5/4}: ratio = ", Simplify[lead/(-(1/2) b^(9/4) V0^(-5/4))]];

Print["=== (7) kernel integral Int_x^V0 [(2x-t)/sqrt(t-x)+(V0-2x)/sqrt(V0-x)] dt = sqrt(V0-x)(V0+2x)/3 ==="];
kerI = Integrate[(2 x - t)/Sqrt[t - x] + (V0 - 2 x)/Sqrt[V0 - x], {t, x, V0}, Assumptions -> {0 < x < V0}];
Print["  diff = ", Simplify[kerI - Sqrt[V0 - x] (V0 + 2 x)/3, Assumptions -> {0 < x < V0}]];
Print["=== (6-fix) Phi''(x_pk) leading via V0=t^4  (expect -1/2 b^{9/4} V0^{-5/4} = -1/2 b^{9/4} t^{-5}) ==="];
Phi = ArcTan[Sqrt[V0/x - 1]] - b^(-3/2) Sqrt[1/x - 1/V0];
Phipp = D[Phi, {x, 2}];
xpk = b^(-3/2) Sqrt[V0];
val = Phipp /. x -> xpk;
(* substitute V0 -> t^4, expand for t->inf *)
valT = val /. V0 -> t^4;
ser = Series[valT, {t, Infinity, 6}] // Normal // Simplify;
Print["  Phi''(x_pk) in t=V0^{1/4}: ", ser];
Print["  leading coeff * (-2 b^{-9/4} t^5) [expect 1] = ",
   Simplify[Limit[ser*(-2 b^(-9/4) t^5), t -> Infinity], Assumptions->{b>0}]];

Print["=== (2c) W(V) correction coefficient in V: expect W - 1/(2V) ~ -1/(2 b^{3/2}) V^{-3/2} ==="];
(* r as function of V: invert V=r^2(r-2)/(br+2) asymptotically. Use V=r^2/b(1+..), r=Sqrt[bV]+.. *)
Vf[r_]:=r(r^2-2r)/(b r+2); Wf[r_]:=(1/Sqrt[r(r-2)])/D[Vf[r],r];
(* series solve r(V) for large V *)
rsol = r /. Solve[Vf[r]==V, r]; 
(* pick the real large branch numerically-guided: expand asymptotically instead *)
rasy = InverseSeries[Series[Vf[r]/.r->1/e,{e,0,4}]] ;
Print["  (skipping symbolic inverse; verifying correction coeff numerically below)"];
Do[
  bb=0.96; 
  Wv = N[Wf[r]/.{r->rr,b->bb}, 30];
  Vv = N[Vf[r]/.{r->rr,b->bb}, 30];
  h1 = (Wv - 1/(2 Vv)) Vv^(3/2);   (* -> -1/(2 b^{3/2}) = *)
  Print["   r=",rr," : (W-1/2V) V^{3/2} = ", N[h1,8], "   -1/(2 b^{3/2})=", N[-1/(2 bb^(3/2)),8]],
  {rr,{50,200,1000,5000}}];
Print["=== (8) FULL closed form: sqrt(x)Phi' = W(V0)(V0-2x)/sqrt(V0-x)+Int_x^V0 W'(t)(2x-t)/sqrt(t-x)dt ==="];
Print["    (identity from IBP, valid for ANY W; test with a generic polynomial-ish W(V)=V^p) "];
(* Phi(x)=sqrt(x) Int_x^V0 W(A)/sqrt(A-x)dA. Take W(A)=A^p, check derivative identity. *)
p=-3/4;  (* arbitrary exponent *)
Wg[A_]:=A^p;
PhiG[x_]:=Sqrt[x] Integrate[Wg[A]/Sqrt[A-x],{A,x,V0},Assumptions->{0<x<V0}];
lhs8=Simplify[Sqrt[x] D[PhiG[x],x],Assumptions->{0<x<V0}];
rhs8=Simplify[Wg[V0](V0-2x)/Sqrt[V0-x]+Integrate[D[Wg[t],t](2x-t)/Sqrt[t-x],{t,x,V0},Assumptions->{0<x<V0}],Assumptions->{0<x<V0}];
Print["  lhs-rhs = ", Simplify[lhs8-rhs8,Assumptions->{0<x<V0}]];

Print["=== (9) quarter-regime exact integrals ==="];
i1=Integrate[(2x-t)/Sqrt[t-x],{t,x,2x},Assumptions->{x>0}];
Print["  Int_x^2x (2x-t)/sqrt(t-x) dt = ",Simplify[i1]," (expect (4/3)x^{3/2}) diff=",Simplify[i1-(4/3)x^(3/2)]];
i2=Integrate[(t-2x)/Sqrt[t-x],{t,2x,V0},Assumptions->{0<2x<V0}];
i2t=2((V0-x)^(3/2)/3 - x Sqrt[V0-x]+(2/3)x^(3/2));
Print["  Int_2x^V0 (t-2x)/sqrt(t-x) dt diff = ",Simplify[i2-i2t,Assumptions->{0<2x<V0}]];

Print["=== (10) single-crossing algebra: at Phi'=0, L''= -3/4(a^2+c^2)+1/2 a c + I''/I ==="];
(* L=logPhi=1/2 log x+1/2 log(V0-x)+log I ; a=1/x,c=1/(V0-x). L'=1/2(a-c)+I'/I; L''=-1/2(a^2+c^2)+I''/I-(I'/I)^2. At L'=0: I'/I=-1/2(a-c). *)
a=1/x; c=1/(V0-x);
Lpp = -1/2 (a^2+c^2) + IPP - (-1/2 (a-c))^2;   (* substitute I'/I=-1/2(a-c), and I''/I=IPP symbol *)
targ = -3/4 (a^2+c^2)+1/2 a c + IPP;
Print["  Lpp - target = ", Simplify[Lpp-targ]];
