(* Independent Mathematica check for probe_conjugate_points.py: curvature of the rail optical metric
   |dx|^2/(2 Omega - C), and conjugate times from L4 (Earth-Moon CR3BP) and around Jupiter.
   Exits with status 1 if any check fails. *)
$Failures = 0;
check[label_, cond_] := (If[!TrueQ[cond], $Failures++]; Print[If[TrueQ[cond], "PASS  ", "FAIL  "], label]);
Om[x_,y_]:=(x^2+y^2)/2+(1-m)/Sqrt[(x+m)^2+y^2]+m/Sqrt[(x-1+m)^2+y^2];
f=2 Om[x,y]-cJ; ph=-Log[f]/2; K=-Exp[-2 ph](D[ph,{x,2}]+D[ph,{y,2}]);
check["curvature K = Lap Omega - 2|grad Omega|^2/(2 Omega - C)",
  Simplify[K-((D[Om[x,y],{x,2}]+D[Om[x,y],{y,2}])-2(D[Om[x,y],x]^2+D[Om[x,y],y]^2)/f)]===0];
fk=2 e+2 mu/Sqrt[x^2+y^2]; pk=-Log[fk]/2; Kk=-Exp[-2 pk](D[pk,{x,2}]+D[pk,{y,2}]);
check["Kepler: K = E mu/(r^2 (E r + mu))",
  FullSimplify[Kk-e mu/((x^2+y^2)(e Sqrt[x^2+y^2]+mu)),{x>0,y>0,mu>0,e>0}]===0];

(* conjugate times: first zero of the Jacobi field, located by bisection on the dense solution *)
firstZero[jf_, tmin_, tmax_] := Module[{ts = Subdivide[tmin, tmax, 4000], v, k},
  v = jf /@ ts; k = FirstPosition[Sign[Most[v]] Sign[Rest[v]], -1, {None}][[1]];
  If[k === None, None, t /. FindRoot[jf[t], {t, ts[[k]], ts[[k + 1]]}, Method -> "Brent", AccuracyGoal -> 14]]];

mv=0.012150585; TU=1/2.6617*^-6;
Omv[x_,y_]=Om[x,y]/.m->mv; gx[x_,y_]=D[Omv[x,y],x]; gy[x_,y_]=D[Omv[x,y],y]; lap[x_,y_]=D[Omv[x,y],{x,2}]+D[Omv[x,y],{y,2}];
x0=1/2-mv; y0=Sqrt[3]/2; c4=2 Omv[x0,y0]-1;
tcL4[deg_] := Module[{sol}, sol = NDSolve[{
   x'[s]==Sqrt[2Omv[x[s],y[s]]-c4] Cos[th[s]], y'[s]==Sqrt[2Omv[x[s],y[s]]-c4] Sin[th[s]],
   th'[s]==-(-gx[x[s],y[s]] Sin[th[s]]+gy[x[s],y[s]] Cos[th[s]])/Sqrt[2Omv[x[s],y[s]]-c4],
   j''[s]==-(lap[x[s],y[s]]-2(gx[x[s],y[s]]^2+gy[x[s],y[s]]^2)/(2Omv[x[s],y[s]]-c4)) j[s],
   x[0]==x0,y[0]==y0,th[0]==deg Degree,j[0]==0,j'[0]==1}, {x,y,th,j}, {s,0,6},
   PrecisionGoal->12, AccuracyGoal->12, MaxSteps->10^7];
  firstZero[(j /. First[sol]), 10^-2, 6] TU/86400.];
t290 = tcL4[290]; t326 = tcL4[326];
Print["L4 ray 290 deg: t_c = ", t290, " d;  L4 ray 326 deg: t_c = ", t326, " d"];
check["L4 ray 290 deg: t_c within 1e-6 d of the Python value 4.262761947", NumberQ[t290] && Abs[t290-4.262761947] < 10^-6];
check["L4 ray 326 deg: t_c within 1e-6 d of the Python value 13.347070586", NumberQ[t326] && Abs[t326-13.347070586] < 10^-6];

(* Jupiter: time in days to keep the event scale of order one *)
mj=1.26686534*^17; RJ=7.1492*^7; EJ=6000.^2/2; r0=100 RJ; thj=Pi+ArcTan[r0,RJ/2]; dsec=86400.;
solJ = NDSolve[{x'[s]==dsec Sqrt[2EJ+2mj/Sqrt[x[s]^2+y[s]^2]] Cos[th[s]], y'[s]==dsec Sqrt[2EJ+2mj/Sqrt[x[s]^2+y[s]^2]] Sin[th[s]],
  th'[s]==-dsec (mj/(x[s]^2+y[s]^2)^(3/2) (x[s] Sin[th[s]]-y[s] Cos[th[s]]))/Sqrt[2EJ+2mj/Sqrt[x[s]^2+y[s]^2]],
  j''[s]==-dsec^2 EJ mj/((x[s]^2+y[s]^2)(EJ Sqrt[x[s]^2+y[s]^2]+mj)) j[s],x[0]==r0,y[0]==0,th[0]==thj,j[0]==0,j'[0]==1},
  {x,y,th,j},{s,0,40},PrecisionGoal->12,AccuracyGoal->12,MaxSteps->10^7];
tJ = firstZero[(j /. First[solJ]), 10^-3, 40];
Print["Jupiter hyperbolic rail (v_inf 6 km/s, offset 0.5 R from 100 R): t_c = ", tJ, " d"];
check["Jupiter: t_c within 1e-6 d of the Python value 9.014412504", NumberQ[tJ] && Abs[tJ-9.014412504] < 10^-6];

Print[If[$Failures === 0, "ALL CHECKS PASSED", ToString[$Failures] <> " CHECK(S) FAILED"]];
If[$Failures =!= 0, Exit[1]];
