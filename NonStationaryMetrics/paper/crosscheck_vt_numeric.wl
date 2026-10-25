(* Mathematica CROSS-CHECK of the v (Vaidya) and t (Thakurta-Kerr, two double roots) branches:
   DIRECT first-order responses (ground truth) and new building blocks, against the independent Python
   references of crosscheck_t_reference.py (ODE and mpmath quadrature), relative tolerance 1e-12.
   Vaidya v (M=1, E=7/5, r0=12): sqrtQ4 = dr/dz through Weierstrass P; delta phi_v(11) = int (d_m F) v dr.
   Thakurta-Kerr t (M=1, a=9/10, E=6/5, r0=20): the drift clock rho_t = P3 + R_Delta/Delta of Paper II,
   eq. clock-t (denominator a constant multiple of Delta, P3 equal to the printed one); then, at the prograde
   and retrograde double roots of Q2, delta phi_t(x) = int_{r0}^{x} (d_E F_t) eta dr, eta' = rho_t/sqrt(R6),
   at x = 19 and 18, with the t-branch numerator N_t = (d_E K) R6 - K (d_E R6)/2 (E symbolic, J fixed).
   Each check prints PASS or FAIL; the script exits with status 1 if any fails.
   An earlier version of the t block used the tau-branch numerator, lost the constant factor of the
   denominator of rho_t, compared with an undefined symbol and left the final integral unevaluated.
   Both double roots have r_d < 0: algebraic genus degenerations, not physical separatrices
   (separatrix_classification.py; Paper II, Proposition prop:exterior-sep); the integrals are real on the
   exterior arcs used.  Run: wolframscript -file crosscheck_vt_numeric.wl *)
$MinPrecision=35; prec=30;
Off[NIntegrate::ncvb];Off[NIntegrate::slwcon];Off[General::stop];Off[NIntegrate::inumr];
Print["===== CROSS-CHECK v (Vaidya) + t (TK) ====="];
fails=0; verdict[ok_]:=(If[!TrueQ[ok],fails++]; If[TrueQ[ok],"PASS","FAIL"]);
r=.;J=.;

(* ============ VAIDYA v (advanced time), M=1 E=7/5 ============ *)
Print["\n--- VAIDYA v (M=1, E=7/5) ---"];
Ev=7/5; Mv=1; r0v=12;
DEv=(Ev^2-1)r+2Mv; Sv=Expand[r(r-2Mv)DEv(r^2(r-2Mv)-J^2 DEv)];
Jcv=J/.Solve[Resultant[Sv,D[Sv,r],r]==0&&J>1,J,Reals][[1]]; Jcvn=N[Jcv,prec];
Svn=Expand[Sv/.J->Jcvn]; a4v=N[Coefficient[Svn,r,6],prec];
rsv=Sort[N[r/.Solve[Svn==0,r],prec],Re[#1]<Re[#2]&];
rdv=Re@First@SelectFirst[Flatten[Table[{rsv[[i]],rsv[[j]]},{i,6},{j,i+1,6}],1],Abs[#[[1]]-#[[2]]]<10^-6&];
ev=Sort[Re/@DeleteCases[rsv,x_/;Abs[x-rdv]<10^-6,1,2]]; {e1,e2,e3,e4}=ev; a4=a4v;
Q4[x_]:=a4(x-e1)(x-e2)(x-e3)(x-e4);
k2=((e3-e2)(e4-e1))/((e4-e2)(e3-e1)); pref=2/Sqrt[(e4-e2)(e3-e1)]/Sqrt[a4];
om1=N[pref EllipticK[k2],prec]; wim=N[pref EllipticK[1-k2],prec];
{g2,g3}=WeierstrassInvariants[{om1,I wim}]; WZ[z_]:=WeierstrassZeta[z,{g2,g3}]; WP[z_]:=WeierstrassP[z,{g2,g3}];
zinf=NIntegrate[1/Sqrt[Q4[x]],{x,e4,Infinity},WorkingPrecision->prec]; sa=Sqrt[a4]; cr=e4-(2/sa)WZ[zinf];
rz[z_]:=cr-(1/sa)(WZ[z-zinf]-WZ[z+zinf]); zr[rv_]:=NIntegrate[1/Sqrt[Q4[x]],{x,e4,rv},WorkingPrecision->prec];
(* new v building block: sqrtQ4(z) = dr/dz = (1/sqrt a4)[P(z-zinf)-P(z+zinf)] *)
 ztest=zr[11]; Print["  Jc=",N[Jcvn,8]," (7.026624); r_d=",N[rdv,6]," (-3.36371)"];
Print["  sqrtQ4=dr/dz check: (1/sqrt a4)(P(z-zinf)-P(z+zinf)) - sqrtQ4(r(z)) = ",
   N[(1/sa)(WP[ztest-zinf]-WP[ztest+zinf])-Sqrt[Q4[11]],8]," (must be ~0)"];
(* DIRECT delta phi_v = int dm F * v dr, v=E U3 + r + 2m ln(r-2m) *)
Sm[x_?NumericQ,mv_]:=x(x-2mv)((Ev^2-1)x+2mv)(x^2(x-2mv)-Jcvn^2((Ev^2-1)x+2mv));
dmF[x_?NumericQ]:=Module[{h=10^-8},((Jcvn((Ev^2-1)x+2(Mv+h)))/Sqrt[Sm[x,Mv+h]]-(Jcvn((Ev^2-1)x+2(Mv-h)))/Sqrt[Sm[x,Mv-h]])/(2h)];
Uv[x_?NumericQ,k_]:=NIntegrate[t^k/((t-rdv)Sqrt[Q4[t]]),{t,r0v,x},WorkingPrecision->prec];
vcl[x_?NumericQ]:=Ev Uv[x,3]+(x-r0v)+2Mv(Log[Abs[x-2Mv]]-Log[Abs[r0v-2Mv]]);
dphiv[x_?NumericQ]:=NIntegrate[dmF[t]vcl[t],{t,r0v,x},WorkingPrecision->prec];
vv=dphiv[11]; refv=0.0906620968390034`20; relv=If[NumericQ[vv],Abs[vv-refv]/Abs[refv],Infinity];
Print["  DIRECT delta phi_v (11) = ",N[vv,16],"   Python reference ",refv,"   rel. diff ",N[relv,3],"  ",verdict[NumericQ[vv]&&relv<10^-12]];

(* ============ TK t (coordinate time), M=1 a=9/10 E=6/5, the two double roots of Q2 ============ *)
Print["\n--- TK t (M=1, a=9/10, E=6/5): direct delta phi at the two double roots of Q2 ---"];
Et=6/5; at=9/10; r0t=20;
DEt=(Et^2-1)r+2; Delt=r^2-2r+at^2;
Q2t=2Et^2 J^2 r-Et^2 J^2 r^2-4Et^2 J at r+2Et^2 at^2 r+Et^2 at^2 r^2+Et^2 r^4+4J^2-4J^2 r+J^2 r^2-8J at+4J at r+4at^2;
R6=Expand[r Q2t DEt]; Kt=r DEt(J(r-2)+2 at)/Delt;
(* t-branch numerator of d_E F_t = N_t/R6^(3/2): E symbolic, J fixed *)
EE=.; DEe=(EE^2-1)r+2;
Q2e=2EE^2 J^2 r-EE^2 J^2 r^2-4EE^2 J at r+2EE^2 at^2 r+EE^2 at^2 r^2+EE^2 r^4+4J^2-4J^2 r+J^2 r^2-8J at+4J at r+4at^2;
R6e=r Q2e DEe; Kte=r DEe(J(r-2)+2 at)/Delt;
Ntsym=Together[D[Kte,EE] R6e-Kte D[R6e,EE]/2];
(* drift clock rho_t (Paper II, eq. clock-t) = P3 + R_Delta/Delta *)
rhoC=Together[(Et^2 r^3-2 at Kt/r)/((r-2)/r)];
ratio=Cancel[Denominator[rhoC]/Delt];
Print["  Denominator[rho_t]/Delta = ",ratio,"  (constant: the pole r = 2M is removable)  ",verdict[FreeQ[ratio,r]]];
P3=PolynomialQuotient[Numerator[rhoC],Denominator[rhoC],r]; RD=Cancel[(rhoC-P3)Delt];
res=Simplify[P3+RD/Delt-rhoC];
Print["  rho_t - (P3 + R_Delta/Delta) = ",res,", R_Delta of degree ",Exponent[RD,r],"  ",verdict[res===0&&PolynomialQ[RD,r]&&Exponent[RD,r]==1]];
resP=Expand[P3-(Et^2 r^3+2Et^2 r^2+4Et^2 r+8Et^2-2(Et^2-1)J at)];
Print["  P3 - (P3 printed in Paper II) = ",resP,"  ",verdict[resP===0]];
Do[
  {Jg,rg,tag,ref19,ref18}=cfg;
  sol=FindRoot[{(Q2t/.J->Jj)==0,(D[Q2t,r]/.J->Jj)==0},{{r,rg},{Jj,Jg}},WorkingPrecision->prec];
  rdt=r/.sol; Jctn=Jj/.sol;
  R6n=Expand[R6/.J->Jctn]; rhon=rhoC/.J->Jctn; Ntn=Ntsym/.EE->Et/.J->Jctn;
  etat[x_?NumericQ]:=NIntegrate[(rhon/.r->tt)/Sqrt[R6n/.r->tt],{tt,r0t,x},WorkingPrecision->prec];
  dEFt[x_?NumericQ]:=(Ntn/.r->x)/(R6n/.r->x)^(3/2);
  dphit[x_?NumericQ]:=NIntegrate[dEFt[t] etat[t],{t,r0t,x},WorkingPrecision->prec];
  Print["  ",tag,": Jc=",N[Jctn,12]," r_d=",N[rdt,9],"   R6 > 0 on [18,20]: ",
    verdict[And@@Table[(R6n/.r->y)>0,{y,18,20,1/10}]]];
  Do[{xx,ref}=xr; vt=dphit[xx]; rel=If[NumericQ[vt],Abs[vt-ref]/Abs[ref],Infinity];
    Print["       DIRECT delta phi_t (",xx,") = ",N[vt,16],"   Python reference ",ref,"   rel. diff ",N[rel,3],"  ",
      verdict[NumericQ[vt]&&rel<10^-12]],
    {xr,{{19,ref19},{18,ref18}}}],
  {cfg,{{19.089,-6.62,"prograde",0.08719608191295993`20,0.3993812410101394`20},
        {-18.671,-6.588,"retrograde",-0.07952527575427483`20,-0.3613357695835678`20}}}];
Print[If[fails==0,"ALL CHECKS PASSED",ToString[fails]<>" CHECK(S) FAILED"]];
Print["\n===== END ====="];
Exit[If[fails==0,0,1]];
