(* CROSS-CHECK Mathematica dei rami v (Vaidya) e t (TK, 2 separatrici):
   delta phi DIRETTO (ground truth) + mattoni nuovi (sqrtQ4=dr/dz; rho_t; z(r_pm)). vs Python. *)
$MinPrecision=35; prec=30;
Off[NIntegrate::ncvb];Off[NIntegrate::slwcon];Off[General::stop];Off[NIntegrate::inumr];
Print["===== CROSS-CHECK v (Vaidya) + t (TK) ====="];
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
(* mattone nuovo v: sqrtQ4(z) = dr/dz = (1/sqrt a4)[P(z-zinf)-P(z+zinf)] *)
 ztest=zr[11]; Print["  Jc=",N[Jcvn,8]," (7.026624); r_d=",N[rdv,6]," (-3.36371)"];
Print["  sqrtQ4=dr/dz check: (1/sqrt a4)(P(z-zinf)-P(z+zinf)) - sqrtQ4(r(z)) = ",
   N[(1/sa)(WP[ztest-zinf]-WP[ztest+zinf])-Sqrt[Q4[11]],8]," (deve ~0)"];
(* delta phi_v DIRETTO = int dm F * v dr, v=E U3 + r + 2m ln(r-2m) *)
Sm[x_?NumericQ,mv_]:=x(x-2mv)((Ev^2-1)x+2mv)(x^2(x-2mv)-Jcvn^2((Ev^2-1)x+2mv));
dmF[x_?NumericQ]:=Module[{h=10^-8},((Jcvn((Ev^2-1)x+2(Mv+h)))/Sqrt[Sm[x,Mv+h]]-(Jcvn((Ev^2-1)x+2(Mv-h)))/Sqrt[Sm[x,Mv-h]])/(2h)];
Uv[x_?NumericQ,k_]:=NIntegrate[t^k/((t-rdv)Sqrt[Q4[t]]),{t,r0v,x},WorkingPrecision->prec];
vcl[x_?NumericQ]:=Ev Uv[x,3]+(x-r0v)+2Mv(Log[Abs[x-2Mv]]-Log[Abs[r0v-2Mv]]);
dphiv[x_?NumericQ]:=NIntegrate[dmF[t]vcl[t],{t,r0v,x},WorkingPrecision->prec];
Print["  delta phi_v DIRETTO (11) = ",N[dphiv[11],9]," (Python 0.09066210)"];

(* ============ TK t (coordinate time), M=1 a=9/10 E=6/5, DUE separatrici ============ *)
Print["\n--- TK t (M=1, a=9/10, E=6/5) ---"];
Et=6/5; at=9/10; r0t=20;
DEt=(Et^2-1)r+2; Delt=r^2-2r+at^2;
Q2t=2Et^2 J^2 r-Et^2 J^2 r^2-4Et^2 J at r+2Et^2 at^2 r+Et^2 at^2 r^2+Et^2 r^4+4J^2-4J^2 r+J^2 r^2-8J at+4J at r+4at^2;
R6=Expand[r Q2t DEt]; Kt=r DEt(J(r-2)+2 at)/Delt;
rhoC=Cancel[Together[(Et^2 r^3-2 at Kt/r)/((r-2)/r)]]; numR=Numerator[rhoC];
P3=Expand[PolynomialQuotient[numR,Delt,r]]; RD=Expand[PolynomialRemainder[numR,Delt,r]];  (* rho = P3 + RD/Delta *)
Print["  Denominator[Cancel rho] = ",Denominator[rhoC]," (deve = Delta)"];
Print["  rho_t = P3 + R_Delta/Delta ? P3 + RD/Delta - rho = ",Simplify[P3+RD/Delt-rho]," (deve 0)"];
Do[
  {Jc0,rg,tag,pyanchor}=cfg; Jctn=N[Jc0,prec];
  R6n=Expand[R6/.J->Jctn]; a4t=N[Coefficient[R6n,r,6],prec];
  (* raffina r_d,Jc con doppia radice di Q2 *)
  sol=FindRoot[{(Q2t/.J->Jj)==0,(D[Q2t,r]/.J->Jj)==0},{{r,rg},{Jj,Jc0}},WorkingPrecision->prec];
  rdt=r/.sol; Jctn=Jj/.sol; R6n=Expand[R6/.J->Jctn]; a4t=N[Coefficient[R6n,r,6],prec];
  Qt4=Cancel[R6n/(r-rdt)^2];
  et=Sort[N[r/.Solve[Qt4==0,r],prec],Re[#1]<Re[#2]&]; {f1,f2,f3,f4}=et;
  Q4t[x_]:=a4t(x-f1)(x-f2)(x-f3)(x-f4);
  k2t=((f3-f2)(f4-f1))/((f4-f2)(f3-f1)); preft=2/Sqrt[(f4-f2)(f3-f1)]/Sqrt[a4t];
  o1=N[preft EllipticK[k2t],prec]; wi=N[preft EllipticK[1-k2t],prec];
  {gg2,gg3}=WeierstrassInvariants[{o1,I wi}]; WZt[z_]:=WeierstrassZeta[z,{gg2,gg3}];
  zit=NIntegrate[1/Sqrt[Q4t[x]],{x,f4,Infinity},WorkingPrecision->prec]; sat=Sqrt[a4t]; crt=f4-(2/sat)WZt[zit];
  rzt[z_]:=crt-(1/sat)(WZt[z-zit]-WZt[z+zit]);
  rpm=1+Sqrt[1-at^2]; zrp=z/.FindRoot[rzt[z]==rpm,{z,I wi/2},WorkingPrecision->prec];
  P3n=P3/.J->Jctn; RDn=RD/.J->Jctn;
  Sn6[x_?NumericQ]:=R6n/.r->x;
  dEFt[x_?NumericQ]:=(Et Jctn x^4 (x-2)^2 ((Et^2-1)x+2))/Sn6[x]^(3/2);  (* N_t/R6^{3/2} *)
  Ut[x_?NumericQ,k_]:=NIntegrate[t^k/((t-rdt)Sqrt[Q4t[t]]),{t,r0t,x},WorkingPrecision->prec];
  etat[x_?NumericQ]:=NIntegrate[((P3n/.r->tt)+(RDn/.r->tt)/(tt^2-2tt+at^2))/((tt-rdt)Sqrt[Q4t[tt]]),{tt,r0t,x},WorkingPrecision->prec];
  dphit[x_?NumericQ]:=NIntegrate[dEFt[t]etat[t],{t,r0t,x},WorkingPrecision->prec];
  Print["  ",tag,": Jc=",N[Jctn,9]," r_d=",N[rdt,6]," z(r+)=",N[zrp,5]];
  Print["       delta phi_t DIRETTO (19) = ",N[dphit[19],9]," (Python ",pyanchor,")"],
  {cfg,{{19.089,-6.62,"prograda","8.71960819"},{-18.671,-6.588,"retrograda","-7.95252758"}}}];
Print["\n===== FINE ====="];
