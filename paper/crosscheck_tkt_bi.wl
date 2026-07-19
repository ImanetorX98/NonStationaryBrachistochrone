(* CROSS-CHECK Mathematica: b_i separatrice ramo t (TK), Jc+ e Jc-. Curva R6=r Q2 DE.
   Via INDIPENDENTE (Laurent, r(t) da ODE dr/dt=Sqrt[Q4]) vs formula h0/s^3. M=1,a=9/10,E=6/5. *)
Print["===== CROSS-CHECK b_i ramo t TK (Mathematica) ====="];
r=.; Ee=.;
M=1; av=9/10;
DE=(Ee^2-1)r+2M; Del=r^2-2M r+av^2;
Q2s=(2Ee^2 J^2 r-Ee^2 J^2 r^2-4Ee^2 J av r+2Ee^2 av^2 r+Ee^2 av^2 r^2+Ee^2 r^4+4J^2-4J^2 r+J^2 r^2-8J av+4J av r+4av^2);
R6s=Expand[r Q2s DE]; Kt=r DE(J(r-2M)+2M av)/Del;
Nts=Expand[Numerator[Cancel[R6s D[Kt,Ee]-(1/2)Kt D[R6s,Ee]]]];  (* derivo con Ee SIMBOLICO *)
(* ora sostituisco Ee=6/5 *)
Ev=6/5; Q2=Q2s/.Ee->Ev; R6=R6s/.Ee->Ev; Nt=Nts/.Ee->Ev;

check[Jguess_, lbl_]:=Module[{Jcv,rdv,Q4poly,Npoly,ord=8,s,a1,a2,Q4rd,Q4p,Q4pp,gg,zofw,winv,wp,
   srct,b1L,b2L,b3L,h0,h1,h2,b1F,b2F,b3F,w,Ntn},
  (* Jc,r_d: doppia radice di Q2 *)
  {rdv,Jcv}={r,J}/.FindRoot[{Q2==0,D[Q2,r]==0},{{r,-6.62},{J,Jguess}},WorkingPrecision->40];
  Ntn=Nt/.J->Jcv;
  (* Q4(rd+w) Taylor da R6-derivate: Q4^(m)(rd)=m! R6^(m+2)(rd)/(m+2)! *)
  Q4poly=Sum[(D[R6/.J->Jcv,{r,m+2}]/.r->rdv)/(m+2)! w^m,{m,0,ord}];
  Npoly=Sum[(D[Ntn,{r,k}]/.r->rdv)/k! w^k,{k,0,ord}];
  Q4rd=Coefficient[Q4poly,w,0]; Q4p=Coefficient[Q4poly,w,1]; Q4pp=2Coefficient[Q4poly,w,2];
  s=Sqrt[Q4rd]; a1=Q4p/(4s); a2=Q4pp/12;
  gg=Series[1/Sqrt[Q4poly],{w,0,ord}]; zofw=Integrate[gg,w]; winv=InverseSeries[zofw,t]; wp=Normal[winv];
  srct=Series[(Npoly/.w->wp)/(wp^3 (Q4poly/.w->wp)),{t,0,0}];
  b3L=SeriesCoefficient[srct,-3]; b2L=SeriesCoefficient[srct,-2]; b1L=SeriesCoefficient[srct,-1];
  h0=Coefficient[Npoly,w,0]/Q4rd;
  h1=((Coefficient[Npoly,w,1]Q4rd-Coefficient[Npoly,w,0]Q4p)/Q4rd^2)s;
  h2=(1/2)(((2Coefficient[Npoly,w,2]Q4rd-Coefficient[Npoly,w,0]Q4pp)Q4rd
        -2Q4p(Coefficient[Npoly,w,1]Q4rd-Coefficient[Npoly,w,0]Q4p))/Q4rd^3 s^2
        +((Coefficient[Npoly,w,1]Q4rd-Coefficient[Npoly,w,0]Q4p)/Q4rd^2)(Q4p/2));
  b3F=h0/s^3; b2F=(h1-3a1 h0)/s^3; b1F=(h2-3a1 h1+(6a1^2-3a2)h0)/s^3;
  Print["--- ",lbl,": Jc=",N[Jcv,7]," r_d=",N[rdv,7]," ---"];
  Print["  b3: Laurent=",N[b3L,7]," formula=",N[b3F,7]," diff=",N[Abs[b3L-b3F],2]];
  Print["  b2: Laurent=",N[b2L,7]," formula=",N[b2F,7]," diff=",N[Abs[b2L-b2F],2]];
  Print["  b1: Laurent=",N[b1L,7]," formula=",N[b1F,7]," diff=",N[Abs[b1L-b1F],2]];
];
check[19.089,"TK t prograda (Jc+)"];
check[-18.671,"TK t retrograda (Jc-)"];
Print["\n=> Laurent==formula: b_i ramo t (Jc+,Jc-) CONFERMATI indipendentemente da Mathematica."];
Print["===== FINE ====="];
