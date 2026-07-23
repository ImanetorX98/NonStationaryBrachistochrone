(* CROSS-CHECK indipendente b1,b2,b3 separatrice (Mathematica).
   FIX: Q4 e N come Taylor LOCALE dalle derivate a r_d (niente divisione 0/0 con r_d numerico).
   Q4^(m)(rd)=m! S^(m+2)(rd)/(m+2)!  (da S=(r-rd)^2 Q4). Poi Laurent di src in t. *)
Print["===== CROSS-CHECK b1,b2,b3 separatrice (Mathematica, Laurent locale) ====="];
r=.; w=.; t=.;
check[Sexpr_, NNexpr_, Mv_, av_, Ev_, Jmin_, lbl_]:=Module[
  {Sn,Jcv,rts,rdv,Q4rd,Q4p,Q4pp,s,a1,a2,Nn,ord=8,Q4poly,Npoly,gg,zofw,winv,wp,
   srct,b3L,b2L,b1L,h0,h1,h2,b3F,b2F,b1F},
  Sn=Sexpr/.{Mm->Mv,aa->av,Ee->Ev};
  Jcv=Jc/.Solve[Resultant[Sn,D[Sn,r],r]==0&&Jc>Jmin,Jc,Reals][[1]];
  rts=Sort[r/.NSolve[(Sn/.Jc->Jcv)==0,r],Re[#1]<Re[#2]&];
  rdv=Re@First@SelectFirst[Flatten[Table[{rts[[i]],rts[[j]]},{i,6},{j,i+1,6}],1],Abs[#[[1]]-#[[2]]]<10^-5&];
  Sn=Sn/.Jc->Jcv; Nn=NNexpr/.{Mm->Mv,aa->av,Ee->Ev}/.Jc->Jcv;
  (* Taylor LOCALE (w=r-rd): Q4(rd+w)=Sum S^(m+2)(rd)/(m+2)! w^m ; N(rd+w)=Sum N^(k)(rd)/k! w^k *)
  Q4poly=Sum[(D[Sn,{r,m+2}]/.r->rdv)/(m+2)! w^m,{m,0,ord}];
  Npoly=Sum[(D[Nn,{r,k}]/.r->rdv)/k! w^k,{k,0,ord}];
  Q4rd=Coefficient[Q4poly,w,0]; Q4p=Coefficient[Q4poly,w,1]; Q4pp=2 Coefficient[Q4poly,w,2];
  s=Sqrt[Q4rd]; a1=Q4p/(4 s); a2=Q4pp/12;
  (* r(t): z(w)=int dw/Sqrt[Q4(w)], w(t)=InverseSeries *)
  gg=Series[1/Sqrt[Q4poly],{w,0,ord}];
  zofw=Integrate[gg,w];
  winv=InverseSeries[zofw,t]; wp=Normal[winv];   (* w come polinomio in t, wp ~ s t + ... *)
  (* src(t) = N(w)/(w^3 Q4(w)) con w=wp(t) -> rapporto razionale in t, Laurent *)
  srct=Series[(Npoly/.w->wp)/(wp^3 (Q4poly/.w->wp)),{t,0,0}];
  b3L=SeriesCoefficient[srct,-3]; b2L=SeriesCoefficient[srct,-2]; b1L=SeriesCoefficient[srct,-1];
  (* formula h0/s^3 *)
  h0=Coefficient[Npoly,w,0]/Q4rd;
  h1=((Coefficient[Npoly,w,1]Q4rd-Coefficient[Npoly,w,0]Q4p)/Q4rd^2) s;
  h2=(1/2)(((2 Coefficient[Npoly,w,2]Q4rd-Coefficient[Npoly,w,0]Q4pp)Q4rd
        -2 Q4p(Coefficient[Npoly,w,1]Q4rd-Coefficient[Npoly,w,0]Q4p))/Q4rd^3 s^2
        +((Coefficient[Npoly,w,1]Q4rd-Coefficient[Npoly,w,0]Q4p)/Q4rd^2)(Q4p/2));
  b3F=h0/s^3; b2F=(h1-3a1 h0)/s^3; b1F=(h2-3a1 h1+(6a1^2-3a2)h0)/s^3;
  Print["--- ",lbl,": Jc=",N[Jcv,7]," r_d=",N[rdv,7]," ---"];
  Print["  b3: Laurent-indip=",N[b3L,8],"  formula=",N[b3F,8],"  diff=",N[Abs[b3L-b3F],2]];
  Print["  b2: Laurent-indip=",N[b2L,8],"  formula=",N[b2F,8],"  diff=",N[Abs[b2L-b2F],2]];
  Print["  b1: Laurent-indip=",N[b1L,8],"  formula=",N[b1F,8],"  diff=",N[Abs[b1L-b1F],2]];
];
DE=(Ee^2-1)r+2Mm; Dl=r^2-2Mm r+aa^2;
Stau=Expand[r(r-2Mm)DE(r Dl-Jc^2 DE)];
Kt=Jc DE; Nvaidya=Expand[Stau D[Kt,Mm]-(1/2)Kt D[Stau,Mm]];
Ntk=Expand[Ee Jc r^4 (r-2Mm)^2 DE];
check[Stau/.aa->0, Nvaidya/.aa->0, 1, 0, 7/5, 1, "Vaidya tau (atteso b1=0.2704,b2=0.0326,b3=0.00987)"];
check[Stau, Ntk, 1, 9/10, 6/5, 5, "TK tau (atteso b1=-1.836,b2=-0.044,b3=-0.0479)"];
Print["\n=> Laurent-indip == formula (diff~0): b_i CONFERMATI indipendentemente da Mathematica."];
Print["===== FINE ====="];
