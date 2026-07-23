(* CROSS-CHECK NUMERICO (Mathematica, WeierstrassP nativo) - TK tau separatrice.
   Valida il caso ROTANTE indipendentemente da mpmath/theta1. vs anchor Python. *)
$MinPrecision=40; prec=35;
Off[NIntegrate::ncvb];Off[NIntegrate::slwcon];Off[General::stop];
Print["===== CROSS-CHECK NUMERICO TK tau (Mathematica) ====="];
M=1; a=9/10; Eq=6/5; r0=20;
r=.;J=.;
DE=(Eq^2-1)r+2M; Delta=r^2-2M r+a^2;
Jc=J/.Solve[Resultant[r Delta-J^2 DE,D[r Delta-J^2 DE,r],r]==0&&J>5,J,Reals][[1]];
Jcn=N[Jc,prec];
Ssym=Expand[r(r-2M)DE(r Delta-J^2 DE)/.J->Jcn];
a4=N[Coefficient[Ssym,r,6],prec];
rootsS=Sort[N[r/.Solve[Ssym==0,r],prec],Re[#1]<Re[#2]&];
rd=Re@First@SelectFirst[Flatten[Table[{rootsS[[i]],rootsS[[j]]},{i,6},{j,i+1,6}],1],Abs[#[[1]]-#[[2]]]<10^-6&];
erts=Sort[Re/@DeleteCases[rootsS,x_/;Abs[x-rd]<10^-6,1,2]];
{e1,e2,e3,e4}=erts; Q4[x_]:=a4 (x-e1)(x-e2)(x-e3)(x-e4);
Print["  Jc = ",N[Jcn,10]," (Python 20.327866)"];
Print["  r_d = ",N[rd,8]," (Python -7.12951)"];
Print["  radici Q4 = ",N[erts,6]," (Python {-4.5455,0,2,16.259})"];
(* reticolo + r(z) via WeierstrassZeta nativo *)
k2=((e3-e2)(e4-e1))/((e4-e2)(e3-e1)); pref=2/Sqrt[(e4-e2)(e3-e1)]/Sqrt[a4];
om1=N[pref EllipticK[k2],prec]; wim=N[pref EllipticK[1-k2],prec];
{g2,g3}=WeierstrassInvariants[{om1,I wim}];
WZ[z_]:=WeierstrassZeta[z,{g2,g3}];
zr[rv_]:=NIntegrate[1/Sqrt[Q4[x]],{x,e4,rv},WorkingPrecision->prec];
zinf=NIntegrate[1/Sqrt[Q4[x]],{x,e4,Infinity},WorkingPrecision->prec];
sa=Sqrt[a4]; cr=e4-(2/sa)WZ[zinf];
rOfz[z_]:=cr-(1/sa)(WZ[z-zinf]-WZ[z+zinf]);
Print["  check r(z): r(z(18))-18 = ",N[rOfz[zr[18]]-18,8]," (deve ~0)"];
(* residui b_n^a (formule chiuse) con N_tau *)
Ntau=Eq Jcn r^4 (r-2M)^2 ((Eq^2-1)r+2M);
Nf[x_]:=Ntau/.r->x; Q4d[x_]:=(D[Q4[y],y]/.y->x); Q4dd[x_]:=(D[Q4[y],{y,2}]/.y->x);
Frat[x_]:=Nf[x]/Q4[x]; sP=Sqrt[Q4[rd]]; a1=Q4d[rd]/(4sP); a2=Q4dd[rd]/12;
h0=Frat[rd]; h1=(D[Frat[x],x]/.x->rd)sP; h2=(1/2)((D[Frat[x],{x,2}]/.x->rd)sP^2+(D[Frat[x],x]/.x->rd)(Q4d[rd]/2));
b1=(h2-3a1 h1+(6a1^2-3a2)h0)/sP^3; b2=(h1-3a1 h0)/sP^3; b3=h0/sP^3;
Print["  residui R_tau a z_d: b1=",N[b1,8]," b2=",N[b2,8]," b3=",N[b3,8]];
(* delta phi DIRETTO = int eta dEF dr, eta=U3-2M U2, dEF=N_tau/S^{3/2} *)
Sn[x_?NumericQ]:=Ssym/.r->x;
U[x_?NumericQ,k_]:=NIntegrate[t^k/((t-rd)Sqrt[Q4[t]]),{t,r0,x},WorkingPrecision->prec];
eta[x_?NumericQ]:=U[x,3]-2M U[x,2];
dEF[x_?NumericQ]:=(Ntau/.r->x)/Sn[x]^(3/2);
Idir[x_?NumericQ]:=NIntegrate[dEF[t]eta[t],{t,r0,x},WorkingPrecision->prec];
Print["  delta phi DIRETTO I(19) = ",N[Idir[19],10]," (Python 0.59265992)"];
Print["  I(18) = ",N[Idir[18],10]," (Python 3.7283815)"];
Print["\n===== FINE CROSS-CHECK NUMERICO TK =====" ];
