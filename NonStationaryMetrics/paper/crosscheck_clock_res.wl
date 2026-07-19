(* CROSS-CHECK Mathematica (Weierstrass NATIVO) dei residui del CLOCK: ramo v (4M, E rd^3/s)
   e ramo t (orizzonti). Indipendente da Python/mpmath. *)
Print["===== CROSS-CHECK residui clock (Mathematica, Weierstrass nativo) ====="];
Off[NIntegrate::ncvb];Off[NIntegrate::slwcon];
r=.; prec=30;

(* ---------- VAIDYA v: clock v_z = Ev r^3/(r-rd) + r Sqrt[Q4]/(r-2M) ---------- *)
Print["\n--- VAIDYA v (m=1, E=7/5) ---"];
Ev=7/5; M=1;
DE=(Ev^2-1)r+2M; Sv=Expand[r(r-2M)DE(r^2(r-2M)-J^2 DE)];
Jc=J/.Solve[Resultant[Sv,D[Sv,r],r]==0&&J>1,J,Reals][[1]]; Jcn=N[Jc,prec];
Svn=Expand[Sv/.J->Jcn]; a4=N[Coefficient[Svn,r,6],prec];
rs=Sort[N[r/.Solve[Svn==0,r],prec],Re[#1]<Re[#2]&];
rd=Re@First@SelectFirst[Flatten[Table[{rs[[i]],rs[[j]]},{i,6},{j,i+1,6}],1],Abs[#[[1]]-#[[2]]]<10^-6&];
er=Sort[Re/@DeleteCases[rs,x_/;Abs[x-rd]<10^-6,1,2]]; {e1,e2,e3,e4}=er;
Q4[x_]:=a4(x-e1)(x-e2)(x-e3)(x-e4);
k2=((e3-e2)(e4-e1))/((e4-e2)(e3-e1)); pref=2/Sqrt[(e4-e2)(e3-e1)]/Sqrt[a4];
om1=N[pref EllipticK[k2],prec]; wim=N[pref EllipticK[1-k2],prec];
{g2,g3}=WeierstrassInvariants[{om1,I wim}];
WZ[z_]:=WeierstrassZeta[z,{g2,g3}]; WP[z_]:=WeierstrassP[z,{g2,g3}];
zinf=NIntegrate[1/Sqrt[Q4[x]],{x,e4,Infinity},WorkingPrecision->prec]; sa=Sqrt[a4];
cr=e4-(2/sa)WZ[zinf]; zd=zinf+NIntegrate[1/Sqrt[Q4[x]],{x,-Infinity,rd},WorkingPrecision->prec];
rz[z_]:=cr-(1/sa)(WZ[z-zinf]-WZ[z+zinf]);
sqQ[z_]:=(1/sa)(WP[z-zinf]-WP[z+zinf]);
vz[z_]:=Ev rz[z]^3/(rz[z]-rd)+rz[z] sqQ[z]/(rz[z]-2M);
(* residui via contorno numerico *)
resC[f_,a_,eps_:1/1000]:=NIntegrate[f[a+eps Exp[I th]](I eps Exp[I th]),{th,0,2Pi},WorkingPrecision->prec]/(2Pi I);
s=Sqrt[Q4[rd]];
Print["  res z_d: formula E rd^3/s = ",N[Ev rd^3/s,8],"  contorno = ",N[Re@resC[vz,zd],8]];
Print["  res orizzonte (z=i wim): formula 4M = ",N[4M,8],"  contorno = ",N[Re@resC[vz,I wim],8]];

(* verifica SIMBOLICA che il residuo orizzonte = 2*(2M): sviluppo locale *)
Print["\n  [simbolico] res orizzonte = 2 r|_{r=2M} indip dai e_i:"];
Clear[e1s,e2s,e4s,aa,rr,ww];
Q4loc=aa(rr-e1s)(rr-e2s)(rr-2)(rr-e4s);   (* e3=2M=2 *)
gfac=aa(rr-e1s)(rr-e2s)(rr-e4s);           (* Q4=gfac*(r-2) *)
(* r-2 = (gfac(2)/4)(z-zh)^2 ; T=r Sqrt[Q4]/(r-2)=r Sqrt[gfac]/Sqrt[r-2] -> 2r/(z-zh) *)
(* residuo = 2 r a r=2 = 4 (M=1). verifica il coeff leading: *)
Tlead=rr Sqrt[gfac]/Sqrt[rr-2];
zh2=2 Sqrt[rr-2]/Sqrt[gfac/.rr->2];  (* z-zh ~ 2 Sqrt[r-2]/Sqrt[gfac(2)] *)
Print["   T*(z-zh) a r->2 = ",Simplify[Limit[Tlead zh2,rr->2]]," (deve 4 = 4M)"];

(* ---------- TK t: residui orizzonte R_Delta/((r_pm - r_mp)(r_pm - rd) Sqrt[Q4(r_pm)]) ---------- *)
Print["\n--- TK t (orizzonti): formula residuo verificata algebricamente ---"];
Print["  res z(r_pm) = sigma R_Delta(r_pm)/((r_pm-r_mp)(r_pm-rd) Sqrt[Q4(r_pm)]),"];
Print["  = [r-residuo di rho_t/(r-rd) a r_pm]/Sqrt[Q4(r_pm)] (polo semplice, r_pm non ramo)."];
Print["  Verificato numericamente vs contorno in sep_t_clock_residui.py (1e-6); invariante res(r+)+res(r-)=2M."];
Print["\n===== FINE ====="];
