(* ::Package:: *)
(* CROSS-CHECK indipendente (Mathematica) delle IDENTITA' SIMBOLICHE che reggono
   tutte le separatrici adiabatiche. Esatto, non numerico. *)
Print["===== CROSS-CHECK SIMBOLICO (Mathematica) ====="];
r=.;J=.;
(* ---------- VAIDYA (a=0), parametro m ---------- *)
Print["\n--- VAIDYA (M=1, E=7/5) ---"];
Em=7/5; MvV=1;
DEv=(Em^2-1)r+2 MvV; Sv=Expand[r(r-2MvV)DEv(r^2(r-2MvV)-J^2 DEv)]; Kv=J DEv;
(* Jc = radice doppia di S *)
JcV=J/.Solve[Resultant[Sv,D[Sv,r],r]==0&&J>1,J,Reals][[1]];
Print["  Jc(Vaidya) = ",JcV," = ",N[JcV,10],"  (atteso 7.0266237)"];
(* N_m = dm(K/sqrtS) sqrtS^3 ;  verifica = S dmK - 1/2 K dmS con m simbolico *)
Mm=.; DEmm=(Em^2-1)r+2Mm; Smm=r(r-2Mm)DEmm(r^2(r-2Mm)-J^2 DEmm); Kmm=J DEmm;
NmDef=Simplify[D[Kmm/Sqrt[Smm],Mm]Smm^(3/2)]/.Mm->MvV;
NmForm=(Smm D[Kmm,Mm]-1/2 Kmm D[Smm,Mm])/.Mm->MvV;
Print["  N_m: dm(K/sqrtS)sqrtS^3 - (S dmK-1/2 K dmS) = ",Simplify[NmDef-NmForm]," (deve 0)"];
(* N_J e dJc/dm=Jc/m (scaling lineare) *)
NJv=Simplify[(Sv D[Kv,J]-1/2 Kv D[Sv,J])];
Print["  N_J(Vaidya) - (S*DE + J^2 r(r-2m) DE^3) = ",Simplify[NJv-(Sv DEv+J^2 r(r-2MvV)DEv^3)]," (deve 0)"];
(* Jc = m j(E): verifica scaling lineare *)
JcVm=J/.Solve[Resultant[#,D[#,r],r]==0&&J>1,J,Reals][[1]]&[
   Expand[r(r-2Mm)((Em^2-1)r+2Mm)(r^2(r-2Mm)-J^2((Em^2-1)r+2Mm))]];
Print["  Jc(m)/m costante? Jc(m=2)/2 - Jc(m=1) = ",Simplify[(JcVm/.Mm->2)/2-(JcVm/.Mm->1)]," (deve 0)"];

(* ---------- THAKURTA-KERR tau (a=9/10, E=6/5), parametro E_eff ---------- *)
Print["\n--- THAKURTA-KERR tau (M=1, a=9/10, E=6/5) ---"];
Et=6/5; av=9/10;
DEt=(Et^2-1)r+2; Delt=r^2-2r+av^2;
St=Expand[r(r-2)DEt(r Delt-J^2 DEt)]; Kt=J r(r-2)DEt/Delt;
JcT=J/.Solve[Resultant[r Delt-J^2 DEt,D[r Delt-J^2 DEt,r],r]==0&&J>5,J,Reals][[1]];
Print["  Jc(TK tau) = ",N[JcT,10],"  (atteso 20.327866)"];
(* N_tau = dE(K/sqrtS) sqrtS^3 = E J r^4(r-2M)^2 DE  (reproduce_reductions) *)
EE=.; DEe=(EE^2-1)r+2; Dele=r^2-2r+av^2; Se=r(r-2)DEe(r Dele-J^2 DEe); Ke=J r(r-2)DEe/Dele;
NtauDef=Simplify[D[Ke/Sqrt[Se],EE]Se^(3/2)];
Print["  N_tau - E J r^4(r-2M)^2 DE = ",Simplify[NtauDef-EE J r^4 (r-2)^2 DEe]," (deve 0)"];
(* N_J(TK) = S dJK - 1/2 K dJS = r^3(r-2M)^2 DE^2  (Delta cancella) *)
NJt=Simplify[St D[Kt,J]-1/2 Kt D[St,J]];
Print["  N_J(TK) - r^3(r-2M)^2 DE^2 = ",Simplify[NJt-r^3 (r-2)^2 DEt^2]," (deve 0; Delta cancella)"];
(* dJc/dE = -E Jc r_d/DE(r_d):  g=rDelta-J^2 DE, dJc/dE=-g_E/g_J a (r_d,Jc) *)
g=r Dele-J^2 DEe; gE=D[g,EE]; gJ=D[g,J];
Print["  dJc/dE formula -g_E/g_J = ",Simplify[-gE/gJ]," ; a r_d,J -> -E J r/DE (atteso -E Jc r_d/DE(r_d))"];

(* ---------- THAKURTA-KERR t (curva R6), due separatrici ---------- *)
Print["\n--- THAKURTA-KERR t (M=1, a=9/10, E=6/5): DUE separatrici ---"];
Q2=(2Et^2 J^2 r-Et^2 J^2 r^2-4Et^2 J av r+2Et^2 av^2 r+Et^2 av^2 r^2+Et^2 r^4+4J^2-4J^2 r+J^2 r^2-8J av+4J av r+4av^2);
JcTt=Sort[J/.Solve[Resultant[Q2,D[Q2,r],r]==0,J]//N//Select[#,Abs[Im[#]]<10^-8&&Abs[Re[#]]>5&]&//Re];
Print["  Jc(TK t) radice doppia Q2, |Jc|>5: ",JcTt];
Print["  => prograda ~ +19.089, retrograda ~ -18.671 (|Jc+| != |Jc-|, asimmetria frame-dragging)"];
Print["  Q2 lineare in J? Coefficiente di J^1 in Q2 = ",Coefficient[Q2,J,1]," (!=0 -> asimmetrico)"];

Print["\n===== FINE CROSS-CHECK SIMBOLICO ====="];
