(* Mathematica CROSS-CHECK of the t branch of Thakurta-Kerr, prograde root (M=1, a=9/10, E=6/5, r0=20).
   At the double root of the t-branch quartic Q2 found from (r, J) = (-6.62, 19.089):
     (1) rho_t = P3 + RD/Delta exactly (the t-branch drift clock of Paper II, eq. clock-t);
     (2) z(r_+) in the Weierstrass uniformisation of the quartic Q4 = R6/(r - r_d)^2, with |r(z(r_+)) - r_+| < 1e-20;
     (3) the direct energy-response integral delta phi_t(x) = int_{r0}^{x} (d_E F_t) eta dr, eta' = rho_t/sqrt(R6)
         (at fixed J: one ingredient of the first-order correction, not the whole non-autonomous correction),
         at x = 19 and 18 by nested quadrature (d_E F_t with the t-branch numerator, E symbolic), against the
         independent Python references of crosscheck_t_reference.py (ODE and mpmath quadrature), relative
         tolerance 1e-12.
   Each check prints PASS or FAIL; the script exits with status 1 if any fails.
   An earlier version compared (3) with hand-copied values 100 times too large, from an unidentified producer.
   Note: r_d < 0, so this J is an algebraic genus degeneration (J_deg), not a physical separatrix
   (separatrix_classification.py; Paper II, Proposition prop:exterior-sep); the integrals are real on the
   exterior arc [18, 20].  Run: wolframscript -file crosscheck_tkt_flat.wl *)
$MinPrecision=35; prec=30;
Off[NIntegrate::ncvb];Off[NIntegrate::slwcon];Off[General::stop];Off[NIntegrate::inumr];
Print["=== TK t prograde: direct delta phi (Mathematica) vs the Python references ==="];
fails=0; verdict[ok_]:=(If[!TrueQ[ok],fails++]; If[TrueQ[ok],"PASS","FAIL"]);
r=.;J=.; Et=6/5; at=9/10; r0t=20;
DEt=(Et^2-1)r+2; Delt=r^2-2r+at^2;
Q2t=2Et^2 J^2 r-Et^2 J^2 r^2-4Et^2 J at r+2Et^2 at^2 r+Et^2 at^2 r^2+Et^2 r^4+4J^2-4J^2 r+J^2 r^2-8J at+4J at r+4at^2;
R6=Expand[r Q2t DEt]; Kt=r DEt(J(r-2)+2 at)/Delt;
(* N_t (t-branch) = dE(Kt/sqrtR6) R6^{3/2}, EE symbolic *)
DEe=(EE^2-1)r+2; Dele=r^2-2r+at^2;
Q2e=2EE^2 J^2 r-EE^2 J^2 r^2-4EE^2 J at r+2EE^2 at^2 r+EE^2 at^2 r^2+EE^2 r^4+4J^2-4J^2 r+J^2 r^2-8J at+4J at r+4at^2;
R6e=r Q2e DEe; Kte=r DEe(J(r-2)+2 at)/Dele;
Ntsym=Expand[Simplify[D[Kte/Sqrt[R6e],EE]R6e^(3/2)]];
rho=Together[(Et^2 r^3-2 at Kt/r)/((r-2)/r)];
P3=Expand[PolynomialQuotient[Numerator[rho],Denominator[rho],r]];
RD=Expand[Cancel[(rho-P3)Delt]];
resrho=Simplify[P3+RD/Delt-rho]; Print["  rho = P3 + RD/Delta ? ",resrho," (must be 0)  ",verdict[resrho===0]];
sol=FindRoot[{Q2t/.J->Jj,D[Q2t,r]/.J->Jj},{{r,-6.62},{Jj,19.089}},WorkingPrecision->prec];
rdt=r/.sol; Jc=Jj/.sol; Print["  Jc=",N[Jc,9]," r_d=",N[rdt,7]];
R6n=Expand[R6/.J->Jc]; a4t=N[Coefficient[R6n,r,6],prec];
rootsR6=Sort[N[r/.Solve[R6n==0,r],prec],Re[#1]<Re[#2]&];
et=Sort[Re/@DeleteCases[rootsR6,x_/;Abs[x-rdt]<1/10,1,2]]; {f1,f2,f3,f4}=et; Q4[x_]:=a4t(x-f1)(x-f2)(x-f3)(x-f4);
Print["  Q4 roots=",N[et,6]];
k2=((f3-f2)(f4-f1))/((f4-f2)(f3-f1)); pf=2/Sqrt[(f4-f2)(f3-f1)]/Sqrt[a4t];
o1=N[pf EllipticK[k2],prec]; wi=N[pf EllipticK[1-k2],prec]; {gg2,gg3}=WeierstrassInvariants[{o1,I wi}];
WZ[z_]:=WeierstrassZeta[z,{gg2,gg3}]; zit=NIntegrate[1/Sqrt[Q4[x]],{x,f4,Infinity},WorkingPrecision->prec];
sat=Sqrt[a4t]; crt=f4-(2/sat)WZ[zit]; rzt[z_]:=crt-(1/sat)(WZ[z-zit]-WZ[z+zit]);
rp=1+Sqrt[1-at^2]; zrp=z/.FindRoot[rzt[z]==rp,{z,0.16+0.57 I},WorkingPrecision->prec];
reszr=Abs[rzt[zrp]-rp]; Print["  z(r+)=",N[zrp,6],"  |r(z(r+))-r+| = ",N[reszr,4],"  ",verdict[NumericQ[reszr]&&reszr<10^-20]];
P3n=P3/.J->Jc; RDn=RD/.J->Jc;
Sn6[x_?NumericQ]:=R6n/.r->x;
Ntc=Ntsym/.EE->Et/.J->Jc; Ntf[x_?NumericQ]:=Ntc/.r->x;
dEFt[x_?NumericQ]:=Ntf[x]/Sn6[x]^(3/2);
etat[x_?NumericQ]:=NIntegrate[((P3n/.r->tt)+(RDn/.r->tt)/(tt^2-2tt+at^2))/((tt-rdt)Sqrt[Q4[tt]]),{tt,r0t,x},WorkingPrecision->prec];
dphit[x_?NumericQ]:=NIntegrate[dEFt[t]etat[t],{t,r0t,x},WorkingPrecision->prec];
(* references: crosscheck_t_reference.py, prograde root *)
cmp[x_,ref_]:=Module[{v=dphit[x],rel},rel=If[NumericQ[v],Abs[v-ref]/Abs[ref],Infinity];
  Print["  delta phi_t DIRECT (",x,") = ",N[v,16],"   Python reference ",ref,"   rel. diff ",N[rel,3],"  ",verdict[NumericQ[v]&&rel<10^-12]]];
cmp[19,0.08719608191295993`20]; cmp[18,0.3993812410101394`20];
Print[If[fails==0,"ALL CHECKS PASSED",ToString[fails]<>" CHECK(S) FAILED"]];
Print["=== END ==="];
Exit[If[fails==0,0,1]];
