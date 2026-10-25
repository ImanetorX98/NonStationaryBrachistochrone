(* Mathematica CROSS-CHECK: TRACKING the separatrix cancels the triple pole at r_d.
   Theorem: N_tot(r_d)=0 modulo the ideal {S(r_d)=0, S'(r_d)=0} (double root).
   N_tot=N+(dJc/dl)N_J. Checked through PolynomialReduce (Groebner basis). *)
Print["===== CROSS-CHECK TRACKING (Mathematica) ====="];
r=.; rd=.; Jc=.;

(* ---------- VAIDYA (a=0, parameter m; dJc/dm=Jc/m scaling) ---------- *)
Print["\n--- VAIDYA tau tracking (dJc/dm=Jc/m) ---"];
DEv=(Ee^2-1)r+2mm;
Sv=Expand[r(r-2mm)DEv(r^2(r-2mm)-Jc^2 DEv)];
Kv=Jc DEv;
Nm=Expand[Sv D[Kv,mm]-(1/2)Kv D[Sv,mm]];
NJ=Expand[Sv D[Kv,Jc]-(1/2)Kv D[Sv,Jc]];
(* m*N_tot = m*N_m + Jc*N_J  (removes 1/m) *)
mNtot=Expand[mm Nm + Jc NJ];
(* at r=rd ; separatrix ideal {Sd,Spd} *)
Sd=Sv/.r->rd; Spd=D[Sv,r]/.r->rd; mNtotd=mNtot/.r->rd;
gb=GroebnerBasis[{Sd,Spd},{rd,Jc}];
rem=PolynomialReduce[mNtotd,gb,{rd,Jc}][[2]];
Print["  m*N_tot(r_d) mod {S(rd),S'(rd)} = ",Simplify[rem]," (must be 0 => tracking cancels the triple pole)"];

(* ---------- TK tau tracking (dJc/dE=-E Jc r_d/DE(r_d)) ---------- *)
Print["\n--- TK tau tracking (dJc/dE=-Ee Jc rd/DE(rd)) ---"];
M=1; DEt=(Ee^2-1)r+2M; Delt=r^2-2M r+aa^2;
St=Expand[r(r-2M)DEt(r Delt-Jc^2 DEt)];
Kt=Jc r(r-2M)DEt/Delt;
Ntau=Ee Jc r^4 (r-2M)^2 DEt;              (* E-derivative source *)
NJt=r^3 (r-2M)^2 DEt^2;                    (* N_J *)
DErd=(Ee^2-1)rd+2M;
dJcdE=-Ee Jc rd/DErd;
(* N_tot = Ntau + dJcdE NJt ; at r=rd ; multiply by DE(rd) to remove the denominator *)
NtotD=Expand[(Ntau/.r->rd)DErd + (-Ee Jc rd)(NJt/.r->rd)];
Sdt=St/.r->rd; Spdt=D[St,r]/.r->rd;
gbt=GroebnerBasis[{Sdt,Spdt},{rd,Jc}];
remt=PolynomialReduce[NtotD,gbt,{rd,Jc}][[2]];
Print["  DE(rd)*N_tot(r_d) mod {S(rd),S'(rd)} = ",Simplify[remt]," (must be 0)"];

Print["\n=> if both are 0: TRACKING cancels the triple pole at r_d SYMBOLICALLY (S'(rd)=0). b3^track=0."];
Print["===== END ====="];
