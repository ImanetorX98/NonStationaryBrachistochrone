(* Independent CROSS-CHECK (Mathematica) of the genus-2 closed form (TK tau, generic J).
   SYMBOLIC check of the coefficients and of the reductions. Parameters M=1, a=9/10, J=5/2, E symbolic. *)
Print["===== CROSS-CHECK genus-2 (Mathematica) ====="];
r=.; E=.;
M=1; av=9/10; J=5/2;
DE=(E^2-1)r+2M; Dl=r^2-2M r+av^2;
S=Expand[r(r-2M)DE(r Dl-J^2 DE)];
K=J r(r-2M)DE/Dl; F=K/Sqrt[S];

(* (1) N_tau = dE(K/sqrtS) sqrtS^3 = E J r^4 (r-2M)^2 DE  (Delta cancels) *)
Ntau=Simplify[D[F,E] S^(3/2)];
Print["(1) N_tau - E J r^4(r-2M)^2 DE = ",Simplify[Ntau - E J r^4 (r-2M)^2 DE]," (must be 0)"];

(* (2) second-kind reduction: 2 N = 2 S A5' - A5 S' + 2 S sum c_k r^k, A5 of degree 5, c_k rational in E *)
N2=E J r^4 (r-2M)^2 DE;
A5=Sum[ac[i] r^i,{i,0,5}]; Mp=Sum[cc[i] r^i,{i,0,4}];
eq=Expand[2 N2 - (2 S D[A5,r] - A5 D[S,r] + 2 S Mp)];
sol=Solve[CoefficientList[eq,r]==0,Join[Table[ac[i],{i,0,5}],Table[cc[i],{i,0,4}]]][[1]];
Print["(2) is the second-kind reduction solvable? ",Length[sol]>0," ; c_k rational in E:"];
Do[Print["    c",i," = ",Together[cc[i]/.sol]],{i,0,4}];

(* (3) g_i = Taylor coefficients of q6^{-1/2} at s=0, q6=s^6 S(1/s) *)
q6=Expand[s^6 (S/.r->1/s)];
gser=Series[1/Sqrt[q6],{s,0,2}];
Print["(3) g_0 = ",Simplify[SeriesCoefficient[gser,0]]," (expected 1/Sqrt[E^2-1])"];
Print["    g_1 = ",Simplify[SeriesCoefficient[gser,1]]," (expected (2E^2-3)/(E^2-1)^(3/2))"];
Print["    check g_0-1/Sqrt[E^2-1]=",Simplify[SeriesCoefficient[gser,0]-1/Sqrt[E^2-1]]];
Print["    check g_1-(2E^2-3)/(E^2-1)^(3/2)=",Simplify[SeriesCoefficient[gser,1]-(2E^2-3)/(E^2-1)^(3/2)]];

(* (4) elementary T_alg: the integrand A5(r^3-2M r^2)/S is RATIONAL (no sqrt) -> elementary integral *)
A5e=A5/.sol;
integ=Together[A5e (r^3-2M r^2)/S];
Print["(4) is the T_alg integrand rational (denominator without sqrt)? ",
   PolynomialQ[Numerator[integ],r]&&PolynomialQ[Denominator[integ],r]];
IntTalg=Integrate[integ/.E->7/5,r];
Print["    is the integral elementary (polynomial+log)? heads of the terms: ",
   Union[Head/@Level[IntTalg,{1}]]," (Log/atan/rational = elementary)"];

(* (5) Q_kj = c_k b_j - c_j b_i, b=(0,0,-2M,1,0); the 7 non-zero ones *)
bb={0,0,-2M,1,0}; cv=Table[cc[i]/.sol,{i,0,4}];
Qkj=Table[Together[cv[[k+1]] bb[[j+1]]-cv[[j+1]] bb[[k+1]]],{k,0,4},{j,0,4}];
nz=Select[Flatten[Table[{k,j,Qkj[[k+1,j+1]]},{k,0,4},{j,k+1,4}],1],#[[3]]=!=0&];
Print["(5) non-zero Q_kj (k,j): ",nz[[All,{1,2}]]," (expected {0,2},{0,3},{1,2},{1,3},{2,3},{2,4},{3,4})"];
Print["    Q02+2M c0 = ",Simplify[(Qkj[[1,3]])-(-2M cv[[1]])]," ; Q34+c4 = ",Simplify[Qkj[[4,5]]+cv[[5]]]];

Print["\n===== END OF THE genus-2 CROSS-CHECK ====="];
