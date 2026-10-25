(* FULLY SYMBOLIC coefficients in (a, Ee, J), M=1 (scale, without loss of generality).
   c_k (source) and P(r) (T_alg) as RATIONAL functions of ALL the parameters, not only Ee. *)
Print["===== Symbolic coefficients in (a, Ee, J), M=1 ====="];
r=.; Ee=.; aa=.; JJ=.;
M=1;
DE=(Ee^2-1)r+2M; Dl=r^2-2M r+aa^2;
S=Expand[r(r-2M)DE(r Dl-JJ^2 DE)];
N2=Ee JJ r^4 (r-2M)^2 DE;
A5=Sum[ac[i] r^i,{i,0,5}]; Mp=Sum[cc[i] r^i,{i,0,4}];
eq=Expand[2 N2-(2 S D[A5,r]-A5 D[S,r]+2 S Mp)];
vars=Join[Table[ac[i],{i,0,5}],Table[cc[i],{i,0,4}]];
sol=Solve[CoefficientList[eq,r]==0,vars][[1]];

Print["\n--- c_k (source), rational in (a,Ee,J) ---"];
Do[Print["  c",i," = ",Together[cc[i]/.sol]],{i,0,4}];

Print["\n--- P(r) (elementary piece of T_alg), ALL coefficients in (a,Ee,J) ---"];
A5e=A5/.sol; integ=Together[A5e (r^3-2M r^2)/S];
Pquo=PolynomialQuotient[Numerator[integ],Denominator[integ],r];
Pint=Integrate[Pquo,r];
Do[Print["  [r^",k,"] = ",Together[Coefficient[Collect[Together[Pint],r],r,k]]],{k,1,3}];

Print["\n--- explicit dependence (non-zero derivative => it depends) ---"];
c2coef=Together[cc[2]/.sol];
Print["  c2 depends on a? ",Simplify[D[c2coef,aa]=!=0],"  on J? ",Simplify[D[c2coef,JJ]=!=0],
      "  on Ee? ",Simplify[D[c2coef,Ee]=!=0]];

Print["\n===== END ====="];
