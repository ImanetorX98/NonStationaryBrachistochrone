(* Coefficienti PIENAMENTE SIMBOLICI in (a, Ee, J), M=1 (WLOG scala).
   c_k (source) e P(r) (T_alg) come funzioni RAZIONALI di TUTTI i parametri, non solo Ee. *)
Print["===== Coefficienti simbolici in (a, Ee, J), M=1 ====="];
r=.; Ee=.; aa=.; JJ=.;
M=1;
DE=(Ee^2-1)r+2M; Dl=r^2-2M r+aa^2;
S=Expand[r(r-2M)DE(r Dl-JJ^2 DE)];
N2=Ee JJ r^4 (r-2M)^2 DE;
A5=Sum[ac[i] r^i,{i,0,5}]; Mp=Sum[cc[i] r^i,{i,0,4}];
eq=Expand[2 N2-(2 S D[A5,r]-A5 D[S,r]+2 S Mp)];
vars=Join[Table[ac[i],{i,0,5}],Table[cc[i],{i,0,4}]];
sol=Solve[CoefficientList[eq,r]==0,vars][[1]];

Print["\n--- c_k (source), razionali in (a,Ee,J) ---"];
Do[Print["  c",i," = ",Together[cc[i]/.sol]],{i,0,4}];

Print["\n--- P(r) (pezzo elementare di T_alg), TUTTI i coefficienti in (a,Ee,J) ---"];
A5e=A5/.sol; integ=Together[A5e (r^3-2M r^2)/S];
Pquo=PolynomialQuotient[Numerator[integ],Denominator[integ],r];
Pint=Integrate[Pquo,r];
Do[Print["  [r^",k,"] = ",Together[Coefficient[Collect[Together[Pint],r],r,k]]],{k,1,3}];

Print["\n--- dipendenza esplicita (derivate non nulle => dipende) ---"];
c2coef=Together[cc[2]/.sol];
Print["  c2 dipende da a? ",Simplify[D[c2coef,aa]=!=0],"  da J? ",Simplify[D[c2coef,JJ]=!=0],
      "  da Ee? ",Simplify[D[c2coef,Ee]=!=0]];

Print["\n===== FINE ====="];
