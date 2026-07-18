(* P(r) dipende da E soltanto, o anche da a,J? Tengo M=1 (scala) e a,Ee,J SIMBOLICI. *)
Print["===== P(r): dipendenza dai parametri (a, Ee, J) ====="];
r=.; Ee=.; aa=.; JJ=.;
M=1;
DE=(Ee^2-1)r+2M; Dl=r^2-2M r+aa^2;
S=Expand[r(r-2M)DE(r Dl-JJ^2 DE)];
N2=Ee JJ r^4 (r-2M)^2 DE;
A5=Sum[ac[i] r^i,{i,0,5}]; Mp=Sum[cc[i] r^i,{i,0,4}];
eq=Expand[2 N2-(2 S D[A5,r]-A5 D[S,r]+2 S Mp)];
sol=Solve[CoefficientList[eq,r]==0,Join[Table[ac[i],{i,0,5}],Table[cc[i],{i,0,4}]]][[1]];
A5e=A5/.sol;
integ=Together[A5e (r^3-2M r^2)/S];
Pquo=PolynomialQuotient[Numerator[integ],Denominator[integ],r];
Pint=Collect[Together[Integrate[Pquo,r]],r];
c3=Simplify[Coefficient[Pint,r,3]];
Print["coeff di r^3 in P(r) (M=1, a,Ee,J simbolici):"];
Print["  ",c3];
Print["\ndipende da a? D[c3,a] == 0 ? ",Simplify[D[c3,aa]===0], "  (se False -> dipende da a)"];
Print["dipende da J? D[c3,J] == 0 ? ",Simplify[D[c3,JJ]===0], "  (se False -> dipende da J)"];
Print["dipende da Ee? D[c3,Ee]==0 ? ",Simplify[D[c3,Ee]===0]];
(* valori a due a diversi, stesso Ee,J *)
Print["\ncoeff r^3 a {Ee=7/5,J=5/2}, a=9/10: ",Simplify[c3/.{aa->9/10,JJ->5/2,Ee->7/5}]];
Print["coeff r^3 a {Ee=7/5,J=5/2}, a=1/2: ",Simplify[c3/.{aa->1/2,JJ->5/2,Ee->7/5}]];
Print["\n===== FINE ====="];
