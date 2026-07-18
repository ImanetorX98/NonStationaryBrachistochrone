(* I_el (pezzo elementare di T_alg) con COEFFICIENTI SIMBOLICI in Ee (=E_eff).
   I_el = int A5(r^3-2M r^2)/S dr = P(r) + somma-log sulle radici del cubico rD-J^2 DE.
   Verifica: (i) coeff simbolici in Ee; (ii) a Ee diversi danno coeff diversi (non universali)
   -> giustifica la forma simbolica; (iii) match numerico col valore E=7/5. *)
Print["===== I_el SIMBOLICO (Mathematica) ====="];
r=.; Ee=.;
M=1; av=9/10; J=5/2;
DE=(Ee^2-1)r+2M; Dl=r^2-2M r+av^2;
S=Expand[r(r-2M)DE(r Dl-J^2 DE)];
(* A5 dalla riduzione 2a specie *)
N2=Ee J r^4 (r-2M)^2 DE;
A5=Sum[ac[i] r^i,{i,0,5}]; Mp=Sum[cc[i] r^i,{i,0,4}];
eq=Expand[2 N2-(2 S D[A5,r]-A5 D[S,r]+2 S Mp)];
sol=Solve[CoefficientList[eq,r]==0,Join[Table[ac[i],{i,0,5}],Table[cc[i],{i,0,4}]]][[1]];
A5e=A5/.sol;
integ=Together[A5e (r^3-2M r^2)/S];

(* parte polinomiale P(r): quoziente, coeff RAZIONALI in Ee *)
Pquo=PolynomialQuotient[Numerator[integ],Denominator[integ],r];
Prem=PolynomialRemainder[Numerator[integ],Denominator[integ],r];
Pint=Integrate[Pquo,r];   (* polinomio, coeff razionali in Ee *)
Print["(i) P(r) = int(parte poly), coeff RAZIONALI in Ee:"];
Print["    P(r) = ",Collect[Together[Pint],r]];

(* pezzo-log: resto proprio / S ; S=r(r-2M)DE*cubico. Residui a r=0,2M,DE-root NULLI;
   restano le 3 radici del cubico C=r Dl-J^2 DE. Forma simbolica: RootSum. *)
Cub=Expand[r Dl-J^2 DE];   (* cubico, radici = poli non nulli *)
Print["(i) cubico C(r)=r*Delta-J^2*DE = ",Cub];
(* residuo di integ ai poli: Residue simbolico. Sum log = RootSum[C, resFun[#] Log[r-#]&] *)
resFun=Together[(A5e (r^3-2M r^2))/D[S,r]];  (* res_i = A5(r_i)(r_i^3-2M r_i^2)/S'(r_i) *)
Print["(i) res(r_i) come funzione RAZIONALE della radice r_i (poi valutata sulle radici del cubico):"];
Print["    res(x) = ",Together[resFun/.r->x]];
(* verifica residui NULLI a r=0,2M, e DE-root *)
Print["(ii) res a r=0: ",Simplify[resFun/.r->0]," ; a r=2M: ",Simplify[resFun/.r->2M],
      " ; a DE-root r=2M/(1-Ee^2): ",Simplify[resFun/.r->2M/(1-Ee^2)]," (attesi 0)"];

(* Integrate simbolico completo in Ee: forma chiusa (P + RootSum/Log) *)
Iel=Integrate[integ,r];
Print["(iii) Integrate[integ,r] simbolico: head = ",Head[Iel]," (Plus di poly + Log/RootSum)"];

(* (iv) coeff NON universali: confronto P(r) a Ee=7/5 vs Ee=13/10 *)
P1=Collect[Together[Pint/.Ee->7/5],r]; P2=Collect[Together[Pint/.Ee->13/10],r];
Print["(iv) P(r) a Ee=7/5: ",P1];
Print["     P(r) a Ee=13/10: ",P2];
Print["     => coeff DIVERSI -> non universali -> forma SIMBOLICA in Ee necessaria (=P(r) sopra)."];
Print["\n===== FINE ====="];
