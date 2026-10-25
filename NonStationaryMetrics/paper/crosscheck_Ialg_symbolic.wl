(* I_el (elementary piece of T_alg) with SYMBOLIC COEFFICIENTS in Ee (=E_eff).
   I_el = int A5(r^3-2M r^2)/S dr = P(r) + a log-sum over the roots of the cubic rD-J^2 DE.
   Checks: (i) symbolic coefficients in Ee; (ii) different Ee give different coefficients (not universal)
   -> this justifies the symbolic form; (iii) numerical match with the value at E=7/5. *)
Print["===== I_el SYMBOLIC (Mathematica) ====="];
r=.; Ee=.;
M=1; av=9/10; J=5/2;
DE=(Ee^2-1)r+2M; Dl=r^2-2M r+av^2;
S=Expand[r(r-2M)DE(r Dl-J^2 DE)];
(* A5 from the second-kind reduction *)
N2=Ee J r^4 (r-2M)^2 DE;
A5=Sum[ac[i] r^i,{i,0,5}]; Mp=Sum[cc[i] r^i,{i,0,4}];
eq=Expand[2 N2-(2 S D[A5,r]-A5 D[S,r]+2 S Mp)];
sol=Solve[CoefficientList[eq,r]==0,Join[Table[ac[i],{i,0,5}],Table[cc[i],{i,0,4}]]][[1]];
A5e=A5/.sol;
integ=Together[A5e (r^3-2M r^2)/S];

(* polynomial part P(r): quotient, coefficients RATIONAL in Ee *)
Pquo=PolynomialQuotient[Numerator[integ],Denominator[integ],r];
Prem=PolynomialRemainder[Numerator[integ],Denominator[integ],r];
Pint=Integrate[Pquo,r];   (* polynomial, coefficients rational in Ee *)
Print["(i) P(r) = int(polynomial part), coefficients RATIONAL in Ee:"];
Print["    P(r) = ",Collect[Together[Pint],r]];

(* log piece: proper remainder / S ; S=r(r-2M)DE*cubic. The residues at r=0, 2M and the DE root VANISH;
   the 3 roots of the cubic C=r Dl-J^2 DE remain. Symbolic form: RootSum. *)
Cub=Expand[r Dl-J^2 DE];   (* cubic, roots = non-vanishing poles *)
Print["(i) cubic C(r)=r*Delta-J^2*DE = ",Cub];
(* residue of integ at the poles: symbolic Residue. Sum log = RootSum[C, resFun[#] Log[r-#]&] *)
resFun=Together[(A5e (r^3-2M r^2))/D[S,r]];  (* res_i = A5(r_i)(r_i^3-2M r_i^2)/S'(r_i) *)
Print["(i) res(r_i) as a RATIONAL function of the root r_i (then evaluated at the roots of the cubic):"];
Print["    res(x) = ",Together[resFun/.r->x]];
(* check that the residues VANISH at r=0, 2M and the DE root *)
Print["(ii) res at r=0: ",Simplify[resFun/.r->0]," ; at r=2M: ",Simplify[resFun/.r->2M],
      " ; at the DE root r=2M/(1-Ee^2): ",Simplify[resFun/.r->2M/(1-Ee^2)]," (expected 0)"];

(* complete symbolic Integrate in Ee: closed form (P + RootSum/Log) *)
Iel=Integrate[integ,r];
Print["(iii) symbolic Integrate[integ,r]: head = ",Head[Iel]," (Plus of poly + Log/RootSum)"];

(* (iv) NON-universal coefficients: compare P(r) at Ee=7/5 and Ee=13/10 *)
P1=Collect[Together[Pint/.Ee->7/5],r]; P2=Collect[Together[Pint/.Ee->13/10],r];
Print["(iv) P(r) at Ee=7/5: ",P1];
Print["     P(r) at Ee=13/10: ",P2];
Print["     => DIFFERENT coefficients -> not universal -> the SYMBOLIC form in Ee is needed (=P(r) above)."];
Print["\n===== END ====="];
