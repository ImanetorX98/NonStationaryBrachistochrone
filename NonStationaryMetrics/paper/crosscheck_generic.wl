(* CROSS-CHECK indipendente (Mathematica) dei coeff GENERICI genus-2: Vaidya e TK-t.
   Mathematica risolve la riduzione 2a specie autonomamente; verifica identita'=0, g_i, c_k num. *)
Print["===== CROSS-CHECK generici genus-2 (Mathematica) ====="];
r=.;

(* ---------- VAIDYA generico (a=0, param m) ---------- *)
Print["\n--- VAIDYA (param m; mm,Ee,JJ simbolici) ---"];
DEv=(Ee^2-1)r+2mm;
Sv=Expand[r(r-2mm)DEv(r^2(r-2mm)-JJ^2 DEv)];
Kv=JJ DEv; Nv=Expand[Sv D[Kv,mm]-(1/2)Kv D[Sv,mm]];
A5=Sum[ac[i]r^i,{i,0,5}]; Mp=Sum[cc[i]r^i,{i,0,4}];
eqv=Expand[2Nv-(2Sv D[A5,r]-A5 D[Sv,r]+2Sv Mp)];
solv=Solve[CoefficientList[eqv,r]==0,Join[Table[ac[i],{i,0,5}],Table[cc[i],{i,0,4}]]][[1]];
Print["  riduzione chiude? ",Length[solv]>0];
Print["  identita' 2N-(...) = ",Simplify[eqv/.solv]," (deve 0)"];
q6v=Expand[s^6 (Sv/.r->1/s)]; gv=Series[1/Sqrt[q6v],{s,0,1}];
Print["  g0 = ",Simplify[SeriesCoefficient[gv,0]]," (atteso 1/Sqrt[Ee^2-1])"];
Print["  g1 = ",Simplify[SeriesCoefficient[gv,1]]," (atteso mm(2Ee^2-3)/(Ee^2-1)^(3/2))"];
Print["  check g1: ",Simplify[SeriesCoefficient[gv,1]-mm(2Ee^2-3)/(Ee^2-1)^(3/2)]];
(* c_k numerici a mm=1,Ee=7/5,JJ=5/2 vs Python *)
np={mm->1,Ee->7/5,JJ->5/2};
Print["  c_k(num) a (1,7/5,5/2): ",Table[N[cc[i]/.solv/.np,6],{i,0,4}]];
Print["    [Python: c0=-0.531,c1=1.979,c2=-0.812,c3=-0.360,c4=0.189 (E=7/5) -- NB diverso perche' param m non E]"];

(* ---------- TK-t generico (M=1, param E) ---------- *)
Print["\n--- TK-t (M=1; aa,Ee,JJ simbolici) ---"];
M=1; DEt=(Ee^2-1)r+2M; Delt=r^2-2M r+aa^2;
Q2=(2Ee^2 JJ^2 r-Ee^2 JJ^2 r^2-4Ee^2 JJ aa r+2Ee^2 aa^2 r+Ee^2 aa^2 r^2+Ee^2 r^4+4JJ^2-4JJ^2 r+JJ^2 r^2-8JJ aa+4JJ aa r+4aa^2);
R6=Expand[r Q2 DEt]; Kt=r DEt(JJ(r-2M)+2M aa)/Delt;
Nt=Together[R6 D[Kt,Ee]-(1/2)Kt D[R6,Ee]];
Print["  N_t denominatore = ",Denominator[Cancel[Nt]]," (deve 1: Delta cancella)"];
Ntp=Expand[Numerator[Cancel[Nt]]];
eqt=Expand[2Ntp-(2R6 D[A5,r]-A5 D[R6,r]+2R6 Mp)];
solt=Solve[CoefficientList[eqt,r]==0,Join[Table[ac[i],{i,0,5}],Table[cc[i],{i,0,4}]]];
Print["  riduzione chiude? ",Length[solt]>0];
If[Length[solt]>0, solt=solt[[1]];
  Print["  identita' 2N_t-(...) = ",Simplify[eqt/.solt]," (deve 0)"]];
q6t=Expand[s^6 (R6/.r->1/s)]; gt=Series[1/Sqrt[q6t],{s,0,1}];
Print["  g0 = ",Simplify[SeriesCoefficient[gt,0]]," (atteso 1/(Ee Sqrt[Ee^2-1]))"];
Print["  g1 = ",Simplify[SeriesCoefficient[gt,1]]," (atteso -1/(Ee(Ee^2-1)^(3/2)))"];
Print["  check g0: ",Simplify[SeriesCoefficient[gt,0]-1/(Ee Sqrt[Ee^2-1])],
      " ; g1: ",Simplify[SeriesCoefficient[gt,1]+1/(Ee(Ee^2-1)^(3/2))]];
Print["\n===== FINE ====="];
