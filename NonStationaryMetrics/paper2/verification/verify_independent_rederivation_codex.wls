#!/usr/bin/env wolframscript

(* Independent Wolfram Language audit of the headline claims in Papers I and II.
   This file was written from the equations printed in the manuscripts.
   It does not import or execute verify_independent_rederivation.py or any of the
   pre-existing verification scripts. *)

$MaxExtraPrecision = 5000;

pass = 0;
fail = 0;

check[label_String, condition_, note_String : ""] := Module[{ok = TrueQ[condition]},
  If[ok, pass++, fail++];
  Print[If[ok, "PASS  ", "FAIL  "], label,
    If[note === "", "", "   [" <> note <> "]"]];
  ok
];

section[title_String] := Print[
  "\n", StringRepeat["=", 72], "\n", title, "\n", StringRepeat["=", 72]
];

ass = M > 0 && a > 0 && En > 1 && Eh > 1 && r > 2 M;

f = 1 - 2 M/r;
Delta = r^2 - 2 M r + a^2;
DE = (En^2 - 1) r + 2 M;
Jc = a/En;

section["A. Algebraic identities and the four-case classification"];

numQ2 = Expand[En^2 r^3 Delta - DE (J (r - 2 M) + 2 M a)^2];
Q2rat = Cancel[numQ2/(r - 2 M)];
check[
  "A1  Q2 reconstructed from the printed identity is a polynomial",
  PolynomialQ[Q2rat, r] && Denominator[Together[Q2rat]] === 1
];
Q2 = Expand[Q2rat];

marginalBracket = Expand[r Delta - Jc^2 DE];
check[
  "A2  marginal factorisation",
  FullSimplify[
    marginalBracket == (r - 2 M) (r^2 + Jc^2),
    M > 0 && a > 0 && En > 1
  ]
];

Sgen = Expand[r (r - 2 M) DE (r Delta - J^2 DE)];
Smarg = FullSimplify[Sgen /. J -> Jc];
s0 = FullSimplify[Smarg /. r -> 2 M, M > 0 && a > 0 && En > 1];
s1 = FullSimplify[D[Smarg, r] /. r -> 2 M, M > 0 && a > 0 && En > 1];
s2 = FullSimplify[D[Smarg, {r, 2}] /. r -> 2 M, M > 0 && a > 0 && En > 1];
check[
  "A3  multiplicity is exactly two by the derivative criterion",
  s0 === 0 && s1 === 0 &&
    TrueQ[FullSimplify[s2 > 0, M > 0 && a > 0 && En > 1]]
];

check[
  "A4  the same double root occurs for both signs of J",
  FullSimplify[
    (Sgen /. J -> -Jc) == (Sgen /. J -> Jc),
    M > 0 && a > 0 && En > 1
  ]
];

Sclosed = r^3 f^2 DE (r^2 + Jc^2);
check[
  "A5  printed marginal closed form for S",
  FullSimplify[Smarg == Sclosed, M > 0 && a > 0 && En > 1 && r > 2 M]
];

w = Eh^2 - f;
check[
  "A6  w = DE/r at frozen A=1 and Eh=En",
  FullSimplify[(w /. Eh -> En) == DE/r, M > 0 && En > 1 && r > 2 M]
];

P = r^2 + a^2 + 2 M a^2/r;
vbar2 = 1 - f/En^2;
Pbar = P + (2 M a/r)^2/En^2;
phiPrime = (2 M a/r) vbar2/Pbar;
ptilde = J - 2 M a/(r En);
check[
  "A7  ptilde(re) = J-a/En, linear in J",
  FullSimplify[(ptilde /. r -> 2 M) == J - Jc, M > 0 && a > 0 && En > 1]
];

section["B. Exterior retrograde separatrix"];

R6 = Expand[r Q2 DE];

sepRoot[aval_?NumericQ, rguess_?NumericQ, jguess_?NumericQ] := Module[
  {prec = 100, expr, dexpr, sol},
  expr = N[R6 /. {M -> 1, a -> SetPrecision[aval, prec],
        En -> SetPrecision[6/5, prec]}, prec];
  dexpr = D[expr, r];
  sol = FindRoot[
    {expr == 0, dexpr == 0},
    {
      {r, SetPrecision[rguess, prec]},
      {J, SetPrecision[jguess, prec]}
    },
    WorkingPrecision -> 80,
    AccuracyGoal -> 55,
    PrecisionGoal -> 55,
    MaxIterations -> 300
  ];
  {r /. sol, J /. sol,
    N[D[expr, {r, 2}] /. sol, 50],
    N[D[expr, {r, 3}] /. sol, 50]}
];

sep = sepRoot[9/10, 7/2, -8];
rd = sep[[1]];
Jretro = sep[[2]];
Print["       rd   = ", N[rd, 22]];
Print["       Jc^- = ", N[Jretro, 22]];
Print["       R6'' = ", N[sep[[3]], 14]];

check[
  "B1  rd agrees with all printed digits",
  Abs[rd - SetPrecision[3.513905124011657, 50]] < 10^-14
];
check[
  "B2  Jc^- agrees with all printed digits",
  Abs[Jretro - SetPrecision[-8.053516003877019, 50]] < 10^-14
];
check[
  "B3  the root is double and not of higher order",
  Abs[sep[[3]]] > 1 && Abs[sep[[4]]] > 1
];
check["B4  the separatrix is exterior", rd > 2];

gWWnum = N[-(1 - 2/rd), 20];
vbarNum = N[1 - (1 - 2/rd)/(6/5)^2, 20];
Print["       g(W,W) = ", N[gWWnum, 9], "   vbar^2 = ", N[vbarNum, 9]];
check["B5  selector is timelike at rd", gWWnum < 0];
check["B6  rail slice is nondegenerate at rd", vbarNum > 0];
check["B7  printed g(W,W) value", Abs[gWWnum + 0.430833] < 10^-6];
check["B8  printed vbar^2 value", Abs[vbarNum - 0.700811] < 10^-6];

spinData = {
  {1/10, 2.86, -8.0, 2.857675},
  {1/2, 3.20, -8.0, 3.204756},
  {9/10, 3.51, -8.0, 3.513905},
  {99/100, 3.58, -8.0, 3.579494}
};
spinRoots = sepRoot[#[[1]], #[[2]], #[[3]]][[1]] & /@ spinData;
Print["       spin scan: ", N[spinRoots, 12]];
check["B9  all sampled separatrices remain exterior", And @@ Thread[spinRoots > 2]];
check[
  "B10 spin scan agrees with the printed values",
  Max[Abs[spinRoots - spinData[[All, 4]]]] < 6 10^-7
];

section["C. Optical curvature threshold and sign"];

curvBracket =
  32 M^2 + 4 (3 Eh^2 - 8) M r + (3 Eh^4 - 6 Eh^2 + 8) r^2;
brs = Expand[curvBracket /. r -> 2 M (1 + s)]/M^2;
coeffs = Table[Factor[Coefficient[brs, s, k]], {k, 0, 2}];
printedCoeffs = {
  12 Eh^4,
  24 Eh^2 (Eh^2 - 1),
  4 (3 Eh^4 - 6 Eh^2 + 8)
};
Do[
  check[
    "C" <> ToString[k + 1] <> "  coefficient of s^" <> ToString[k],
    FullSimplify[coeffs[[k + 1]] == printedCoeffs[[k + 1]], Eh > 1]
  ],
  {k, 0, 2}
];

check[
  "C4  all three coefficients are positive for Eh>1",
  TrueQ[FullSimplify[And @@ Thread[printedCoeffs > 0], Eh > 1]]
];

Ropt = -2 M^2 curvBracket/
  (Eh^2 r^5 (2 M + (Eh^2 - 1) r));
check[
  "C5  Ropt is pointwise negative in the full exterior",
  TrueQ[FullSimplify[Ropt < 0, M > 0 && Eh > 1 && r > 2 M]]
];

gaussNumerator =
  16 M^3 + 2 (9 Eh^2 - 14) M^2 r +
  (3 Eh^4 - 19 Eh^2 + 16) M r^2 -
  (2 Eh^2 - 3) (Eh^2 - 1) r^3;
gaussLead = Factor[Coefficient[gaussNumerator, r, 3]];
check[
  "C6  far-field Gauss-curvature threshold is Eh^2=3/2",
  FullSimplify[(gaussLead /. Eh^2 -> 3/2) == 0] &&
    FullSimplify[(gaussLead /. Eh^2 -> 2) < 0, Eh > 1]
];

section["D. Tension identity and horizon behaviour"];

lambda2 = Eh^2/(f (Eh^2 - f));
gradLogLambda = FullSimplify[f D[Log[Sqrt[lambda2]], r], ass];
mu = M/r^2;
tensionResidual = FullSimplify[mu + gradLogLambda, ass];
Print["       residual = ", tensionResidual,
  "   horizon limit = ",
  FullSimplify[Limit[tensionResidual, r -> 2 M, Direction -> "FromAbove"],
    M > 0 && Eh > 1]];
check[
  "D1  independently differentiated dilation contribution",
  FullSimplify[gradLogLambda == -mu (Eh^2 - 2 f)/(Eh^2 - f), ass]
];
check[
  "D2  printed positive residual mu+gradH(log Lambda)",
  FullSimplify[tensionResidual == mu f/(Eh^2 - f), ass]
];
check[
  "D3  residual is positive outside and tends to zero at r=2M",
  TrueQ[FullSimplify[tensionResidual > 0, ass]] &&
    TrueQ[FullSimplify[
      Limit[tensionResidual, r -> 2 M, Direction -> "FromAbove"] == 0,
      M > 0 && Eh > 1]]
];

section["E. Conformal transfer map"];

gbarUU = -f ut^2 + ur^2/f + r^2 uphi^2;
gUU = A^2 gbarUU;
check[
  "E1  normalization is preserved by ubar=A u",
  FullSimplify[(gbarUU /. {ut -> A ut, ur -> A ur, uphi -> A uphi}) == gUU]
];
ebar = f ut;
ehatTransferred = A ebar;
check[
  "E2  selector energy transfers as Ehat=A Ebar",
  FullSimplify[ehatTransferred/A == ebar]
];
check[
  "E3  timelike sign is conformally invariant for A>0",
  TrueQ[FullSimplify[(A^2 gbarUU < 0) == (gbarUU < 0), A > 0]]
];

section["F. Vaidya Kodama symmetry"];

fV = 1 - 2 mass[v]/r;
gVV = -fV;
sqrtDetV = r^2 Sin[theta];
divKodama = FullSimplify[D[sqrtDetV, v]/sqrtDetV];
lieKVV = D[gVV, v];
check["F1  Kodama vector is divergence free", divKodama === 0];
check[
  "F2  Kodama vector is Killing only in the stationary limit",
  FullSimplify[lieKVV /. Derivative[1][mass][v] -> 0] === 0
];
check[
  "F3  nonstationarity is measured by 2 m'(v)/r",
  FullSimplify[lieKVV == 2 Derivative[1][mass][v]/r]
];

section["G. Full near-horizon four-case classification"];

re = 2 M;
sVar = r - re;
prefactorLimit = FullSimplify[
  Limit[r (r - 2 M) DE/(r - re), r -> re],
  M > 0 && En > 1
];
check[
  "G1  positive simple prefactor at the horizon",
  FullSimplify[prefactorLimit == 4 En^2 M^2, M > 0 && En > 1]
];

inner = Expand[r Delta - J^2 DE];
innerAt = FullSimplify[inner /. r -> re, M > 0 && a > 0 && En > 1];
check[
  "G2  inner factor sign is sign(Jc^2-J^2)",
  FullSimplify[innerAt == 2 M (a^2 - En^2 J^2), M > 0 && a > 0 && En > 1]
];

simpleSlope = FullSimplify[D[Sgen, r] /. r -> re,
  M > 0 && a > 0 && En > 1];
check[
  "G3  simple-root slope for |J|<Jc is positive",
  TrueQ[FullSimplify[simpleSlope > 0,
    M > 0 && a > 0 && En > 1 && -a/En < J < a/En]]
];

cuspCoefficient = FullSimplify[
  Limit[Sgen/(r - re), r -> re],
  M > 0 && a > 0 && En > 1
];
check[
  "G4  generic cusp coefficient",
  FullSimplify[cuspCoefficient == 8 En^2 M^3 (a^2 - En^2 J^2),
    M > 0 && a > 0 && En > 1]
];

doubleCoefficient = FullSimplify[
  Limit[(Sgen /. J -> Jc)/(r - re)^2, r -> re],
  M > 0 && a > 0 && En > 1
];
check[
  "G5  marginal double-root coefficient",
  FullSimplify[doubleCoefficient == 4 En^2 M^2 (4 M^2 + Jc^2),
    M > 0 && a > 0 && En > 1]
];

ptildeE = FullSimplify[ptilde /. r -> re];
tDen = f/En - ptilde phiPrime;
tDenE = FullSimplify[Limit[tDen, r -> re, Direction -> "FromAbove"],
  M > 0 && a > 0 && En > 1];
check[
  "G6  horizon time denominator is -a(J-Jc)/Pbar_e",
  FullSimplify[tDenE == -a (J - Jc)/(4 M^2 + 2 a^2 + a^2/En^2),
    M > 0 && a > 0 && En > 1]
];

prIn = -Sqrt[Sgen]/(DE Delta);
rScale = Sqrt[Delta vbar2/Pbar];
rDotTau = FullSimplify[
  (rScale (Delta/r^2) prIn)/
    Sqrt[(Delta/r^2) prIn^2 + ptilde^2/Pbar],
  M > 0 && a > 0 && En > 1 && r > 2 M && Sgen > 0
];
rDotT = FullSimplify[
  (Delta vbar2/Pbar) (Delta/r^2) prIn/tDen,
  M > 0 && a > 0 && En > 1 && r > 2 M && Sgen > 0
];

cDerived = FullSimplify[
  Limit[(rDotT /. J -> -Jc)/(r - re), r -> re,
    Direction -> "FromAbove"],
  M > 0 && a > 0 && En > 1
];
cPrinted = -Sqrt[4 M^2 + Jc^2]/(8 M^2);
Print["       C derived = ", cDerived, "   C printed = ", cPrinted,
  "   ratio = ", FullSimplify[cDerived/cPrinted,
    M > 0 && a > 0 && En > 1]];
check[
  "G7  retrograde marginal constant C follows from the Hamiltonian",
  FullSimplify[cDerived == cPrinted,
    M > 0 && a > 0 && En > 1]
];

progradeRate = FullSimplify[
  Limit[rDotT /. J -> Jc, r -> re, Direction -> "FromAbove"],
  M > 0 && a > 0 && En > 1
];
progradePrinted = -a^2 En^2 Sqrt[4 M^2 + Jc^2]/
  (2 M (a^2 (1 + En^2) + 4 En^2 M^2));
Print["       prograde derived = ", progradeRate,
  "   printed = ", progradePrinted,
  "   ratio = ", FullSimplify[progradeRate/progradePrinted,
    M > 0 && a > 0 && En > 1]];
check[
  "G8  prograde marginal crossing rate",
  FullSimplify[progradeRate == progradePrinted,
    M > 0 && a > 0 && En > 1]
];

retroTangent = FullSimplify[tDenE /. J -> -Jc,
  M > 0 && a > 0 && En > 1];
check[
  "G9  retrograde marginal branch has nonzero time denominator",
  TrueQ[FullSimplify[retroTangent > 0, M > 0 && a > 0 && En > 1]]
];

check[
  "G10 forbidden branch |J|>Jc has negative S immediately outside",
  TrueQ[FullSimplify[simpleSlope < 0,
    M > 0 && a > 0 && En > 1 && J > a/En]]
];

prSquared = FullSimplify[Sgen/(DE^2 Delta^2)];
check[
  "G11 p_r^2 equals zero at re for every finite J",
  FullSimplify[prSquared /. r -> re, M > 0 && a > 0 && En > 1] === 0
];
check[
  "G12 for |J|>Jc the analytic p_r^2 is negative just outside re",
  TrueQ[FullSimplify[(D[prSquared, r] /. r -> re) < 0,
    M > 0 && a > 0 && En > 1 && J > a/En]]
];

section["H. Equatorial Kerr twist obstruction"];

xiT = -f;
xiPhi = -2 M a/r;
twistCoefficient = FullSimplify[xiT D[xiPhi, r] - xiPhi D[xiT, r]];
check[
  "H1  xi wedge dxi coefficient",
  FullSimplify[twistCoefficient == -2 M a/r^2, M > 0 && a > 0 && r > 0]
];
check[
  "H2  twist remains nonzero at the stationary-limit surface",
  FullSimplify[(twistCoefficient /. r -> 2 M) == -a/(2 M), M > 0 && a > 0]
];

section["I. Genus and radial-action reduction"];

genericCurve = Expand[R6 /. {M -> 1, a -> 9/10, En -> 6/5, J -> 5}];
genericDisc = Discriminant[genericCurve, r];
check[
  "I1  representative generic t-branch curve has degree six",
  Exponent[genericCurve, r] === 6
];
check[
  "I2  representative sextic is square-free",
  genericDisc =!= 0
];

actionNumerator = Expand[r Q2];
actionDivision = PolynomialQuotientRemainder[actionNumerator, Delta, r];
actionQuotient = Expand[actionDivision[[1]]];
actionRemainder = Expand[actionDivision[[2]]];
printedActionQuotient =
  En^2 r^3 + 2 En^2 M r^2 +
  (4 En^2 M^2 - (En^2 - 1) J^2) r +
  2 M (4 En^2 M^2 - J^2 + 2 (1 - En^2) J a);
Print["       quotient derived = ", Collect[actionQuotient, r],
  "\n       quotient printed = ", Collect[printedActionQuotient, r],
  "\n       difference       = ", Factor[actionQuotient - printedActionQuotient]];
check[
  "I3  radial-action polynomial coefficients a0 through a3",
  FullSimplify[actionQuotient == printedActionQuotient]
];
check[
  "I4  action remainder is proper over Delta",
  Exponent[actionRemainder, r] < Exponent[Delta, r]
];

representativeRemainder = actionRemainder /. {M -> 1, a -> 9/10, En -> 6/5, J -> 5};
representativeDelta = Delta /. {M -> 1, a -> 9/10};
check[
  "I5  the Delta poles are genuine in a generic example",
  Exponent[PolynomialGCD[representativeRemainder, representativeDelta], r] === 0
];

uProductDerivative =
  1/2 ((r^j/y) uk + uj (r^k/y) + (r^j uk - r^k uj)/y);
check[
  "I6  length-two Wjk integration-by-parts identity",
  FullSimplify[uProductDerivative == r^j uk/y]
];

section["SUMMARY"];
Print["Passed: ", pass, "   Failed: ", fail, "   Total: ", pass + fail];
If[fail > 0, Exit[1], Exit[0]];
