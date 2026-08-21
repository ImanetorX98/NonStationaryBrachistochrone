# Editorial notes — Applications

## FLRW

FLRW works well as a sanity check and degenerate base case. The physical message is simple: spatial homogeneity leaves a centred circular indicatrix, all displayed clock branches share the straight comoving spatial path, and expansion produces freezing without branch splitting.

The section is so short that two full-page-scale figures are disproportionate. Combine the kinematic curves and the variational check into one multi-panel figure, and qualify “minimizes both functionals” as a displayed perturbation-family check unless a global certificate is explicitly invoked.

## Vaidya

This is the strongest application and should dominate the main narrative. The sequence Kodama invariant → indicatrix/Hamiltonians → costate/mechanical map → plunge/timing/bounce → first-order response is compelling.

The presentation currently alternates between three levels:

- exact controlled-rail structure;
- numerical/frozen phenomenology;
- highly technical algebraic-curve closure.

Signal these levels visually and verbally. A reader should never have to infer whether a result is proved, conditional on a frozen-clock comparison, or numerical.

## Negative-rate branch

The manuscript now correctly calls `m'<0` in the ingoing metric a formal negative-rate continuation rather than physical evaporation. Preserve that correction. Nevertheless, figures 7–9 devote a large fraction of the main text to the formal branch and disclaimers. This weakens the positive result. Keep one compact comparison if it supports the asymmetry claim; move the detailed “no inversion” diagnostic and its long caveat to supplementary material.

## Closed-form section

Section 4.4 changes register abruptly and occupies several dense pages. The main paper should retain:

- the first-order expansion and the physical meaning of on-shell versus off-shell terms;
- the universal Vaidya source `S_D=[r p_r]-lambda`;
- one schematic statement of the Abelian/iterated-integral basis;
- the true-flow log–log test showing first-order closure;
- branch and domain qualifications.

Move the full polynomial coefficient system, degeneration-family formulae, transcendence discussion, theta/nome construction and numerical-basis compression to appendices or a technical supplement. Appendix B in particular reads like a separate mathematical note because its double root lies at negative radius and is explicitly not a physical separatrix.

## Reproducibility

The DOI, script map and residual tests are excellent editorial assets. Repair the cramped Table A2 and ensure script names, hashes, parameters, tolerances and the exact archived release agree. A reader should not need the GitHub development branch to reproduce the version in the PDF.
