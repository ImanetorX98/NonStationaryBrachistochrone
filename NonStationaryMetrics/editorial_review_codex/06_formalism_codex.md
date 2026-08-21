# Editorial notes — Formalism

## What works well

- The protocol is now explicit: fixed launch event, fixed spatial target, free arrival clock.
- The manuscript correctly states that the ordinary branch Hamiltonian vanishes only at the free endpoint in the non-autonomous problem, while the extended Hamiltonian is the conserved object.
- The distinction between the conserved axial costate and mechanical angular momentum is made explicitly.
- The control-force identity and “thrust, not free fall” paragraph improve physical honesty.
- The existence/normality assumptions and the HJB/global-optimality qualification are materially better than in earlier drafts.

## Editorial restructuring needed

The general section is too front-loaded. Pages 3–10 ask a broad CQG reader to absorb the indicatrix, four degeneracy loci, the exact endpoint protocol, PMP normality, HJB viscosity solutions, Maxwell sets, an ergosphere disclaimer, the selector hierarchy and Perlick equivalence before seeing the first dynamical result.

Recommended split:

1. Main text: protocol, indicatrix, one boxed theorem statement, one paragraph on endpoint transversality, selector hierarchy and physical control cost.
2. Appendix: direct-method proof, normality proof, HJB verification details and cut-locus discussion.

## Internal claim consistency

The body carefully says generic non-stationary trajectories are PMP extremals with local checks, not globally HJB-certified minimisers. Several figure captions nevertheless say that a curve “minimizes” the functional after testing a one-parameter perturbation. Change those captions to “has a local minimum within the displayed perturbation family” or “passes this variational sanity check.”

Similarly, “well-posed” should always be tied to the compact timelike-selector domain and the stated boundary data. Avoid an unqualified use in the abstract or conclusion.

## Theorem placement and numbering

The central first-order result, Theorem I.5, is buried at the start of Appendix C (page 32). Lemma I.6 appears earlier in Appendix B (page 29), producing the reading order I.4, I.6, I.5. Move the statement of Theorem I.5 into the main Vaidya section or immediately after the general principle, retain its derivation in an appendix, and restore numerical order.

## Paper-I scope leakage

Table 1 includes the Thakurta conformal factor, Kerr spin, conformally normalised charges and radial momentum. Figure 1 gives the companion Thakurta–Kerr indicatrix, and multiple paragraphs revisit ergosphere and conformal-momentum caveats. A short comparison is useful, but the present density makes Paper I feel dependent on Paper II. Keep one forward-looking row/paragraph; remove most Paper-II-specific symbols from the Paper-I notation table.
