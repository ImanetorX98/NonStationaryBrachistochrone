# CQG editorial notes — latest Paper I build

Date: 4 August 2026  
Manuscript reviewed: `paper1/paper1.pdf` (34 pages)  
Source checked: `paper1/paper1.tex`  
Mandate: simulate a CQG Associate Editor; assess editorial readiness rather than certify every mathematical derivation.

## Executive judgement

The manuscript is now clearly within *Classical and Quantum Gravity*'s scope and is strong enough to justify external specialist review. The previous editorial blockers have largely been resolved: the direct relativistic-brachistochrone literature is cited; Theorem I.5 is in the main text; the off-shell response is quantified by the slope-two test; the physical capture threshold and algebraic degeneration are separated; Table A2 is legible; Figure 1 no longer imports the Thakurta–Kerr panel; all 43 references are cited; the conclusions are uninterrupted; and the PDF is 34 pages with zero undefined references, zero duplicate labels and zero overfull boxes.

The current recommendation is nevertheless **minor but mandatory editorial revision before submission**, not immediate upload. The remaining issues are narrow enough to fix without changing the scientific architecture, but several touch central definitions and reproducibility claims. A referee would likely identify them quickly.

## Mandatory corrections

1. **Undefined symbol in the Kodama-invariant derivation.** Section 4.1 begins with
   `-g(T,T)=(f u^v-1)^2/Ehat^2`, but `T` is never defined anywhere in the manuscript. Since the next statement is simply the controlled constraint `-u·K=Ehat`, either define `T` precisely and show the implication, or remove the auxiliary identity if it is unnecessary.

2. **Energy notation changes without definition.** Sections 1–4.3 use `Ehat`; section 4.4 and Appendices A–B switch to `E` in the sextic, `DE`, numerical examples and degeneration formulas. No statement `E≡Ehat` is present. Define the identification once or use one symbol throughout.

3. **Adiabatic parameter is inconsistent.** Table 1 says `epsilon = mdot`, whereas Appendix C and Figure 9 define the dimensionless rate `epsilon = M mdot/m`. The expansion `m(v)=m0(1+epsilon s+...)` also requires an explicit dimensionless clock (for example `s=v/M`) or a declared `M=1` convention. Distinguish the dimensional/coordinate rate `mdot` from the small dimensionless bookkeeping parameter `epsilon` everywhere.

4. **The off-shell source is denoted both `S_D` and `S`.** Equation (22) defines `S_D`, but the paragraph immediately below and the Figure 9 caption call it `S`. This collides with the sextic `S`. Use `S_D` consistently.

5. **Control/evolution parametrisation needs one precise statement.** Section 2.1 calls `s∈{t,tau,eta,v}` the evolution parameter and declares the state to be only `x=(r,phi)`. For a Vaidya proper-time cost, coefficients depend on the advanced time `v`; either `v` must be the path parameter or it must be included as a state. The HJB equation also minimizes over `u`, even though the control was named `theta` and `u` already denotes the four-velocity (and later retarded time). State clearly: path parameter, optimized clock/cost one-form, full state in non-autonomous branches, and control symbol.

6. **A self-containedness claim is presently false.** Appendix C says the genus-one nome series `g^(1)` is “printed in Section 4.4”, but no such series appears in the manuscript. Section 4.4 also says the full Kronecker–Eisenstein evaluation is in Appendix B, while Appendix B gives the block structure but not a complete evaluation prescription. Either print the missing definitions/conventions or narrow the claim to “represented in the text and evaluated by the archived routine”.

7. **The negative-radius degeneration is still called a separatrix in Appendix C.** The manuscript correctly distinguishes `J_pen` from the algebraic `J_deg` at `r_d<0`, but the final Appendix-C paragraph is titled “Separatrix limit” and repeatedly says “separatrix correction/form”. Rename it “algebraic genus-degeneration limit” for Vaidya; reserve separatrix for an accessible dynamical boundary.

8. **The abstract and introduction say “spatial-curvature effects of the spatially flat FLRW base”.** This is semantically misleading. The paper actually contrasts homogeneous expansion/no spatial gradient with Vaidya radial gradients and mass flow. Replace “spatial curvature” with “homogeneous expansion” or “spatial-gradient structure”, as appropriate.

9. **Companion-paper citation does not match its deposited title.** Reference [12] and the cover letter cite “Brachistochrones in conformal Kerr (Thakurta) spacetimes: ergosphere trichotomy, separatrices, and adiabatic response”, while DOI `10.5281/zenodo.21740000` is registered as “Relativistic Brachistochrones in Conformal Kerr Spacetimes: Adiabatic Theory, Ergosphere Transitions and Hyperelliptic Closed Forms”. The citation must match the public record, or the public record must be updated through a new version. Protect `Thakurta` capitalization if that word remains.

10. **The cited preprint/archive snapshots predate the final build.** The cover letter calls DOI `10.5281/zenodo.21739998` the present paper, but that record contains `paper1-10.pdf`, not this 34-page final version. The code DOI `10.5281/zenodo.21707378` is a valid v1.0 release, but the final Figure-1 and Figure-9 production wrappers were created after it. Before submission, create final immutable versions (or remove the outdated “present paper is DOI …” sentence) and make manuscript, cover letter, code snapshot, titles and version numbers mutually consistent.

## Strong recommendations

- Reduce the roughly two dozen visible Paper-II forward references to a small stable set: abstract, introduction, one formalism bridge, conclusion and bibliography. Figure 1's caption should not say that a Thakurta–Kerr image is “retained separately”; a journal reader cannot see that asset.
- Preserve Theorem I.5, equation (22), Appendix C, the turning-point non-uniformity caveat, the conjectural labels and Figure 9. Do not reopen the completed off-shell derivation merely for editorial compression.
- Add an explicit abbreviation introduction for PMP and HJB; spell out BVP before use. Use `theta`, not `u`, for the HJB control.
- Complete the visible bibliography entries for Caponio–Javaloyes–Sánchez and Brown–Levin. The source `.bib` contains their titles and the AMS DOI, but the current IOP `.bst` output suppresses useful metadata. Reference [31] should visibly read “Multiple elliptic polylogarithms, arXiv:1110.6917”.
- Perform a verified DOI pass. IOP style asks for a permanent identifier when available. Do not infer identifiers.
- Align the software citation title with the DataCite/Zenodo record or update the record metadata. The code DOI itself is active, findable, versioned v1.0 and MIT licensed.
- Make captions more self-contained. Figures 4–8 leave essential parameter values mainly in small in-plot titles. Move the necessary parameters to captions and shorten plot titles.
- Export remaining line plots as vector PDF when practical and distinguish every multi-curve result by line style and marker as well as colour. Figures 1 and 9 already meet this standard.
- Remove the last Paper-II implementation details from the general theorem and appendices unless they are necessary to prove a Paper-I statement.
- Consider a final language pass for long sentences and repeated parenthetical caveats in section 4.4 and Appendices B–C. The paper is rigorous, but the density remains the principal desk-review risk for CQG's broad readership.

## Section-by-section status

### Title

Accurate and searchable, although long. Changing it now would require synchronising the Zenodo record and cover letter. Keep it unless the author wants a shorter strategic title.

### Abstract

Under 300 words, self-contained, no citations or equation numbers, and clear about conditional global minimisation. Must correct the “spatial-curvature/FLRW” wording and qualify Kodama with spherical symmetry.

### Introduction

Now cites the direct stationary brachistochrone lineage and states the non-autonomous novelty more defensibly. The conceptual arc is good. Replace “spatial curvature” and keep the companion-paper description to one sentence.

### Section 2

The compact-indicatrix domain, normality argument, HJB qualification and Maxwell/conjugate-point caveats are editorial strengths. The main remaining issue is precise separation of path parameter, clock cost and augmented state. Notation `u` is overloaded.

### Section 3

Effective base case. Figure 2 correctly labels the numerical perturbation as local rather than as a global proof. The analytic assumptions behind `arg min T_t = arg min T_tau = arg min T_eta` should remain explicit.

### Section 4.1

Strong controlled-versus-Noether explanation and useful costate/mechanical-angular-momentum distinction. Fix undefined `T`, define the prime/parameter convention once, and remove unnecessary Paper-II paragraphs.

### Sections 4.2–4.3

Physically readable, with appropriate caveats on the formal negative-rate branch. The scope of numerical claims and parameter choices should be carried by captions rather than plot titles. “Absence of inversion” must always remain qualified as a frozen-clock comparison.

### Section 4.4

Central contribution is visible and Figure 9 is now excellent. The section remains very dense and mixes physical result, special-function taxonomy, algebraic degeneration and companion-paper comparison. Editorial compression is advisable, but the mandatory action is to repair notation and the false `g^(1)`/self-containedness statement.

### Conclusions

The “Proved / Conditional and numerical / Outlook” structure is one of the paper's strongest features. Keep it. No figures now interrupt it. End on Paper I's Vaidya result before the brief companion-paper sentence.

### Appendices

Appendix A is now legible and useful. Appendix B is mathematically substantial but concerns an off-physical-arc degeneration; its physical role must remain modest. Appendix C is central and should stay, but rename the final Vaidya degeneration paragraph and repair the reconstruction claim.

### Back matter

Data statement, funding statement, competing-interests declaration, ORCID and AI disclosure are present. The AI disclosure gives model/version, purposes, verification and author responsibility, consistent in substance with current IOP policy. Submission should use single-anonymous review unless an anonymised manuscript package is deliberately prepared.

## Visual and production audit

- 34 pages; target achieved.
- Zero undefined references/citations.
- Zero duplicate labels.
- Zero overfull boxes.
- One visible References heading.
- PDF title, author, subject and keywords populated.
- Fonts embedded.
- No clipped figures, broken glyphs, overlapping tables or artificial near-empty pages.
- Figure 1 is legible at page width and contains only FLRW/Vaidya.
- Figure 9 is vector, grayscale-safe and quantitatively self-documenting.
- Remaining figures are usable but visually heterogeneous and mostly raster/colour-led.

## Simulated decision

**Decision: invite resubmission after minor mandatory editorial revision, then send to external review.**

This is no longer a “major revision before review” manuscript. Its scientific narrative and claim hierarchy are substantially ready. I would not recommend immediate acceptance, and this editorial audit cannot substitute for specialist verification of the optimal-control proof, extended-Hamiltonian signs, or hyperelliptic reductions. After the ten mandatory consistency/snapshot corrections, however, desk-review risk becomes low to moderate and the paper should be submitted to CQG as a Research Paper.
