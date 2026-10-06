# Paper III — formal verification in progress

This directory preserves the formal proof developments for Paper III. It is
separate from the sources of Papers I and II.

## Vaidya: exclusion of Maxwell points before π (via B)

The Lean project is in
[`vaidya_no_maxwell_before_pi/`](vaidya_no_maxwell_before_pi/README.md).
Its name identifies the target theorem, **not a completed certification of that
theorem**. The current checkpoint contains **248 verified public lemmas** and
51 Lean source files, including the axiom audit.

Certified components include the scalar Sturm comparison, the declared ODE
escape/lingering lemmas, and finite-time dependence and escape stability with
the uniform Lipschitz bound derived from a compact regular reference. The
coordinate metric inverse, direction ODE with a constructed local time inverse,
propagation of the null constraint and angular-momentum conservation for the
declared nonstationary Hamiltonian flow are now certified. Local existence
under C¹ mass, local uniqueness under continuous mass, and continuation with
gluing for trajectories confined to regular compact state–mass tubes are also
verified. For positive nondecreasing mass, finite exterior evolution now
derives compact confinement and the nonradial angular strip from the launch,
and gives capture or longer exterior evolution at a finite endpoint.
The maximal exterior coordinate solution is now constructed by gluing all
compatible local extensions. A finite maximal endpoint must be capture; strict
outgoing-cone entry instead produces a future-global DirectionFlow and radial
escape, without assuming global existence. The cone launch theorem also
constructs its initial local segment and cone margins internally.
The converse coordinate bridge is now verified: the direction solution
constructs its future-null momenta, angular primitive, affine clock and inverse,
and satisfies all affine Hamiltonian equations on a neighborhood. This also
applies to the constructed maximal exterior curve. The mass-rate time force is
derived rather than assumed; future affine completeness is not asserted.
The coordinate Levi-Civita connection is now derived from actual metric
derivatives and proved torsion-free, metric-compatible and unique. Independently
derived Hamiltonian acceleration agrees with its geodesic equation; actual
second derivatives, nullness and future orientation hold on a constructed affine
neighborhood, also for the maximal exterior curve. Differentiation of the
geodesic equation is now verified under explicit family-derivative hypotheses.
Connection differentiability is derived, with C² physical mass a sufficient
condition along differentiable sampled coordinates. Angular translations give
a concrete verified family, but do not supply the launch-vanishing screen field.
Mixed derivatives and the actual affine derivatives of the launch-variation
field are now certified for supplied jointly C³ geodesic families, including
Fermat null-geodesic families. Launch from a fixed event forces the variation
field to vanish there. Joint C³ regularity of the constructed launch family
is still unproved and is not implied here by C² mass. Launch-family regularity,
the geometric screen/curvature bridge,
properness, covering/fiber count and
the final Maxwell bridge remain to be formalized. These results do not assert
escape for every exterior launch.

Lean **4.24.0** and mathlib **v4.24.0** are pinned by the toolchain and dependency
manifest. From the project directory, after installing the pinned toolchain:

```sh
lake exe cache get
python3 verify.py
```

The script compiles the project and audits every public theorem, allowing only
Lean/mathlib's standard axioms. The `verification/` directory records the build,
axiom audit and source hashes. Downloaded dependencies and build artifacts in
`.lake/` are regenerated and excluded from version control.

Continue development in this directory. The earlier copy under the ignored
CQG audit directory is a historical checkpoint, not the canonical source.
