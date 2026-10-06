# Paper III — formal verification in progress

This directory preserves the formal proof developments for Paper III. It is
separate from the sources of Papers I and II.

## Vaidya: exclusion of Maxwell points before π (via B)

The Lean project is in
[`vaidya_no_maxwell_before_pi/`](vaidya_no_maxwell_before_pi/README.md).
Its name identifies the target theorem, **not a completed certification of that
theorem**. The current checkpoint contains **142 verified public lemmas** and
34 Lean source files, including the axiom audit.

Certified components include the scalar Sturm comparison, the declared ODE
escape/lingering lemmas, and finite-time dependence and escape stability with
the uniform Lipschitz bound derived from a compact regular reference. The
coordinate metric inverse, direction ODE with a constructed local time inverse,
propagation of the null constraint and angular-momentum conservation for the
declared nonstationary Hamiltonian flow are now certified. Local existence
under C¹ mass, local uniqueness under continuous mass, and continuation with
gluing for trajectories confined to regular compact state–mass tubes are also
verified. Geometric identification, maximal/global exterior evolution,
properness, covering/fiber count and the final Maxwell bridge remain to be
formalized.

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
