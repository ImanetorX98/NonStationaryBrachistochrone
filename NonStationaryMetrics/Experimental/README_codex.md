# Experimental proof-of-concept generators

This directory contains deterministic, reproducible design studies for three
laboratory validation layers discussed in
`PIANO_VALIDAZIONE_SPERIMENTALE_PAPER_I_II_codex.md`.

They are deliberately separated from the manuscript figures.  Nothing here is
claimed to be measured data, a construction blueprint, or a direct realization
of the Einstein dynamics of Vaidya or Thakurta--Kerr.

## Quick start

From the repository root:

```bash
python3 Experimental/run_all_codex.py
```

Dependencies already used by the project:

- Python 3.10 or newer;
- NumPy;
- SciPy;
- Matplotlib.

Every generator writes PNG, vector PDF and SVG files to
`Experimental/output/`.  That output directory is ignored by Git; the scripts
and this document are not ignored.

The final runner also executes `validate_pocs_codex.py`.  A failed endpoint,
an indicatrix inconsistency, excessive Hamiltonian drift, a non-finite BEC
solution, or a missing render makes the command exit with an error.

## What each script does

### `active_particle_poc_codex.py`

Produces a four-panel design study for the most direct controlled-rail test.

- The velocity oval is not a generic decorative ellipse.  Its centre and
  semi-axes are computed from the equatorial Thakurta--Kerr indicatrix
  
  \[
  \varphi'_0=\frac{2Ma}{r}\frac{\bar v^2}{\bar P},\qquad
  R^2=\frac{\Delta\bar v^2}{\bar P},
  \]
  
  and from the quadratic form in the Paper-II conformal-time Hamiltonian.
- A tracked agent is driven between two fixed spatial endpoints with a
  realizable receding-horizon feedback law while
  \(A(t)=A_0\exp(\varepsilon t)\) is ramped.
- Two intermediate optical reference gates keep this deliberately local
  controller away from the central exclusion disk.  They are controller
  waypoints, not additional endpoints or variational constraints.
- The trajectory is a candidate closed-loop path, not a PMP/HJB certificate of
  global time optimality.
- The log--log panel is explicitly labelled as a pre-registered design target.
  It uses the verified Paper-I orders and slopes, not invented laboratory data:
  \(s_{\rm on}=1.003152\), \(s_{\rm full}=1.999997\).

This is the platform to develop first if the goal is to test the complete
on-shell plus off-shell correction.  The next rigorous upgrade is to replace
the feedback heuristic by the actual branch Hamiltonian boundary-value solver
and then replay that command on hardware.

The labels in millimetres and seconds use the provisional display map
\(L_0=1\,\mathrm{mm}\), \(T_0=1\,\mathrm{s}\).  This makes the geometry
readable but is not yet a hardware calibration; measured propulsion speed,
actuator bandwidth and particle size must set the final \((L_0,T_0)\).

### `draining_vortex_poc_codex.py`

Produces a three-panel design study for a shallow-water draining vortex.

- The frozen model is
  \(u_r=-D/r\), \(u_\varphi=C/r\).
- The ray solver integrates the positive-frequency nondispersive Hamiltonian
  \(\omega=\mathbf{k}\!\cdot\!\mathbf{u}+c|\mathbf{k}|\).
- The analytic reference radii are
  
  \[
  r_h=\frac{D}{c},\qquad
  r_e=\frac{\sqrt{C^2+D^2}}{c},\qquad
  r_{\rm LR}^{\pm}=\frac{\sqrt{2\sqrt{C^2+D^2}
  (\sqrt{C^2+D^2}\mp C)}}{c}.
  \]
- A smooth pump/circulation ramp displays the predicted movement of the two
  light rings, horizon and ergosurface.

The ideal profile must not be fitted directly to a pump command.  A real test
needs PIV inference of \(u_r(r,\varphi,t)\) and
\(u_\varphi(r,\varphi,t)\), followed by ray integration in the measured field.
This platform probes null/acoustic Kerr-like phenomenology and separatrix drift;
it does not directly test a massive controlled rail.

### `ring_bec_poc_codex.py`

Produces a three-panel design study for a slowly breathing toroidal condensate.
The numerical observable is a transparent linear-phonon model,

\[
q_m''+\gamma q_m'+\left(\frac{m c_s}{R(t)}\right)^2q_m=0,
\]

with a smooth programmable expansion of \(R(t)\).  The script generates a
toy-model density-wave kymograph and the corresponding phonon redshift.

This platform establishes a genuinely time-dependent effective background but
remains a wave analogue.  A direct off-shell rail test would require a massive
quasiparticle or an independently controlled material agent.

#### Audit of the present kymograph

The current render is numerically consistent with the equation stated above,
but it is **not yet a quantitative prediction for the measured density of a
breathing ring BEC**.  It currently identifies the oscillator amplitudes with

\[
\delta n_{\rm toy}(\theta,t)
 =\sum_m q_m(t)\cos(m\theta+\phi_m).
\]

That is a useful deterministic illustration of standing-mode interference and
temporal redshift.  A thin-ring BEC comparison must instead evolve the phonon
phase amplitude.  In the notation used by Banik et al.,

\[
\ddot{\delta\phi}_m+
\left[2\gamma(t)+\gamma_H\frac{\dot R}{R}\right]
\dot{\delta\phi}_m+\omega_m^2(t)\delta\phi_m=0,
\]

where

\[
\omega_m(t)=\frac{m c_\theta(t)}{R(t)},\qquad
c_\theta(t)=c_{\theta,i}
\left(\frac{R(t)}{R_i}\right)^{-\alpha/2}.
\]

The ideal thin-ring scaling gives \(\gamma_H=\alpha\), but the experimental
value must be fitted rather than imposed: the measured Hubble attenuation in
the cited experiment does not exactly equal that ideal prediction.  The
observable density amplitude is reconstructed from the phase through

\[
\delta n_m(t)=-\frac{\hbar}{g_{\rm eff}}
R^\alpha(t)\dot{\delta\phi}_m(t),
\qquad
\delta n_{1D}(\theta,t)=
\sum_m\delta n_m(t)\sin(m\theta+\delta\theta_m).
\]

Thus the present script still lacks the Hubble term, the radius dependence of
the sound speed, and the phase-to-density reconstruction.  Until those are
implemented, its panel title should be read as *toy-model phonon kymograph*,
not as synthetic laboratory data.

The current ramp is genuinely slow within its own simplified model.  It has

\[
\max_t\left|\frac{\dot R}{R}\right|=0.0522857,
\qquad
\max_t\frac{|\dot R/R|}{\omega_m}=
0.0100\;(m=6),\quad 0.00667\;(m=9).
\]

These numbers diagnose the present toy parameters only and must be recomputed
after dimensional calibration and after introducing \(c_\theta(t)\).

The use of \(m=6\) and \(m=9\) makes the multimode pattern visually legible but
is not yet experimentally certified.  A controlled first comparison should
start from the nearly pure \(m=1\) excitation used by Banik et al.  Higher modes
are admissible only after checking the hydrodynamic condition
\(k_m\xi=m\xi/R\ll1\), separation from radial modes, and linear response.

Finally, the horizontal coordinate is the comoving angle \(\theta\).  Its mode
number \(m\) therefore remains fixed during expansion: the redshift appears as
a reduced temporal frequency, not as horizontal stretching.  The physical
wavelength \(2\pi R(t)/m\) does grow with the ring.

Primary references: S. Eckel et al.,
[*Phys. Rev. X* **8**, 021021 (2018)](https://doi.org/10.1103/PhysRevX.8.021021),
and S. Banik et al.,
[*Phys. Rev. Lett.* **128**, 090401 (2022)](https://doi.org/10.1103/PhysRevLett.128.090401).

The required implementation upgrade is therefore:

1. choose measured or preregistered \(R(t)\), \(\alpha\), \(c_{\theta,i}\),
   \(Q(R)\), and \(\gamma_H\);
2. integrate \(\delta\phi_m\), not a variable labelled directly as density;
3. reconstruct both absolute \(\delta n_{1D}\) and normalized
   \(\delta n_{1D}/n_{0,1D}\);
4. begin with \(m=1\), then add only hydrodynamically validated higher modes;
5. compare synthetic time slices with the same radial integration and imaging
   pipeline used for laboratory data.

### `experimental_style_codex.py` and `run_all_codex.py`

The first keeps paper-compatible styles and saves three formats.  The second is
the reproducible entry point for all generators.

### `validate_pocs_codex.py`

Checks fixed-endpoint arrival and annular clearance for both the frozen and
ramped rails, the exact ellipse identity, the DBT radius ordering and
stationary Hamiltonian conservation, the BEC frequency ratio, and all nine
rendered files.  These are numerical sanity checks, not experimental evidence.

## Scientific status and next upgrades

The present scripts establish realistic geometry, observables and falsifiable
calibration steps.  They do **not** yet establish that a laboratory field is the
inverse image of the full Paper-I/Paper-II Hamiltonian.

Recommended next steps, in order:

1. define the dimensional map \((L_0,T_0)\) and actuator limits;
2. replace ideal fields by measured transfer functions, latency and noise;
3. solve the fixed-endpoint PMP/HJB problem with the actual Hamiltonian;
4. export the resulting command table for hardware-in-the-loop replay;
5. preregister the residual norm, the \(\varepsilon\) interval and the slope fit;
6. only then design the physical bill of materials and uncertainty budget.
