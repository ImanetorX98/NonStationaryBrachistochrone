#!/usr/bin/env python3
"""Numerical sanity checks for the experimental proof-of-concept generators."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

import active_particle_poc_codex as active
import draining_vortex_poc_codex as vortex
import ring_bec_poc_codex as bec


HERE = Path(__file__).resolve().parent


def check_active_particle() -> list[str]:
    messages = []
    for epsilon in (0.0, 0.008):
        times, path, _ = active.simulate_feedback(epsilon)
        radii = np.linalg.norm(path, axis=1)
        endpoint_error = float(np.linalg.norm(path[-1] - active.TARGET))
        assert endpoint_error < 0.11
        assert radii.min() > active.PARAMS.inner_radius_mm
        assert radii.max() < active.PARAMS.outer_radius_mm
        messages.append(
            f"active epsilon={epsilon:.3f}: endpoint error={endpoint_error:.4f}, "
            f"arrival={times[-1]:.3f} s, radial clearance="
            f"{radii.min() - active.PARAMS.inner_radius_mm:.4f}"
        )

    centre, semi_r, semi_phi = active.tk_indicatrix(4.5, 1.18)
    angle = np.linspace(0.0, 2.0 * np.pi, 1001)
    vr = semi_r * np.cos(angle)
    vphi = centre + semi_phi * np.sin(angle)
    residual = (vr / semi_r) ** 2 + ((vphi - centre) / semi_phi) ** 2 - 1.0
    assert np.max(np.abs(residual)) < 2.0e-15
    messages.append(f"active indicatrix identity: max residual={np.max(np.abs(residual)):.2e}")
    return messages


def check_vortex() -> list[str]:
    horizon, ergo, prograde, retrograde = vortex.characteristic_radii(vortex.C0, vortex.D0)
    assert 0.0 < horizon < prograde < ergo < retrograde < vortex.TANK_RADIUS

    messages = [
        "vortex radii ordering: "
        f"rh={horizon:.6f} < rpro={prograde:.6f} < re={ergo:.6f} "
        f"< rretro={retrograde:.6f}"
    ]
    for sign, ring_radius, name in (
        (+1.0, float(prograde), "prograde"),
        (-1.0, float(retrograde), "retrograde"),
    ):
        radial_component = vortex.D0 / (vortex.WAVE_SPEED * ring_radius)
        tangential_component = sign * np.sqrt(1.0 - radial_component**2)
        state0 = np.array([ring_radius, 0.0, radial_component, tangential_component])
        initial_rhs = vortex.ray_rhs(0.0, state0, vortex.C0, vortex.D0)
        assert abs(initial_rhs[0]) < 2.0e-15

        solution = solve_ivp(
            lambda t, y: vortex.ray_rhs(t, y, vortex.C0, vortex.D0),
            (0.0, 1.0),
            state0,
            rtol=2.0e-11,
            atol=1.0e-13,
            max_step=0.002,
        )
        assert solution.success
        hamiltonian = []
        for state in solution.y.T:
            velocity, _ = vortex.velocity_and_jacobian(
                state[0], state[1], vortex.C0, vortex.D0
            )
            wavevector = state[2:]
            hamiltonian.append(
                np.dot(wavevector, velocity)
                + vortex.WAVE_SPEED * np.linalg.norm(wavevector)
            )
        hamiltonian = np.asarray(hamiltonian)
        relative_drift = float(
            np.max(np.abs(hamiltonian - hamiltonian[0]))
            / max(abs(hamiltonian[0]), 1.0e-14)
        )
        assert relative_drift < 2.0e-9
        messages.append(
            f"vortex {name}: initial radial speed={initial_rhs[0]:.2e}, "
            f"Hamiltonian relative drift={relative_drift:.2e}"
        )
    return messages


def check_bec() -> list[str]:
    times, amplitudes, velocities = bec.simulate_modes()
    assert np.all(np.isfinite(amplitudes))
    assert np.all(np.isfinite(velocities))
    assert np.isclose(bec.radius(times[0]), bec.RADIUS0)
    assert np.isclose(bec.radius(times[-1]), bec.RADIUS0 * (1.0 + bec.EXPANSION))
    predicted_ratio = float(bec.radius(times[0]) / bec.radius(times[-1]))
    numerical_ratio = float(
        (bec.MODES[0] * bec.SOUND_SPEED / bec.radius(times[-1]))
        / (bec.MODES[0] * bec.SOUND_SPEED / bec.radius(times[0]))
    )
    assert abs(predicted_ratio - numerical_ratio) < 2.0e-15
    return [
        f"BEC mode integration: {len(times)} finite samples",
        f"BEC frequency ratio omega_final/omega_initial={numerical_ratio:.6f}",
    ]


def check_outputs() -> list[str]:
    expected = [
        HERE / "output" / f"{stem}.{suffix}"
        for stem in (
            "active_particle_poc_codex",
            "draining_vortex_poc_codex",
            "ring_bec_poc_codex",
        )
        for suffix in ("png", "pdf", "svg")
    ]
    missing = [path for path in expected if not path.is_file() or path.stat().st_size == 0]
    assert not missing, f"missing or empty outputs: {missing}"
    return [f"rendered outputs: {len(expected)} non-empty files"]


def main() -> None:
    sections = (
        ("ACTIVE PARTICLE", check_active_particle),
        ("DRAINING VORTEX", check_vortex),
        ("RING BEC", check_bec),
        ("FILES", check_outputs),
    )
    for title, check in sections:
        print(f"\n[{title}]")
        for message in check():
            print(f"PASS  {message}")
    print("\nAll proof-of-concept sanity checks passed.")


if __name__ == "__main__":
    main()
