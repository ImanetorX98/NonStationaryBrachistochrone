#!/usr/bin/env python3
"""
Core numerical routines for TOV interior brachistochrones.

The module works in geometrized units G = c = 1.  The default equation of
state is a relativistic polytrope,

    p = K rho^Gamma,    eps = rho + p / (Gamma - 1),

where rho is the rest-mass density and eps is the total energy density.
"""

from __future__ import annotations

import math
import csv
import os
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np


GAUSS_ORDER = 220
GAUSS_X, GAUSS_W = np.polynomial.legendre.leggauss(GAUSS_ORDER)
MEV_FM3_TO_KM2 = 1.323409798288609e-6


@dataclass(frozen=True)
class PolytropicEOS:
    kappa: float = 100.0
    gamma: float = 2.0
    name: str = "relativistic polytrope"

    def rho_from_pressure(self, pressure: float | np.ndarray) -> float | np.ndarray:
        pressure_arr = np.maximum(pressure, 0.0)
        return (pressure_arr / self.kappa) ** (1.0 / self.gamma)

    def pressure_from_density(self, density: float | np.ndarray) -> float | np.ndarray:
        density_arr = np.maximum(density, 0.0)
        return self.kappa * density_arr**self.gamma

    def epsilon_from_pressure(self, pressure: float | np.ndarray) -> float | np.ndarray:
        rho = self.rho_from_pressure(pressure)
        return rho + np.maximum(pressure, 0.0) / (self.gamma - 1.0)


def _normalize_column_name(name: str) -> str:
    return name.strip().lower().replace(" ", "_").replace("-", "_")


def _read_table_rows(path: str) -> Tuple[List[str], List[List[str]]]:
    with open(path, "r", encoding="utf-8") as handle:
        lines = [line.strip() for line in handle if line.strip() and not line.lstrip().startswith("#")]
    if not lines:
        raise ValueError(f"EOS table is empty: {path}")

    delimiter = "," if "," in lines[0] else None
    if delimiter == ",":
        parsed = list(csv.reader(lines))
    else:
        parsed = [line.split() for line in lines]
    header = [_normalize_column_name(item) for item in parsed[0]]
    rows = parsed[1:]
    return header, rows


def _find_column(header: Sequence[str], requested: Optional[str], candidates: Sequence[str]) -> Optional[int]:
    if requested:
        requested_key = _normalize_column_name(requested)
        if requested_key not in header:
            raise ValueError(f"Column {requested!r} not found in EOS table. Available columns: {', '.join(header)}")
        return list(header).index(requested_key)
    for name in candidates:
        key = _normalize_column_name(name)
        if key in header:
            return list(header).index(key)
    return None


def _as_positive_array(rows: Sequence[Sequence[str]], column: int, label: str) -> np.ndarray:
    values = np.array([float(row[column]) for row in rows], dtype=float)
    if np.any(~np.isfinite(values)):
        raise ValueError(f"Column {label} contains non-finite values.")
    if np.any(values < 0.0):
        raise ValueError(f"Column {label} contains negative values.")
    return values


def _log_interp_extrap(x: np.ndarray | float, xp: np.ndarray, fp: np.ndarray) -> np.ndarray | float:
    x_arr = np.asarray(x, dtype=float)
    result = np.zeros_like(x_arr, dtype=float)
    positive = x_arr > 0.0
    if np.any(positive):
        log_x = np.log(x_arr[positive])
        log_xp = np.log(xp)
        log_fp = np.log(fp)
        log_y = np.interp(log_x, log_xp, log_fp)
        low = log_x < log_xp[0]
        high = log_x > log_xp[-1]
        if np.any(low):
            slope = (log_fp[1] - log_fp[0]) / (log_xp[1] - log_xp[0])
            log_y[low] = log_fp[0] + slope * (log_x[low] - log_xp[0])
        if np.any(high):
            slope = (log_fp[-1] - log_fp[-2]) / (log_xp[-1] - log_xp[-2])
            log_y[high] = log_fp[-1] + slope * (log_x[high] - log_xp[-1])
        result[positive] = np.exp(log_y)
    if np.isscalar(x):
        return float(result)
    return result


@dataclass(frozen=True)
class TabulatedEOS:
    density: np.ndarray
    pressure: np.ndarray
    epsilon: np.ndarray
    name: str = "tabulated EOS"

    @classmethod
    def from_csv(
        cls,
        path: str,
        density_column: Optional[str] = None,
        pressure_column: Optional[str] = None,
        epsilon_column: Optional[str] = None,
        name: Optional[str] = None,
    ) -> "TabulatedEOS":
        header, rows = _read_table_rows(path)
        density_index = _find_column(
            header,
            density_column,
            ("density", "rho", "rho0", "rest_mass_density", "energy_density", "epsilon", "eps"),
        )
        pressure_index = _find_column(header, pressure_column, ("pressure", "p", "press"))
        epsilon_index = _find_column(
            header,
            epsilon_column,
            ("epsilon", "eps", "energy_density", "total_energy_density", "e"),
        )
        if density_index is None or pressure_index is None:
            raise ValueError(
                "EOS table must contain density/rho and pressure/p columns, "
                f"or pass explicit column names. Available columns: {', '.join(header)}"
            )

        density = _as_positive_array(rows, density_index, "density")
        pressure = _as_positive_array(rows, pressure_index, "pressure")
        epsilon = (
            _as_positive_array(rows, epsilon_index, "epsilon")
            if epsilon_index is not None
            else density.copy()
        )

        positive = (density > 0.0) & (pressure > 0.0) & (epsilon > 0.0)
        density = density[positive]
        pressure = pressure[positive]
        epsilon = epsilon[positive]
        if len(density) < 4:
            raise ValueError("EOS table needs at least four positive rows.")

        order = np.argsort(density)
        density = density[order]
        pressure = pressure[order]
        epsilon = epsilon[order]
        keep_density = np.concatenate([[True], np.diff(density) > 0.0])
        density = density[keep_density]
        pressure = pressure[keep_density]
        epsilon = epsilon[keep_density]

        pressure_order = np.argsort(pressure)
        pressure = pressure[pressure_order]
        density_for_pressure = density[pressure_order]
        epsilon_for_pressure = epsilon[pressure_order]
        keep_pressure = np.concatenate([[True], np.diff(pressure) > 0.0])

        table_name = name or os.path.splitext(os.path.basename(path))[0]
        return cls(
            density=density_for_pressure[keep_pressure],
            pressure=pressure[keep_pressure],
            epsilon=epsilon_for_pressure[keep_pressure],
            name=f"tabulated:{table_name}",
        )

    @classmethod
    def from_compose_ns(
        cls,
        nb_path: str,
        thermo_path: str,
        name: Optional[str] = None,
    ) -> "TabulatedEOS":
        nb, first_index = read_compose_nb(nb_path)
        mn_mev, thermo_rows = read_compose_thermo(thermo_path)
        density = []
        pressure = []
        epsilon = []
        for row in thermo_rows:
            inb = int(round(row[1]))
            nb_index = inb - first_index
            if nb_index < 0 or nb_index >= len(nb):
                continue
            baryon_density = nb[nb_index]
            q1 = row[3]
            q7 = row[9]
            density.append(baryon_density)
            pressure.append(baryon_density * q1 * MEV_FM3_TO_KM2)
            epsilon.append(baryon_density * mn_mev * (1.0 + q7) * MEV_FM3_TO_KM2)

        density_arr = np.asarray(density, dtype=float)
        pressure_arr = np.asarray(pressure, dtype=float)
        epsilon_arr = np.asarray(epsilon, dtype=float)
        positive = (density_arr > 0.0) & (pressure_arr > 0.0) & (epsilon_arr > 0.0)
        density_arr = density_arr[positive]
        pressure_arr = pressure_arr[positive]
        epsilon_arr = epsilon_arr[positive]
        if len(density_arr) < 4:
            raise ValueError("CompOSE ns files did not contain enough positive thermodynamic rows.")

        order = np.argsort(density_arr)
        table_name = name or os.path.basename(os.path.dirname(os.path.abspath(nb_path))) or "compose_ns"
        return cls(
            density=density_arr[order],
            pressure=pressure_arr[order],
            epsilon=epsilon_arr[order],
            name=f"CompOSE:{table_name}",
        )

    def pressure_from_density(self, density: float | np.ndarray) -> float | np.ndarray:
        return _log_interp_extrap(density, self.density, self.pressure)

    def rho_from_pressure(self, pressure: float | np.ndarray) -> float | np.ndarray:
        return _log_interp_extrap(pressure, self.pressure, self.density)

    def epsilon_from_pressure(self, pressure: float | np.ndarray) -> float | np.ndarray:
        return _log_interp_extrap(pressure, self.pressure, self.epsilon)


def read_compose_nb(path: str) -> Tuple[np.ndarray, int]:
    with open(path, "r", encoding="utf-8") as handle:
        values = [float(token) for token in handle.read().split()]
    if len(values) < 4:
        raise ValueError(f"CompOSE density grid is too short: {path}")
    first_index = int(round(values[0]))
    last_index = int(round(values[1]))
    expected = last_index - first_index + 1
    if expected > 0 and len(values) >= expected + 2:
        return np.asarray(values[2 : 2 + expected], dtype=float), first_index
    return np.asarray(values, dtype=float), 0


def read_compose_thermo(path: str) -> Tuple[float, List[List[float]]]:
    with open(path, "r", encoding="utf-8") as handle:
        lines = [line.strip() for line in handle if line.strip() and not line.lstrip().startswith("#")]
    if len(lines) < 2:
        raise ValueError(f"CompOSE thermodynamic table is too short: {path}")
    header = [float(token) for token in lines[0].split()]
    mn_mev = float(header[0])
    rows = []
    for line in lines[1:]:
        values = [float(token) for token in line.split()]
        if len(values) >= 11:
            rows.append(values)
    return mn_mev, rows


@dataclass
class TOVSolution:
    r: np.ndarray
    pressure: np.ndarray
    epsilon: np.ndarray
    mass: np.ndarray
    phi_raw: np.ndarray
    lapse: np.ndarray
    radial_metric: np.ndarray
    radius: float
    total_mass: float
    compactness: float
    surface_lapse: float
    central_pressure: float
    eos: PolytropicEOS | TabulatedEOS

    def interp(self, values: np.ndarray, radii: np.ndarray | float) -> np.ndarray | float:
        return np.interp(radii, self.r, values)

    def sample(self, n: int = 1200, stop: float = 1.0) -> np.ndarray:
        stop = min(max(stop, 0.0), 1.0)
        return np.linspace(0.0, stop * self.radius, n)


def tov_rhs(radius: float, state: np.ndarray, eos: PolytropicEOS | TabulatedEOS) -> np.ndarray:
    mass, pressure, phi = state
    if radius <= 0.0 or pressure <= 0.0:
        return np.zeros(3)
    epsilon = float(eos.epsilon_from_pressure(pressure))
    metric_denom = radius * (radius - 2.0 * mass)
    if metric_denom <= 0.0:
        raise RuntimeError("TOV integration reached r <= 2m; compactness is too high.")
    common = mass + 4.0 * math.pi * radius**3 * pressure
    dmass = 4.0 * math.pi * radius**2 * epsilon
    dpressure = -(epsilon + pressure) * common / metric_denom
    dphi = common / metric_denom
    return np.array([dmass, dpressure, dphi], dtype=float)


def rk4_step(radius: float, state: np.ndarray, step: float, eos: PolytropicEOS | TabulatedEOS) -> np.ndarray:
    k1 = tov_rhs(radius, state, eos)
    k2 = tov_rhs(radius + 0.5 * step, state + 0.5 * step * k1, eos)
    k3 = tov_rhs(radius + 0.5 * step, state + 0.5 * step * k2, eos)
    k4 = tov_rhs(radius + step, state + step * k3, eos)
    return state + step * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0


def integrate_tov(
    central_pressure: float,
    eos: PolytropicEOS | TabulatedEOS | None = None,
    dr: float = 2.0e-3,
    max_radius: float = 80.0,
    pressure_floor_fraction: float = 1.0e-10,
) -> TOVSolution:
    if eos is None:
        eos = PolytropicEOS()
    if central_pressure <= 0.0:
        raise ValueError("central_pressure must be positive.")
    if dr <= 0.0:
        raise ValueError("dr must be positive.")

    pressure_floor = pressure_floor_fraction * central_pressure
    epsilon_c = float(eos.epsilon_from_pressure(central_pressure))
    start_radius = min(dr, 1.0e-4)
    start_mass = 4.0 * math.pi * epsilon_c * start_radius**3 / 3.0
    state = np.array([start_mass, central_pressure, 0.0], dtype=float)

    radii: List[float] = [0.0, start_radius]
    masses: List[float] = [0.0, start_mass]
    pressures: List[float] = [central_pressure, central_pressure]
    phis: List[float] = [0.0, 0.0]

    radius = start_radius
    while radius < max_radius and state[1] > pressure_floor:
        previous_radius = radius
        previous_state = state.copy()
        next_state = rk4_step(radius, state, dr, eos)
        next_radius = radius + dr

        if next_state[1] <= pressure_floor:
            fraction = (pressure_floor - previous_state[1]) / (next_state[1] - previous_state[1])
            fraction = min(max(float(fraction), 0.0), 1.0)
            surface_radius = previous_radius + fraction * dr
            surface_state = previous_state + fraction * (next_state - previous_state)
            radii.append(surface_radius)
            masses.append(float(surface_state[0]))
            pressures.append(0.0)
            phis.append(float(surface_state[2]))
            break

        radius = next_radius
        state = next_state
        radii.append(radius)
        masses.append(float(state[0]))
        pressures.append(float(state[1]))
        phis.append(float(state[2]))
    else:
        raise RuntimeError("TOV surface was not reached before max_radius.")

    r = np.array(radii, dtype=float)
    pressure = np.array(pressures, dtype=float)
    mass = np.array(masses, dtype=float)
    phi_raw = np.array(phis, dtype=float)
    epsilon = np.asarray(eos.epsilon_from_pressure(pressure), dtype=float)

    radius = float(r[-1])
    total_mass = float(mass[-1])
    compactness = total_mass / radius
    surface_lapse = 1.0 - 2.0 * compactness
    if surface_lapse <= 0.0:
        raise RuntimeError("Surface lies inside its Schwarzschild radius.")

    phi_shift = 0.5 * math.log(surface_lapse) - float(phi_raw[-1])
    lapse = np.exp(2.0 * (phi_raw + phi_shift))
    radial_metric = np.ones_like(r)
    mask = r > 0.0
    radial_metric[mask] = 1.0 / np.maximum(1.0 - 2.0 * mass[mask] / r[mask], 1.0e-14)

    return TOVSolution(
        r=r,
        pressure=pressure,
        epsilon=epsilon,
        mass=mass,
        phi_raw=phi_raw,
        lapse=lapse,
        radial_metric=radial_metric,
        radius=radius,
        total_mass=total_mass,
        compactness=compactness,
        surface_lapse=surface_lapse,
        central_pressure=central_pressure,
        eos=eos,
    )


def local_speed(sol: TOVSolution, radii: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    lapse = np.asarray(sol.interp(sol.lapse, radii), dtype=float)
    gamma = np.sqrt(sol.surface_lapse / lapse)
    speed = np.sqrt(np.maximum(1.0 - lapse / sol.surface_lapse, 0.0))
    return speed, gamma


def effective_index(sol: TOVSolution, radii: np.ndarray | float, mode: str) -> np.ndarray | float:
    lapse = np.asarray(sol.interp(sol.lapse, radii), dtype=float)
    gap = np.maximum(sol.surface_lapse - lapse, 1.0e-15)
    if mode == "coordinate":
        value = 1.0 / np.sqrt(np.maximum(lapse * gap, 1.0e-300))
    elif mode == "proper":
        value = np.sqrt(lapse / gap)
    else:
        raise ValueError("mode must be 'coordinate' or 'proper'.")
    if np.isscalar(radii):
        return float(value)
    return value


def _integral_nodes() -> Tuple[np.ndarray, np.ndarray]:
    nodes = 0.25 * math.pi * (GAUSS_X + 1.0)
    weights = 0.25 * math.pi * GAUSS_W
    return nodes, weights


def _phi_integrand_u(sol: TOVSolution, mode: str, rstar: float, u: np.ndarray) -> np.ndarray:
    radius = rstar + (sol.radius - rstar) * np.sin(u) ** 2
    dr_du = (sol.radius - rstar) * np.sin(2.0 * u)
    nstar = float(effective_index(sol, rstar, mode))
    angular_momentum = nstar * rstar
    n_values = np.asarray(effective_index(sol, radius, mode), dtype=float)
    b_values = np.asarray(sol.interp(sol.radial_metric, radius), dtype=float)
    denominator_sq = np.maximum((n_values * radius) ** 2 - angular_momentum**2, 1.0e-300)
    with np.errstate(divide="ignore", invalid="ignore"):
        dphi_dr = angular_momentum * np.sqrt(b_values) / (radius * np.sqrt(denominator_sq))
    return np.nan_to_num(dphi_dr * dr_du, nan=0.0, posinf=0.0, neginf=0.0)


def half_opening_angle(sol: TOVSolution, mode: str, rstar: float) -> float:
    if not (0.0 < rstar < sol.radius):
        raise ValueError("rstar must lie inside the star.")
    nodes, weights = _integral_nodes()
    return float(np.sum(weights * _phi_integrand_u(sol, mode, rstar, nodes)))


def solve_turning_radius(
    sol: TOVSolution,
    mode: str,
    delta: float,
    grid_size: int = 260,
) -> float:
    if not (0.0 < delta < math.pi):
        raise ValueError("delta must satisfy 0 < delta < pi for this solver.")
    target = 0.5 * delta
    q_grid = np.linspace(1.0e-5, 0.999999, grid_size)
    r_grid = q_grid * sol.radius
    values = np.array([half_opening_angle(sol, mode, float(r)) for r in r_grid])
    residuals = values - target

    candidates: List[Tuple[float, float]] = []
    for i in range(len(r_grid) - 1):
        if not (math.isfinite(residuals[i]) and math.isfinite(residuals[i + 1])):
            continue
        if residuals[i] == 0.0:
            return float(r_grid[i])
        if residuals[i] * residuals[i + 1] <= 0.0:
            candidates.append((float(r_grid[i]), float(r_grid[i + 1])))

    if not candidates:
        max_delta = 2.0 * float(np.nanmax(values))
        raise RuntimeError(
            f"No turning radius for delta={math.degrees(delta):.3f} deg in mode={mode}; "
            f"largest sampled aperture is {math.degrees(max_delta):.3f} deg."
        )

    low, high = candidates[-1]
    f_low = half_opening_angle(sol, mode, low) - target
    for _ in range(80):
        mid = 0.5 * (low + high)
        f_mid = half_opening_angle(sol, mode, mid) - target
        if abs(f_mid) < 1.0e-11 or abs(high - low) < 1.0e-11 * sol.radius:
            return mid
        if f_low * f_mid <= 0.0:
            high = mid
        else:
            low = mid
            f_low = f_mid
    return 0.5 * (low + high)


def travel_time_over_radius(sol: TOVSolution, mode: str, rstar: float) -> float:
    if not (0.0 < rstar < sol.radius):
        raise ValueError("rstar must lie inside the star.")
    nodes, weights = _integral_nodes()
    radius = rstar + (sol.radius - rstar) * np.sin(nodes) ** 2
    dr_du = (sol.radius - rstar) * np.sin(2.0 * nodes)
    nstar = float(effective_index(sol, rstar, mode))
    angular_momentum = nstar * rstar
    n_values = np.asarray(effective_index(sol, radius, mode), dtype=float)
    b_values = np.asarray(sol.interp(sol.radial_metric, radius), dtype=float)
    denominator_sq = np.maximum((n_values * radius) ** 2 - angular_momentum**2, 1.0e-300)
    integrand = n_values**2 * radius * np.sqrt(b_values) / np.sqrt(denominator_sq)
    integral = float(np.sum(weights * np.nan_to_num(integrand * dr_du, nan=0.0, posinf=0.0)))
    if mode == "coordinate":
        return 2.0 * math.sqrt(sol.surface_lapse) * integral / sol.radius
    if mode == "proper":
        return 2.0 * integral / sol.radius
    raise ValueError("mode must be 'coordinate' or 'proper'.")


def curve_branch(
    sol: TOVSolution,
    mode: str,
    rstar: float,
    n_points: int = 1800,
) -> Tuple[np.ndarray, np.ndarray]:
    u = np.linspace(0.0, 0.5 * math.pi, n_points)
    safe_u = u.copy()
    if n_points < 4:
        raise ValueError("n_points must be at least 4.")
    safe_u[0] = safe_u[1]
    safe_u[-1] = safe_u[-2]
    radius = rstar + (sol.radius - rstar) * np.sin(u) ** 2
    integrand = _phi_integrand_u(sol, mode, rstar, safe_u)
    integrand[0] = integrand[1]
    integrand[-1] = 0.0
    phi = np.zeros_like(u)
    phi[1:] = np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(u))
    return radius, phi


def full_curve_xy(
    sol: TOVSolution,
    mode: str,
    delta: float,
    n_points: int = 1800,
) -> Tuple[np.ndarray, np.ndarray, float, float]:
    rstar = solve_turning_radius(sol, mode, delta)
    radius, phi = curve_branch(sol, mode, rstar, n_points)
    left_r = radius[::-1]
    left_phi = -phi[::-1]
    right_r = radius[1:]
    right_phi = phi[1:]
    full_r = np.concatenate([left_r, right_r]) / sol.radius
    full_phi = np.concatenate([left_phi, right_phi])
    x = full_r * np.cos(full_phi)
    y = full_r * np.sin(full_phi)
    return x, y, rstar / sol.radius, 2.0 * float(phi[-1])


def aperture_scan(
    sol: TOVSolution,
    delta_degrees: Sequence[float],
) -> Dict[str, np.ndarray]:
    delta_rad = np.deg2rad(np.asarray(delta_degrees, dtype=float))
    q_coordinate = np.full_like(delta_rad, np.nan)
    q_proper = np.full_like(delta_rad, np.nan)
    time_coordinate = np.full_like(delta_rad, np.nan)
    time_proper = np.full_like(delta_rad, np.nan)

    for i, delta in enumerate(delta_rad):
        for mode, q_arr, time_arr in (
            ("coordinate", q_coordinate, time_coordinate),
            ("proper", q_proper, time_proper),
        ):
            try:
                rstar = solve_turning_radius(sol, mode, float(delta))
                q_arr[i] = rstar / sol.radius
                time_arr[i] = travel_time_over_radius(sol, mode, rstar)
            except Exception:
                q_arr[i] = np.nan
                time_arr[i] = np.nan

    return {
        "delta_deg": np.asarray(delta_degrees, dtype=float),
        "q_coordinate": q_coordinate,
        "q_proper": q_proper,
        "time_coordinate": time_coordinate,
        "time_proper": time_proper,
    }


def compactness_sequence(
    central_pressures: Iterable[float],
    eos: PolytropicEOS | TabulatedEOS | None = None,
    dr: float = 2.0e-3,
) -> List[TOVSolution]:
    if eos is None:
        eos = PolytropicEOS()
    solutions = []
    for pressure in central_pressures:
        solutions.append(integrate_tov(float(pressure), eos=eos, dr=dr))
    return solutions
