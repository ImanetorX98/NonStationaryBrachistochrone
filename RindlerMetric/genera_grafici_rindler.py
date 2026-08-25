#!/usr/bin/env python3
"""
Generate publication figures for the Rindler brachistochrone paper.

The script uses only numpy plus optional reportlab for PDF export.  SVG and CSV
files are always produced; PDF files are produced when reportlab is available.
Most plotted lengths are dimensionless, using a reference Newtonian maximum
depth H0 as the length scale.
"""

import argparse
import csv
import math
import os
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

try:
    from reportlab.lib import colors
    from reportlab.pdfgen import canvas

    HAS_REPORTLAB = True
except Exception:
    HAS_REPORTLAB = False

try:
    from PIL import Image, ImageDraw, ImageFont

    HAS_PIL = True
except Exception:
    HAS_PIL = False


Color = str
Point = Tuple[float, float]
GAUSS_ORDER = 160
GAUSS_NODES, GAUSS_WEIGHTS = np.polynomial.legendre.leggauss(GAUSS_ORDER)
PNG_DPI: int = 300


PALETTE: Dict[str, Color] = {
    "black": "#111111",
    "gray": "#666666",
    "light_gray": "#d8d8d8",
    "newton": "#4d4d4d",
    "coord": "#0072B2",
    "proper": "#D55E00",
    "nlo": "#009E73",
    "nnlo": "#CC79A7",
    "rho1": "#0072B2",
    "rho2": "#009E73",
    "rho3": "#D55E00",
    "rho4": "#CC79A7",
}


@dataclass
class Series:
    x: np.ndarray
    y: np.ndarray
    label: str
    color: Color = "#111111"
    width: float = 2.0
    style: str = "solid"


@dataclass
class Panel:
    rect: Tuple[float, float, float, float]
    xlim: Tuple[float, float]
    ylim: Tuple[float, float]
    xlabel: str
    ylabel: str
    title: str
    series: List[Series] = field(default_factory=list)
    legend: bool = True
    equal_aspect: bool = False
    xticks: Optional[Sequence[float]] = None
    yticks: Optional[Sequence[float]] = None


def parse_rho_list(text: str) -> List[float]:
    values = []
    for raw in text.split(","):
        item = raw.strip()
        if item:
            values.append(float(item))
    if not values:
        raise ValueError("rho list is empty")
    return values


def rho_tag(rho: float) -> str:
    return f"{rho:.3f}".rstrip("0").rstrip(".").replace(".", "p")


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def cumulative_trapezoid(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    out = np.zeros_like(x, dtype=float)
    out[1:] = np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(x))
    return out


def elliptic_f_minus_e(phi: np.ndarray, m: float) -> np.ndarray:
    sin2 = np.sin(phi) ** 2
    root = np.sqrt(np.maximum(1.0 - m * sin2, 1.0e-15))
    integrand = 1.0 / root - root
    return cumulative_trapezoid(integrand, phi)


def elliptic_f_minus_e_scalar(phi: float, m: float) -> float:
    if phi <= 0.0:
        return 0.0
    nodes = 0.5 * phi * (GAUSS_NODES + 1.0)
    weights = 0.5 * phi * GAUSS_WEIGHTS
    root_arg = 1.0 - m * np.sin(nodes) ** 2
    if np.any(root_arg <= 0.0):
        return float("nan")
    root = np.sqrt(root_arg)
    return float(np.sum(weights * (1.0 / root - root)))


def exact_curve(rho: float, mode: str, n: int = 3000) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    if not (0.0 < rho < 1.0):
        raise ValueError("rho must satisfy 0 < rho < 1")
    eta = 1.0 - rho
    one_minus_eta2 = 1.0 - eta * eta
    if mode == "coordinate":
        m = one_minus_eta2 / (eta * eta)
    elif mode == "proper":
        m = one_minus_eta2
    else:
        raise ValueError("mode must be 'coordinate' or 'proper'")

    phi = np.linspace(0.0, math.pi, n)
    theta = 2.0 * phi
    f_minus_e = elliptic_f_minus_e(phi, m)
    x = eta / rho * f_minus_e
    y = -(1.0 - np.sqrt(np.maximum(1.0 - one_minus_eta2 * np.sin(phi) ** 2, 0.0))) / rho
    return theta, x, y


def exact_curve_fixed_scale(
    rho: float,
    mode: str,
    rho_ref: float,
    phi_end: float,
    n: int = 3000,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    if not (0.0 < rho < 1.0):
        raise ValueError("rho must satisfy 0 < rho < 1")
    eta = 1.0 - rho
    one_minus_eta2 = 1.0 - eta * eta
    if mode == "coordinate":
        m = one_minus_eta2 / (eta * eta)
    elif mode == "proper":
        m = one_minus_eta2
    else:
        raise ValueError("mode must be 'coordinate' or 'proper'")

    phi = np.linspace(0.0, phi_end, n)
    f_minus_e = elliptic_f_minus_e(phi, m)
    exact_endpoint = elliptic_f_minus_e_scalar(phi_end, m)
    if math.isfinite(exact_endpoint) and abs(f_minus_e[-1]) > 1.0e-15:
        f_minus_e *= exact_endpoint / f_minus_e[-1]
    x = eta / rho_ref * f_minus_e
    y = -(1.0 - np.sqrt(np.maximum(1.0 - one_minus_eta2 * np.sin(phi) ** 2, 0.0))) / rho_ref
    return phi, x, y


def newtonian_curve(theta: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    x = 0.5 * (theta - np.sin(theta))
    y = -0.5 * (1.0 - np.cos(theta))
    return x, y


def endpoint_from_newtonian(theta_end: float) -> Tuple[float, float]:
    x_end = 0.5 * (theta_end - math.sin(theta_end))
    depth_end = 0.5 * (1.0 - math.cos(theta_end))
    return x_end, depth_end


def phi_for_fixed_depth(rho: float, rho_ref: float, depth_end: float, ascending: bool) -> float:
    a_depth = rho_ref * depth_end
    numerator = 1.0 - (1.0 - a_depth) ** 2
    denominator = 1.0 - (1.0 - rho) ** 2
    if denominator <= 0.0:
        return float("nan")
    ratio = numerator / denominator
    if ratio < -1.0e-12 or ratio > 1.0 + 1.0e-12:
        return float("nan")
    base = math.asin(math.sqrt(min(1.0, max(0.0, ratio))))
    return math.pi - base if ascending else base


def exact_endpoint_x_fixed_scale(
    rho: float,
    mode: str,
    rho_ref: float,
    depth_end: float,
    ascending: bool,
) -> Tuple[float, float]:
    phi_end = phi_for_fixed_depth(rho, rho_ref, depth_end, ascending)
    if not math.isfinite(phi_end):
        return float("nan"), float("nan")
    eta = 1.0 - rho
    one_minus_eta2 = 1.0 - eta * eta
    if mode == "coordinate":
        m = one_minus_eta2 / (eta * eta)
    elif mode == "proper":
        m = one_minus_eta2
    else:
        raise ValueError("mode must be 'coordinate' or 'proper'")
    f_minus_e = elliptic_f_minus_e_scalar(phi_end, m)
    if not math.isfinite(f_minus_e):
        return float("nan"), float("nan")
    return eta / rho_ref * f_minus_e, phi_end


def solve_rho_for_fixed_endpoint(
    mode: str,
    rho_ref: float,
    x_end: float,
    depth_end: float,
    ascending: bool,
) -> Tuple[float, float]:
    if not (0.0 < rho_ref < 1.0):
        raise ValueError("rho_ref must satisfy 0 < rho_ref < 1")
    min_rho = max(rho_ref * depth_end, 1.0e-8) + 1.0e-8
    max_rho = 0.285 if mode == "coordinate" else 0.95

    def residual(rho: float) -> float:
        x_value, _phi = exact_endpoint_x_fixed_scale(rho, mode, rho_ref, depth_end, ascending)
        if not math.isfinite(x_value):
            return float("nan")
        return x_value - x_end

    grid = np.linspace(min_rho, max_rho, 500)
    previous_rho = None
    previous_value = None
    candidates = []
    for rho in grid:
        value = residual(float(rho))
        if not math.isfinite(value):
            continue
        if previous_value is not None and value * previous_value <= 0.0:
            candidates.append((previous_rho, float(rho)))
        previous_rho = float(rho)
        previous_value = value

    if not candidates:
        raise RuntimeError(f"Could not bracket a fixed-endpoint solution for mode={mode!r}.")

    low, high = min(candidates, key=lambda pair: abs(0.5 * (pair[0] + pair[1]) - rho_ref))
    f_low = residual(low)
    for _ in range(80):
        mid = 0.5 * (low + high)
        f_mid = residual(mid)
        if not math.isfinite(f_mid):
            high = mid
            continue
        if abs(f_mid) < 1.0e-12:
            low = high = mid
            break
        if f_low * f_mid <= 0.0:
            high = mid
        else:
            low = mid
            f_low = f_mid

    rho_solution = 0.5 * (low + high)
    _x_solution, phi_solution = exact_endpoint_x_fixed_scale(
        rho_solution, mode, rho_ref, depth_end, ascending
    )
    return rho_solution, phi_solution


def pn_curve(theta: np.ndarray, rho: float, mode: str, order: str) -> Tuple[np.ndarray, np.ndarray]:
    x = pn_x_over_H(theta, rho, mode, order)
    y = -0.5 * (1.0 - np.cos(theta))
    return x, y


def pn_x_over_H(theta: np.ndarray, rho: float, mode: str, order: str) -> np.ndarray:
    base = theta - np.sin(theta)
    x = 0.5 * base
    if order in ("nlo", "nnlo"):
        if mode == "coordinate":
            x = x + rho * 5.0 / 8.0 * base
        elif mode == "proper":
            x = x - rho * 3.0 / 8.0 * base
        else:
            raise ValueError("mode must be 'coordinate' or 'proper'")
    if order == "nnlo":
        if mode == "coordinate":
            x = x + rho**2 * (
                47.0 / 32.0 * theta
                - 111.0 / 64.0 * np.sin(theta)
                + 17.0 / 128.0 * np.sin(2.0 * theta)
            )
        elif mode == "proper":
            x = x + rho**2 * (
                -1.0 / 32.0 * theta
                + 1.0 / 64.0 * np.sin(theta)
                + 1.0 / 128.0 * np.sin(2.0 * theta)
            )
    return x


def pn_endpoint_residual(
    theta_end: float,
    mode: str,
    order: str,
    rho_ref: float,
    x_end: float,
    depth_end: float,
) -> Tuple[float, float]:
    depth_over_h = 0.5 * (1.0 - math.cos(theta_end))
    if depth_over_h <= 0.0:
        return float("nan"), float("nan")
    h_over_h0 = depth_end / depth_over_h
    rho_h = rho_ref * h_over_h0
    x_over_h = float(pn_x_over_H(np.array([theta_end]), rho_h, mode, order)[0])
    return h_over_h0 * x_over_h - x_end, h_over_h0


def solve_pn_for_fixed_endpoint(
    mode: str,
    order: str,
    rho_ref: float,
    theta_reference: float,
    x_end: float,
    depth_end: float,
) -> Tuple[float, float]:
    ascending = theta_reference > math.pi
    if ascending:
        theta_min = math.pi + 1.0e-6
        theta_max = 2.0 * math.pi - 1.0e-6
    else:
        theta_min = 1.0e-6
        theta_max = math.pi - 1.0e-6

    def residual(theta_value: float) -> float:
        value, _h_over_h0 = pn_endpoint_residual(
            theta_value, mode, order, rho_ref, x_end, depth_end
        )
        return value

    grid = np.linspace(theta_min, theta_max, 800)
    previous_theta = None
    previous_value = None
    candidates = []
    for theta_value in grid:
        value = residual(float(theta_value))
        if not math.isfinite(value):
            continue
        if previous_value is not None and value * previous_value <= 0.0:
            candidates.append((previous_theta, float(theta_value)))
        previous_theta = float(theta_value)
        previous_value = value

    if not candidates:
        raise RuntimeError(f"Could not bracket endpoint PN solution for {mode=} and {order=}.")

    low, high = min(candidates, key=lambda pair: abs(0.5 * (pair[0] + pair[1]) - theta_reference))
    f_low = residual(low)
    for _ in range(80):
        mid = 0.5 * (low + high)
        f_mid = residual(mid)
        if abs(f_mid) < 1.0e-13:
            low = high = mid
            break
        if f_low * f_mid <= 0.0:
            high = mid
        else:
            low = mid
            f_low = f_mid

    theta_solution = 0.5 * (low + high)
    _residual_value, h_over_h0 = pn_endpoint_residual(
        theta_solution, mode, order, rho_ref, x_end, depth_end
    )
    return h_over_h0, theta_solution


def pn_curve_fixed_endpoint(
    mode: str,
    order: str,
    rho_ref: float,
    theta_reference: float,
    x_end: float,
    depth_end: float,
    n: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    h_over_h0, theta_end = solve_pn_for_fixed_endpoint(
        mode, order, rho_ref, theta_reference, x_end, depth_end
    )
    theta = np.linspace(0.0, theta_end, n)
    rho_h = rho_ref * h_over_h0
    x = h_over_h0 * pn_x_over_H(theta, rho_h, mode, order)
    y = -h_over_h0 * 0.5 * (1.0 - np.cos(theta))
    return theta, x, y, h_over_h0


def effective_indices(rho: float, n: int = 1200) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    z = np.linspace(0.01, 1.0, n)
    delta = rho * z
    lapse = 1.0 - delta
    denom = np.sqrt(np.maximum(1.0 - lapse * lapse, 1.0e-15))
    n_coord = 1.0 / (lapse * denom)
    n_proper = lapse / denom
    n_newton = 1.0 / np.sqrt(2.0 * delta)
    return z, n_coord, n_proper, n_newton


def save_csv(path: str, columns: Sequence[Tuple[str, np.ndarray]]) -> None:
    ensure_dir(os.path.dirname(path))
    length = len(columns[0][1])
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow([name for name, _data in columns])
        for index in range(length):
            writer.writerow([f"{float(data[index]):.12g}" for _name, data in columns])


def nice_ticks(lo: float, hi: float, target: int = 6) -> List[float]:
    if not math.isfinite(lo) or not math.isfinite(hi) or hi <= lo:
        return [lo, hi]
    span = hi - lo
    raw = span / max(target - 1, 1)
    exp = math.floor(math.log10(raw))
    frac = raw / 10**exp
    if frac <= 1.5:
        step = 1.0 * 10**exp
    elif frac <= 3.0:
        step = 2.0 * 10**exp
    elif frac <= 7.0:
        step = 5.0 * 10**exp
    else:
        step = 10.0 * 10**exp
    start = math.ceil(lo / step) * step
    ticks = []
    value = start
    while value <= hi + 0.5 * step:
        if value >= lo - 1.0e-12:
            ticks.append(0.0 if abs(value) < 1.0e-12 else value)
        value += step
    return ticks


def fmt_tick(value: float) -> str:
    if abs(value) < 1.0e-12:
        return "0"
    if abs(value) >= 100 or abs(value) < 0.01:
        return f"{value:.1e}"
    if abs(value - round(value)) < 1.0e-9:
        return str(int(round(value)))
    return f"{value:.2f}".rstrip("0").rstrip(".")


def dash_pattern(style: str) -> Optional[List[float]]:
    if style == "dash":
        return [8.0, 5.0]
    if style == "dot":
        return [2.0, 4.0]
    if style == "dashdot":
        return [8.0, 4.0, 2.0, 4.0]
    return None


def hex_to_rgb01(value: str) -> Tuple[float, float, float]:
    value = value.strip().lstrip("#")
    return tuple(int(value[i : i + 2], 16) / 255.0 for i in (0, 2, 4))


def adjusted_limits(panel: Panel, width: float, height: float) -> Tuple[Tuple[float, float], Tuple[float, float]]:
    xlim = panel.xlim
    ylim = panel.ylim
    if not panel.equal_aspect:
        return xlim, ylim
    xrange = xlim[1] - xlim[0]
    yrange = ylim[1] - ylim[0]
    if xrange <= 0.0 or yrange <= 0.0 or width <= 0.0 or height <= 0.0:
        return xlim, ylim
    data_ratio = xrange / yrange
    box_ratio = width / height
    if data_ratio > box_ratio:
        new_yrange = xrange / box_ratio
        center = 0.5 * (ylim[0] + ylim[1])
        ylim = (center - 0.5 * new_yrange, center + 0.5 * new_yrange)
    else:
        new_xrange = yrange * box_ratio
        center = 0.5 * (xlim[0] + xlim[1])
        xlim = (center - 0.5 * new_xrange, center + 0.5 * new_xrange)
    return xlim, ylim


def transform(
    x: np.ndarray,
    y: np.ndarray,
    rect: Tuple[float, float, float, float],
    xlim: Tuple[float, float],
    ylim: Tuple[float, float],
) -> Tuple[np.ndarray, np.ndarray]:
    left, bottom, width, height = rect
    sx = left + (x - xlim[0]) / (xlim[1] - xlim[0]) * width
    sy = bottom + (y - ylim[0]) / (ylim[1] - ylim[0]) * height
    return sx, sy


def _pil_load_font(size_pt: float) -> "ImageFont.FreeTypeFont":
    from PIL import ImageFont

    size_px = max(1, round(size_pt * PNG_DPI / 72))
    for path in (
        "/System/Library/Fonts/Helvetica.ttc",
        "/System/Library/Fonts/HelveticaNeue.ttc",
        "/Library/Fonts/Arial.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/liberation/LiberationSans-Regular.ttf",
    ):
        if os.path.exists(path):
            try:
                return ImageFont.truetype(path, size_px)
            except Exception:
                pass
    return ImageFont.load_default()


def _draw_polyline_pil(
    draw: "ImageDraw.ImageDraw",
    xs: np.ndarray,
    ys: np.ndarray,
    fill: Tuple[int, int, int],
    width: int,
    dash: Optional[List[float]],
) -> None:
    n = len(xs)
    if n < 2:
        return
    if dash is None:
        draw.line(list(zip(xs.tolist(), ys.tolist())), fill=fill, width=width)
        return
    on = True
    pat_idx = 0
    pat_pos = 0.0
    pat_n = len(dash)
    for i in range(n - 1):
        x1, y1 = float(xs[i]), float(ys[i])
        x2, y2 = float(xs[i + 1]), float(ys[i + 1])
        seg = math.hypot(x2 - x1, y2 - y1)
        if seg < 0.5:
            continue
        dx, dy = (x2 - x1) / seg, (y2 - y1) / seg
        pos = 0.0
        while pos < seg - 1e-9:
            avail = dash[pat_idx % pat_n] - pat_pos
            step = min(avail, seg - pos)
            if on:
                sx, sy = x1 + dx * pos, y1 + dy * pos
                ex, ey = x1 + dx * (pos + step), y1 + dy * (pos + step)
                draw.line([(sx, sy), (ex, ey)], fill=fill, width=width)
            pos += step
            pat_pos += step
            if pat_pos >= dash[pat_idx % pat_n] - 1e-9:
                pat_pos = 0.0
                pat_idx += 1
                on = not on


def svg_polyline(points: Iterable[Point]) -> str:
    return " ".join(f"{x:.2f},{y:.2f}" for x, y in points)


def write_svg(path: str, panels: Sequence[Panel], size: Tuple[int, int]) -> None:
    width, height = size
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        "<style>",
        "text{font-family:Helvetica,Arial,sans-serif;fill:#111}",
        ".small{font-size:13px}.label{font-size:15px}.title{font-size:16px;font-weight:600}",
        "</style>",
        f'<rect x="0" y="0" width="{width}" height="{height}" fill="white"/>',
    ]

    for panel_index, panel in enumerate(panels):
        left = panel.rect[0] * width
        bottom = panel.rect[1] * height
        panel_width = panel.rect[2] * width
        panel_height = panel.rect[3] * height
        rect_svg = (left, height - bottom - panel_height, panel_width, panel_height)
        xlim, ylim = adjusted_limits(panel, panel_width, panel_height)
        xticks = list(panel.xticks) if panel.xticks is not None else nice_ticks(*xlim)
        yticks = list(panel.yticks) if panel.yticks is not None else nice_ticks(*ylim)
        clip_id = f"clip{panel_index}"

        parts.append(
            f'<clipPath id="{clip_id}"><rect x="{rect_svg[0]:.2f}" y="{rect_svg[1]:.2f}" '
            f'width="{rect_svg[2]:.2f}" height="{rect_svg[3]:.2f}"/></clipPath>'
        )

        for tick in xticks:
            px, _py = transform(
                np.array([tick]), np.array([ylim[0]]), (rect_svg[0], rect_svg[1], rect_svg[2], rect_svg[3]), xlim, ylim
            )
            x = px[0]
            parts.append(
                f'<line x1="{x:.2f}" y1="{rect_svg[1]:.2f}" x2="{x:.2f}" '
                f'y2="{rect_svg[1] + rect_svg[3]:.2f}" stroke="#e7e7e7" stroke-width="1"/>'
            )
            parts.append(
                f'<line x1="{x:.2f}" y1="{rect_svg[1] + rect_svg[3]:.2f}" '
                f'x2="{x:.2f}" y2="{rect_svg[1] + rect_svg[3] + 5:.2f}" stroke="#111" stroke-width="1"/>'
            )
            parts.append(
                f'<text class="small" x="{x:.2f}" y="{rect_svg[1] + rect_svg[3] + 21:.2f}" text-anchor="middle">'
                f"{fmt_tick(tick)}</text>"
            )

        for tick in yticks:
            _px, py = transform(
                np.array([xlim[0]]), np.array([tick]), (rect_svg[0], rect_svg[1], rect_svg[2], rect_svg[3]), xlim, ylim
            )
            y = rect_svg[1] + rect_svg[3] - (py[0] - rect_svg[1])
            parts.append(
                f'<line x1="{rect_svg[0]:.2f}" y1="{y:.2f}" x2="{rect_svg[0] + rect_svg[2]:.2f}" '
                f'y2="{y:.2f}" stroke="#e7e7e7" stroke-width="1"/>'
            )
            parts.append(
                f'<line x1="{rect_svg[0] - 5:.2f}" y1="{y:.2f}" x2="{rect_svg[0]:.2f}" '
                f'y2="{y:.2f}" stroke="#111" stroke-width="1"/>'
            )
            parts.append(
                f'<text class="small" x="{rect_svg[0] - 9:.2f}" y="{y + 4:.2f}" text-anchor="end">'
                f"{fmt_tick(tick)}</text>"
            )

        parts.append(
            f'<rect x="{rect_svg[0]:.2f}" y="{rect_svg[1]:.2f}" width="{rect_svg[2]:.2f}" '
            f'height="{rect_svg[3]:.2f}" fill="none" stroke="#111" stroke-width="1.1"/>'
        )

        for item in panel.series:
            x, y = transform(item.x, item.y, (rect_svg[0], rect_svg[1], rect_svg[2], rect_svg[3]), xlim, ylim)
            y = rect_svg[1] + rect_svg[3] - (y - rect_svg[1])
            dash = dash_pattern(item.style)
            dash_attr = "" if dash is None else f' stroke-dasharray="{",".join(str(d) for d in dash)}"'
            parts.append(
                f'<polyline points="{svg_polyline(zip(x, y))}" fill="none" stroke="{item.color}" '
                f'stroke-width="{item.width:.2f}" stroke-linejoin="round" stroke-linecap="round"'
                f'{dash_attr} clip-path="url(#{clip_id})"/>'
            )

        parts.append(
            f'<text class="title" x="{rect_svg[0] + 0.5 * rect_svg[2]:.2f}" y="{rect_svg[1] - 14:.2f}" '
            f'text-anchor="middle">{panel.title}</text>'
        )
        parts.append(
            f'<text class="label" x="{rect_svg[0] + 0.5 * rect_svg[2]:.2f}" '
            f'y="{rect_svg[1] + rect_svg[3] + 48:.2f}" text-anchor="middle">{panel.xlabel}</text>'
        )
        ylabel_x = rect_svg[0] - 58
        ylabel_y = rect_svg[1] + 0.5 * rect_svg[3]
        parts.append(
            f'<text class="label" x="{ylabel_x:.2f}" y="{ylabel_y:.2f}" text-anchor="middle" '
            f'transform="rotate(-90 {ylabel_x:.2f} {ylabel_y:.2f})">{panel.ylabel}</text>'
        )

        if panel.legend and panel.series:
            max_label = max(len(item.label) for item in panel.series)
            leg_width = max(120, 48 + 7 * max_label)
            leg_height = 22 + 20 * len(panel.series)
            leg_x = rect_svg[0] + rect_svg[2] - leg_width - 12
            leg_y = rect_svg[1] + 12
            parts.append(
                f'<rect x="{leg_x:.2f}" y="{leg_y:.2f}" width="{leg_width:.2f}" height="{leg_height:.2f}" '
                f'fill="white" stroke="#cfcfcf" stroke-width="1" opacity="0.94"/>'
            )
            for i, item in enumerate(panel.series):
                y0 = leg_y + 18 + 20 * i
                dash = dash_pattern(item.style)
                dash_attr = "" if dash is None else f' stroke-dasharray="{",".join(str(d) for d in dash)}"'
                parts.append(
                    f'<line x1="{leg_x + 11:.2f}" y1="{y0:.2f}" x2="{leg_x + 41:.2f}" y2="{y0:.2f}" '
                    f'stroke="{item.color}" stroke-width="{item.width:.2f}" stroke-linecap="round"{dash_attr}/>'
                )
                parts.append(
                    f'<text class="small" x="{leg_x + 49:.2f}" y="{y0 + 4:.2f}">{item.label}</text>'
                )

    parts.append("</svg>")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(parts))


def pdf_line(canvas_obj, points: Sequence[Point]) -> None:
    if len(points) < 2:
        return
    path = canvas_obj.beginPath()
    path.moveTo(points[0][0], points[0][1])
    for x, y in points[1:]:
        path.lineTo(x, y)
    canvas_obj.drawPath(path, stroke=1, fill=0)


def write_pdf(path: str, panels: Sequence[Panel], size: Tuple[float, float]) -> None:
    if not HAS_REPORTLAB:
        return
    width, height = size
    doc = canvas.Canvas(path, pagesize=(width, height))
    doc.setTitle(os.path.splitext(os.path.basename(path))[0])
    doc.setFillColor(colors.white)
    doc.rect(0, 0, width, height, stroke=0, fill=1)

    for panel in panels:
        left = panel.rect[0] * width
        bottom = panel.rect[1] * height
        panel_width = panel.rect[2] * width
        panel_height = panel.rect[3] * height
        xlim, ylim = adjusted_limits(panel, panel_width, panel_height)
        xticks = list(panel.xticks) if panel.xticks is not None else nice_ticks(*xlim)
        yticks = list(panel.yticks) if panel.yticks is not None else nice_ticks(*ylim)

        def tx(values: np.ndarray) -> np.ndarray:
            return left + (values - xlim[0]) / (xlim[1] - xlim[0]) * panel_width

        def ty(values: np.ndarray) -> np.ndarray:
            return bottom + (values - ylim[0]) / (ylim[1] - ylim[0]) * panel_height

        doc.setStrokeColor(colors.HexColor("#e7e7e7"))
        doc.setLineWidth(0.5)
        for tick in xticks:
            x = float(tx(np.array([tick]))[0])
            doc.line(x, bottom, x, bottom + panel_height)
        for tick in yticks:
            y = float(ty(np.array([tick]))[0])
            doc.line(left, y, left + panel_width, y)

        doc.setStrokeColor(colors.HexColor("#111111"))
        doc.setLineWidth(0.8)
        doc.rect(left, bottom, panel_width, panel_height, stroke=1, fill=0)
        doc.setFont("Helvetica", 8.5)
        for tick in xticks:
            x = float(tx(np.array([tick]))[0])
            doc.line(x, bottom, x, bottom - 3)
            doc.drawCentredString(x, bottom - 14, fmt_tick(tick))
        for tick in yticks:
            y = float(ty(np.array([tick]))[0])
            doc.line(left - 3, y, left, y)
            doc.drawRightString(left - 5, y - 3, fmt_tick(tick))

        doc.setFont("Helvetica-Bold", 10.5)
        doc.drawCentredString(left + 0.5 * panel_width, bottom + panel_height + 16, panel.title)
        doc.setFont("Helvetica", 9.5)
        doc.drawCentredString(left + 0.5 * panel_width, bottom - 31, panel.xlabel)
        doc.saveState()
        doc.translate(left - 42, bottom + 0.5 * panel_height)
        doc.rotate(90)
        doc.drawCentredString(0, 0, panel.ylabel)
        doc.restoreState()

        for item in panel.series:
            doc.setStrokeColor(colors.Color(*hex_to_rgb01(item.color)))
            doc.setLineWidth(item.width * 0.75)
            pattern = dash_pattern(item.style)
            if pattern:
                doc.setDash(pattern)
            else:
                doc.setDash()
            xs = tx(item.x)
            ys = ty(item.y)
            points = list(zip(xs.tolist(), ys.tolist()))
            pdf_line(doc, points)
        doc.setDash()

        if panel.legend and panel.series:
            max_label = max(len(item.label) for item in panel.series)
            leg_width = max(88, 34 + 4.9 * max_label)
            leg_height = 16 + 13.5 * len(panel.series)
            leg_x = left + panel_width - leg_width - 8
            leg_y = bottom + panel_height - leg_height - 8
            doc.setFillColor(colors.white)
            doc.setStrokeColor(colors.HexColor("#cfcfcf"))
            doc.setLineWidth(0.5)
            doc.rect(leg_x, leg_y, leg_width, leg_height, stroke=1, fill=1)
            doc.setFont("Helvetica", 8)
            doc.setFillColor(colors.HexColor("#111111"))
            for i, item in enumerate(panel.series):
                y0 = leg_y + leg_height - 12 - 13.5 * i
                doc.setStrokeColor(colors.Color(*hex_to_rgb01(item.color)))
                doc.setLineWidth(item.width * 0.75)
                pattern = dash_pattern(item.style)
                if pattern:
                    doc.setDash(pattern)
                else:
                    doc.setDash()
                doc.line(leg_x + 7, y0, leg_x + 27, y0)
                doc.setDash()
                doc.setFillColor(colors.HexColor("#111111"))
                doc.drawString(leg_x + 33, y0 - 3, item.label)

    doc.showPage()
    doc.save()


def write_png(path: str, panels: Sequence[Panel], size_pt: Tuple[float, float], dpi: int = PNG_DPI) -> None:
    if not HAS_PIL:
        return
    from PIL import Image, ImageDraw

    s = dpi / 72.0
    wp = round(size_pt[0] * s)
    hp = round(size_pt[1] * s)

    img = Image.new("RGB", (wp, hp), (255, 255, 255))
    draw = ImageDraw.Draw(img)
    font_sm = _pil_load_font(8.5)
    font_lbl = _pil_load_font(9.5)
    font_ttl = _pil_load_font(10.5)
    black = (17, 17, 17)
    grid_col = (231, 231, 231)

    for panel in panels:
        pl = panel.rect[0] * wp
        pw = panel.rect[2] * wp
        ph = panel.rect[3] * hp
        pt = hp - (panel.rect[1] + panel.rect[3]) * hp

        xlim, ylim = adjusted_limits(panel, pw, ph)
        xticks = list(panel.xticks) if panel.xticks is not None else nice_ticks(*xlim)
        yticks = list(panel.yticks) if panel.yticks is not None else nice_ticks(*ylim)

        def _tx(v: np.ndarray, _pl: float = pl, _pw: float = pw, _xl: Tuple = xlim) -> np.ndarray:
            return _pl + (v - _xl[0]) / (_xl[1] - _xl[0]) * _pw

        def _ty(v: np.ndarray, _pt: float = pt, _ph: float = ph, _yl: Tuple = ylim) -> np.ndarray:
            return _pt + _ph - (v - _yl[0]) / (_yl[1] - _yl[0]) * _ph

        for tick in xticks:
            x = float(_tx(np.array([tick]))[0])
            draw.line([(x, pt), (x, pt + ph)], fill=grid_col, width=1)
        for tick in yticks:
            y = float(_ty(np.array([tick]))[0])
            draw.line([(pl, y), (pl + pw, y)], fill=grid_col, width=1)

        for item in panel.series:
            xs_px = _tx(item.x)
            ys_px = _ty(item.y)
            rgb01 = hex_to_rgb01(item.color)
            col = (round(rgb01[0] * 255), round(rgb01[1] * 255), round(rgb01[2] * 255))
            lw = max(1, round(item.width * 0.75 * s))
            raw_dash = dash_pattern(item.style)
            scaled_dash = [d * s for d in raw_dash] if raw_dash else None
            _draw_polyline_pil(draw, xs_px, ys_px, col, lw, scaled_dash)

        draw.rectangle([(pl, pt), (pl + pw, pt + ph)], outline=black, width=1)

        for tick in xticks:
            x = float(_tx(np.array([tick]))[0])
            draw.line([(x, pt + ph), (x, pt + ph + s * 3)], fill=black, width=1)
            label = fmt_tick(tick)
            bb = draw.textbbox((0, 0), label, font=font_sm)
            tw, th = bb[2] - bb[0], bb[3] - bb[1]
            draw.text((x - tw / 2, pt + ph + s * 5), label, font=font_sm, fill=black)

        for tick in yticks:
            y = float(_ty(np.array([tick]))[0])
            draw.line([(pl - s * 3, y), (pl, y)], fill=black, width=1)
            label = fmt_tick(tick)
            bb = draw.textbbox((0, 0), label, font=font_sm)
            tw, th = bb[2] - bb[0], bb[3] - bb[1]
            draw.text((pl - s * 5 - tw, y - th / 2), label, font=font_sm, fill=black)

        bb = draw.textbbox((0, 0), panel.xlabel, font=font_lbl)
        tw, th = bb[2] - bb[0], bb[3] - bb[1]
        draw.text((pl + pw / 2 - tw / 2, pt + ph + s * 31 - th / 2), panel.xlabel, font=font_lbl, fill=black)

        from PIL import Image as _PILImg
        bb = draw.textbbox((0, 0), panel.ylabel, font=font_lbl)
        tw, th = bb[2] - bb[0], bb[3] - bb[1]
        tmp = _PILImg.new("RGBA", (round(tw) + 4, round(th) + 4), (0, 0, 0, 0))
        tdraw = ImageDraw.Draw(tmp)
        tdraw.text((2, 2), panel.ylabel, font=font_lbl, fill=black)
        tmp = tmp.rotate(90, expand=True)
        img.paste(tmp, (round(pl - s * 58), round(pt + ph / 2 - tmp.height / 2)), tmp)

        if panel.title:
            bb = draw.textbbox((0, 0), panel.title, font=font_ttl)
            tw, th = bb[2] - bb[0], bb[3] - bb[1]
            draw.text((pl + pw / 2 - tw / 2, pt - s * 14 - th), panel.title, font=font_ttl, fill=black)

        if panel.legend and panel.series:
            max_lbl_len = max(len(item.label) for item in panel.series)
            leg_w = max(s * 88, s * 34 + 4.9 * s * max_lbl_len)
            leg_h = s * 16 + len(panel.series) * s * 13.5
            leg_x = pl + pw - leg_w - s * 8
            leg_y = pt + s * 8
            draw.rectangle(
                [(leg_x, leg_y), (leg_x + leg_w, leg_y + leg_h)],
                fill=(255, 255, 255),
                outline=(207, 207, 207),
                width=1,
            )
            for i, item in enumerate(panel.series):
                y0 = leg_y + s * 12 + i * s * 13.5
                rgb01 = hex_to_rgb01(item.color)
                col = (round(rgb01[0] * 255), round(rgb01[1] * 255), round(rgb01[2] * 255))
                lw = max(1, round(item.width * 0.75 * s))
                raw_dash = dash_pattern(item.style)
                scaled_dash = [d * s for d in raw_dash] if raw_dash else None
                _draw_polyline_pil(
                    draw,
                    np.array([leg_x + s * 7, leg_x + s * 27]),
                    np.array([y0, y0]),
                    col,
                    lw,
                    scaled_dash,
                )
                bb = draw.textbbox((0, 0), item.label, font=font_sm)
                th = bb[3] - bb[1]
                draw.text((leg_x + s * 33, y0 - th / 2), item.label, font=font_sm, fill=black)

    img.save(path, "PNG", dpi=(dpi, dpi))


def save_figure(
    outdir: str,
    basename: str,
    panels: Sequence[Panel],
    svg_size: Tuple[int, int],
    pdf_size: Tuple[float, float],
) -> None:
    ensure_dir(outdir)
    write_svg(os.path.join(outdir, f"{basename}.svg"), panels, svg_size)
    if HAS_REPORTLAB:
        write_pdf(os.path.join(outdir, f"{basename}.pdf"), panels, pdf_size)
    if HAS_PIL:
        write_png(os.path.join(outdir, f"{basename}.png"), panels, pdf_size)


def limits_from_series(series: Sequence[Series], pad_x: float = 0.04, pad_y: float = 0.08) -> Tuple[Tuple[float, float], Tuple[float, float]]:
    xs = np.concatenate([item.x[np.isfinite(item.x)] for item in series])
    ys = np.concatenate([item.y[np.isfinite(item.y)] for item in series])
    xmin, xmax = float(np.min(xs)), float(np.max(xs))
    ymin, ymax = float(np.min(ys)), float(np.max(ys))
    xrange = xmax - xmin
    yrange = ymax - ymin
    return (
        (xmin - pad_x * xrange, xmax + pad_x * xrange),
        (ymin - pad_y * yrange, ymax + pad_y * yrange),
    )


def exact_fixed_endpoint_curve_for_mode(
    mode: str,
    rho_ref: float,
    theta_end_deg: float,
    n: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    theta_end = math.radians(theta_end_deg)
    ascending = theta_end > math.pi
    x_end, depth_end = endpoint_from_newtonian(theta_end)
    rho_solution, phi_end = solve_rho_for_fixed_endpoint(
        mode, rho_ref, x_end, depth_end, ascending
    )
    phi, x, y = exact_curve_fixed_scale(rho_solution, mode, rho_ref, phi_end, n)
    return phi, x, y, rho_solution


def branch_depth_x(phi: np.ndarray, x: np.ndarray, y: np.ndarray) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
    depth = -y
    turn_index = int(np.argmax(depth))
    branches = {
        "descending": (depth[: turn_index + 1], x[: turn_index + 1]),
    }
    if turn_index < len(depth) - 1:
        branches["ascending"] = (depth[turn_index:][::-1], x[turn_index:][::-1])
    return branches


def same_depth_separation(
    rho_ref: float,
    theta_end_deg: float,
    n: int,
    samples: int = 900,
) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
    phi_coord, x_coord, y_coord, _rho_coord = exact_fixed_endpoint_curve_for_mode(
        "coordinate", rho_ref, theta_end_deg, n
    )
    phi_proper, x_proper, y_proper, _rho_proper = exact_fixed_endpoint_curve_for_mode(
        "proper", rho_ref, theta_end_deg, n
    )
    coord_branches = branch_depth_x(phi_coord, x_coord, y_coord)
    proper_branches = branch_depth_x(phi_proper, x_proper, y_proper)
    out = {}
    for branch in sorted(set(coord_branches) & set(proper_branches)):
        coord_depth, coord_x = coord_branches[branch]
        proper_depth, proper_x = proper_branches[branch]
        lower = max(float(np.min(coord_depth)), float(np.min(proper_depth)))
        upper = min(float(np.max(coord_depth)), float(np.max(proper_depth)))
        if upper <= lower:
            continue
        depth_grid = np.linspace(lower, upper, samples)
        coord_interp = np.interp(depth_grid, coord_depth, coord_x)
        proper_interp = np.interp(depth_grid, proper_depth, proper_x)
        out[branch] = (depth_grid, coord_interp - proper_interp)
    return out


def figure_exact_comparison(
    outdir: str,
    datadir: str,
    rho: float,
    theta_end_deg: float,
    n: int,
) -> None:
    theta_end = math.radians(theta_end_deg)
    ascending = theta_end > math.pi
    x_end, depth_end = endpoint_from_newtonian(theta_end)

    theta = np.linspace(0.0, theta_end, n)
    x_newton, y_newton = newtonian_curve(theta)

    rho_coord, phi_coord_end = solve_rho_for_fixed_endpoint(
        "coordinate", rho, x_end, depth_end, ascending
    )
    rho_proper, phi_proper_end = solve_rho_for_fixed_endpoint(
        "proper", rho, x_end, depth_end, ascending
    )
    phi_coord, x_coord, y_coord = exact_curve_fixed_scale(
        rho_coord, "coordinate", rho, phi_coord_end, n
    )
    phi_proper, x_proper, y_proper = exact_curve_fixed_scale(
        rho_proper, "proper", rho, phi_proper_end, n
    )

    series = [
        Series(x_newton, y_newton, "Newtonian", PALETTE["newton"], 2.1, "dot"),
        Series(x_coord, y_coord, "coordinate time", PALETTE["coord"], 2.3, "solid"),
        Series(x_proper, y_proper, "proper time", PALETTE["proper"], 2.3, "dash"),
    ]
    xlim, ylim = limits_from_series(series, pad_x=0.035, pad_y=0.10)
    panel = Panel(
        rect=(0.105, 0.18, 0.84, 0.68),
        xlim=xlim,
        ylim=ylim,
        xlabel="x / H0",
        ylabel="y / H0",
        title=f"Fixed endpoints, a H0 = {rho:g}, theta_B = {theta_end_deg:g} deg",
        series=series,
        equal_aspect=False,
    )
    base = f"fig_01_confronto_curve_exact_rho_{rho_tag(rho)}"
    save_figure(outdir, base, [panel], (1100, 520), (7.4 * 72, 3.5 * 72))
    save_csv(
        os.path.join(datadir, f"{base}.csv"),
        [
            ("theta", theta),
            ("x_newton_over_H0", x_newton),
            ("y_newton_over_H0", y_newton),
            ("phi_coordinate", phi_coord),
            ("x_coordinate_over_H0", x_coord),
            ("y_coordinate_over_H0", y_coord),
            ("rho_coordinate_H", np.full_like(theta, rho_coord)),
            ("phi_proper", phi_proper),
            ("x_proper_over_H0", x_proper),
            ("y_proper_over_H0", y_proper),
            ("rho_proper_H", np.full_like(theta, rho_proper)),
            ("x_endpoint_over_H0", np.full_like(theta, x_end)),
            ("depth_endpoint_over_H0", np.full_like(theta, depth_end)),
        ],
    )


def figure_separation(
    outdir: str,
    datadir: str,
    rho: float,
    theta_end_deg: float,
    rhos: Sequence[float],
    n: int,
) -> None:
    branch_sep = same_depth_separation(rho, theta_end_deg, n)
    left_series = [
        Series(depth, delta_x, branch, color, 2.1, style)
        for (branch, (depth, delta_x)), color, style in zip(
            branch_sep.items(),
            [PALETTE["coord"], PALETTE["proper"]],
            ["solid", "dash"],
        )
    ]

    rho_grid = np.linspace(0.005, max(max(rhos), 0.24), 90)
    max_sep = []
    for rho_value in rho_grid:
        separations = same_depth_separation(float(rho_value), theta_end_deg, n, samples=500)
        if not separations:
            max_sep.append(float("nan"))
            continue
        values = np.concatenate([np.abs(delta_x) for _depth, delta_x in separations.values()])
        max_sep.append(float(np.nanmax(values)))
    max_sep = np.array(max_sep)
    right_series = [
        Series(rho_grid, max_sep, "fixed endpoints", PALETTE["coord"], 2.2),
    ]

    left_xlim, left_ylim = limits_from_series(left_series, pad_x=0.02, pad_y=0.12)
    right_xlim, right_ylim = limits_from_series(right_series, pad_x=0.03, pad_y=0.06)
    panels = [
        Panel(
            rect=(0.085, 0.18, 0.40, 0.66),
            xlim=left_xlim,
            ylim=left_ylim,
            xlabel="h / H0",
            ylabel="(x_t - x_tau) / H0",
            title=f"Same-depth separation, a H0 = {rho:g}",
            series=left_series,
        ),
        Panel(
            rect=(0.595, 0.18, 0.345, 0.66),
            xlim=right_xlim,
            ylim=right_ylim,
            xlabel="a H0",
            ylabel="max |x_t - x_tau| / H0",
            title="Endpoint-fixed growth",
            series=right_series,
        ),
    ]
    base = "fig_02_separazione_al_variare_di_rho"
    save_figure(outdir, base, panels, (1200, 560), (8.0 * 72, 3.75 * 72))

    columns = []
    for item in left_series:
        columns.append((f"h_over_H0_{item.label}", item.x))
        columns.append((f"delta_x_over_H0_{item.label}", item.y))
    save_csv(os.path.join(datadir, f"{base}_pointwise.csv"), columns)
    save_csv(
        os.path.join(datadir, f"{base}_max.csv"),
        [("aH0", rho_grid), ("max_abs_delta_x_over_H0_exact_fixed_endpoint", max_sep)],
    )


def figure_exact_vs_pn(
    outdir: str,
    datadir: str,
    rho: float,
    theta_end_deg: float,
    mode: str,
    n: int,
) -> None:
    theta_end = math.radians(theta_end_deg)
    ascending = theta_end > math.pi
    x_end, depth_end = endpoint_from_newtonian(theta_end)
    theta_newton = np.linspace(0.0, theta_end, n)
    x_newton, y_newton = newtonian_curve(theta_newton)

    rho_exact, phi_exact_end = solve_rho_for_fixed_endpoint(
        mode, rho, x_end, depth_end, ascending
    )
    phi_exact, x_exact, y_exact = exact_curve_fixed_scale(
        rho_exact, mode, rho, phi_exact_end, n
    )
    theta_nlo, x_nlo, y_nlo, h_nlo = pn_curve_fixed_endpoint(
        mode, "nlo", rho, theta_end, x_end, depth_end, n
    )
    theta_nnlo, x_nnlo, y_nnlo, h_nnlo = pn_curve_fixed_endpoint(
        mode, "nnlo", rho, theta_end, x_end, depth_end, n
    )
    mode_label = "coordinate time" if mode == "coordinate" else "proper time"
    mode_tag = "tempo_coordinato" if mode == "coordinate" else "tempo_proprio"
    main_color = PALETTE["coord"] if mode == "coordinate" else PALETTE["proper"]
    series = [
        Series(x_exact, y_exact, "exact elliptic", main_color, 2.4),
        Series(x_newton, y_newton, "Newtonian", PALETTE["newton"], 1.9, "dot"),
        Series(x_nlo, y_nlo, "NLO", PALETTE["nlo"], 2.0, "dash"),
        Series(x_nnlo, y_nnlo, "NNLO", PALETTE["nnlo"], 2.0, "dashdot"),
    ]
    xlim, ylim = limits_from_series(series, pad_x=0.035, pad_y=0.10)
    panel = Panel(
        rect=(0.105, 0.18, 0.84, 0.68),
        xlim=xlim,
        ylim=ylim,
        xlabel="x / H0",
        ylabel="y / H0",
        title=f"Fixed-endpoint exact vs PN, {mode_label}",
        series=series,
    )
    base = f"fig_03_exact_vs_pn_{mode_tag}_rho_{rho_tag(rho)}"
    if mode == "proper":
        base = f"fig_04_exact_vs_pn_{mode_tag}_rho_{rho_tag(rho)}"
    save_figure(outdir, base, [panel], (1100, 520), (7.4 * 72, 3.5 * 72))
    save_csv(
        os.path.join(datadir, f"{base}.csv"),
        [
            ("theta_newton", theta_newton),
            ("x_newton_over_H0", x_newton),
            ("y_newton_over_H0", y_newton),
            ("phi_exact", phi_exact),
            ("x_exact_over_H0", x_exact),
            ("y_exact_over_H0", y_exact),
            ("rho_exact_H", np.full_like(theta_newton, rho_exact)),
            ("theta_nlo", theta_nlo),
            ("x_nlo_over_H0", x_nlo),
            ("y_nlo_over_H0", y_nlo),
            ("H_nlo_over_H0", np.full_like(theta_newton, h_nlo)),
            ("theta_nnlo", theta_nnlo),
            ("x_nnlo_over_H0", x_nnlo),
            ("y_nnlo_over_H0", y_nnlo),
            ("H_nnlo_over_H0", np.full_like(theta_newton, h_nnlo)),
            ("x_endpoint_over_H0", np.full_like(theta_newton, x_end)),
            ("depth_endpoint_over_H0", np.full_like(theta_newton, depth_end)),
        ],
    )


def figure_effective_indices(outdir: str, datadir: str, rho: float, n: int) -> None:
    z, n_coord, n_proper, n_newton = effective_indices(rho, n)
    series = [
        Series(z, n_coord / n_newton, "n_t / n_N", PALETTE["coord"], 2.2),
        Series(z, n_proper / n_newton, "n_tau / n_N", PALETTE["proper"], 2.2, "dash"),
        Series(z, np.ones_like(z), "Newtonian", PALETTE["newton"], 1.8, "dot"),
    ]
    xlim, ylim = limits_from_series(series, pad_x=0.02, pad_y=0.06)
    panel = Panel(
        rect=(0.105, 0.18, 0.84, 0.68),
        xlim=(0.0, 1.0),
        ylim=ylim,
        xlabel="h / H0",
        ylabel="effective index / Newtonian index",
        title=f"Relative Fermat weights, a H0 = {rho:g}",
        series=series,
    )
    base = f"fig_05_indici_effettivi_rho_{rho_tag(rho)}"
    save_figure(outdir, base, [panel], (1100, 600), (7.4 * 72, 4.0 * 72))
    save_csv(
        os.path.join(datadir, f"{base}.csv"),
        [
            ("h_over_H", z),
            ("n_coordinate", n_coord),
            ("n_proper", n_proper),
            ("n_newtonian_leading", n_newton),
            ("n_coordinate_over_newton", n_coord / n_newton),
            ("n_proper_over_newton", n_proper / n_newton),
        ],
    )


def write_readme(
    outdir: str,
    datadir: str,
    rho: float,
    rhos: Sequence[float],
    theta_end_deg: float,
) -> None:
    theta_end = math.radians(theta_end_deg)
    x_end, depth_end = endpoint_from_newtonian(theta_end)
    path = os.path.join(outdir, "README.md")
    lines = [
        "# Rindler brachistochrone figures",
        "",
        "Generated by `RindlerMetric/genera_grafici_rindler.py`.",
        "",
        "Main parameters:",
        f"- Reference acceleration scale `a H0`: `{rho:g}`",
        f"- Figure 01 common endpoint: `theta_B = {theta_end_deg:g} deg`, `L/H0 = {x_end:.8g}`, `D/H0 = {depth_end:.8g}`",
        f"- Separation rho values: `{', '.join(f'{value:g}' for value in rhos)}`",
        "",
        "Files:",
        "- `fig_01_*`: fixed-endpoint exact coordinate-time/proper-time curves vs Newtonian cycloid.",
        "- `fig_02_*`: same-depth coordinate/proper separation with fixed endpoints.",
        "- `fig_03_*`: endpoint-fixed exact coordinate-time curve vs NLO/NNLO expansion.",
        "- `fig_04_*`: endpoint-fixed exact proper-time curve vs NLO/NNLO expansion.",
        "- `fig_05_*`: normalized effective refractive indices.",
        "",
        f"Numerical source data are in `{os.path.relpath(datadir, outdir)}`.",
        "",
        "`H0` is the reference Newtonian maximum depth associated with the chosen endpoint.",
    ]
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Generate Rindler brachistochrone publication figures.")
    parser.add_argument("--outdir", default=os.path.join("RindlerMetric", "figures"), help="Output figure directory.")
    parser.add_argument("--datadir", default=None, help="Output CSV directory. Defaults to OUTDIR/data.")
    parser.add_argument("--rho", type=float, default=0.10, help="Reference weak-field parameter.")
    parser.add_argument(
        "--rho-list",
        default="0.02,0.05,0.10,0.20",
        help="Comma-separated rho values for the separation figure.",
    )
    parser.add_argument(
        "--theta-end-deg",
        type=float,
        default=270.0,
        help="Newtonian reference endpoint angle for the fixed-endpoint comparison.",
    )
    parser.add_argument("--samples", type=int, default=3000, help="Number of samples per curve.")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    outdir = args.outdir
    datadir = args.datadir or os.path.join(outdir, "data")
    rhos = parse_rho_list(args.rho_list)
    ensure_dir(outdir)
    ensure_dir(datadir)

    figure_exact_comparison(outdir, datadir, args.rho, args.theta_end_deg, args.samples)
    figure_separation(outdir, datadir, args.rho, args.theta_end_deg, rhos, args.samples)
    figure_exact_vs_pn(outdir, datadir, args.rho, args.theta_end_deg, "coordinate", args.samples)
    figure_exact_vs_pn(outdir, datadir, args.rho, args.theta_end_deg, "proper", args.samples)
    figure_effective_indices(outdir, datadir, args.rho, args.samples)
    write_readme(outdir, datadir, args.rho, rhos, args.theta_end_deg)

    fmts = ["SVG", "CSV"]
    if HAS_REPORTLAB:
        fmts.insert(0, "PDF")
    if HAS_PIL:
        fmts.insert(0, f"PNG ({PNG_DPI} dpi)")
    print(f"Generated Rindler figures in {outdir} ({', '.join(fmts)}).")
    if not HAS_REPORTLAB:
        print("PDF export skipped: reportlab not installed.")
    if not HAS_PIL:
        print("PNG export skipped: Pillow not installed.")


if __name__ == "__main__":
    main()
