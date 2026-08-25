#!/usr/bin/env python3
"""
Generate APS-style publication figures for TOV interior brachistochrones.

Outputs are vector SVG and PDF files plus CSV source data.  The script does not
require matplotlib; it uses numpy for numerics and reportlab for PDF export.
"""

from __future__ import annotations

import argparse
import csv
import html
import math
import os
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

from tov_brachistochrone import (
    PolytropicEOS,
    TabulatedEOS,
    TOVSolution,
    aperture_scan,
    compactness_sequence,
    effective_index,
    full_curve_xy,
    integrate_tov,
    local_speed,
)

try:
    from reportlab.lib import colors
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.ttfonts import TTFont
    from reportlab.pdfgen import canvas

    HAS_REPORTLAB = True
except Exception:
    HAS_REPORTLAB = False


Color = str


PALETTE: Dict[str, Color] = {
    "black": "#111111",
    "gray": "#666666",
    "light_gray": "#d9d9d9",
    "pressure": "#0072B2",
    "density": "#009E73",
    "mass": "#333333",
    "lapse": "#CC79A7",
    "metric": "#E69F00",
    "coordinate": "#0072B2",
    "proper": "#D55E00",
    "surface": "#8a8a8a",
    "accent": "#009E73",
}
M_SUN_KM = 1.4766250614046494
GREEK_FONT_NAME = "Helvetica"
GREEK_FONT_PATH = "/System/Library/Fonts/Supplemental/STIXTwoText.ttf"

if HAS_REPORTLAB and os.path.exists(GREEK_FONT_PATH):
    pdfmetrics.registerFont(TTFont("STIXTwoText", GREEK_FONT_PATH))
    GREEK_FONT_NAME = "STIXTwoText"


@dataclass
class Series:
    x: np.ndarray
    y: np.ndarray
    label: str
    color: Color = "#111111"
    width: float = 1.45
    style: str = "solid"


@dataclass
class Inset:
    rect: Tuple[float, float, float, float]
    xlim: Tuple[float, float]
    ylim: Tuple[float, float]
    label: str
    series: List[Series] = field(default_factory=list)
    show_yticks: bool = False
    connect_zoom: bool = False


@dataclass
class Panel:
    rect: Tuple[float, float, float, float]
    xlim: Tuple[float, float]
    ylim: Tuple[float, float]
    xlabel: str
    ylabel: str
    series: List[Series] = field(default_factory=list)
    tag: str = ""
    tag_loc: str = "upper left"
    legend: bool = True
    legend_location: str = "upper right"
    inset: Optional[Inset] = None
    equal_aspect: bool = False
    xticks: Optional[Sequence[float]] = None
    yticks: Optional[Sequence[float]] = None


def parse_float_list(text: str) -> List[float]:
    values = []
    for raw in text.split(","):
        item = raw.strip()
        if item:
            values.append(float(item))
    if not values:
        raise ValueError("empty numeric list")
    return values


def default_table_path() -> str:
    script_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(script_dir, "eos_tables", "sample_polytrope_table.csv")


def default_compose_dir() -> str:
    script_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(script_dir, "eos_tables", "compose", "qmc_rmf_1")


def build_eos(args: argparse.Namespace) -> PolytropicEOS | TabulatedEOS:
    if args.eos == "polytrope":
        return PolytropicEOS(kappa=args.kappa, gamma=args.gamma)
    if args.eos == "compose":
        compose_dir = args.compose_dir or default_compose_dir()
        nb_path = args.compose_nb or os.path.join(compose_dir, "eos.nb.ns")
        thermo_path = args.compose_thermo or os.path.join(compose_dir, "eos.thermo.ns")
        return TabulatedEOS.from_compose_ns(nb_path, thermo_path, name=args.eos_name)
    table_path = args.eos_table or default_table_path()
    return TabulatedEOS.from_csv(
        table_path,
        density_column=args.density_column,
        pressure_column=args.pressure_column,
        epsilon_column=args.epsilon_column,
        name=args.eos_name,
    )


def central_pressure_from_args(args: argparse.Namespace, eos: PolytropicEOS | TabulatedEOS) -> float:
    if args.central_density is not None:
        return float(eos.pressure_from_density(args.central_density))
    return float(args.central_pressure)


def sequence_pressures_from_args(args: argparse.Namespace, eos: PolytropicEOS | TabulatedEOS) -> List[float]:
    if args.sequence_density:
        return [float(eos.pressure_from_density(value)) for value in parse_float_list(args.sequence_density)]
    return parse_float_list(args.sequence_pc)


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def save_csv(path: str, columns: Sequence[Tuple[str, np.ndarray]]) -> None:
    ensure_dir(os.path.dirname(path))
    length = len(columns[0][1])
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow([name for name, _data in columns])
        for row_index in range(length):
            row = []
            for _name, data in columns:
                value = data[row_index]
                if isinstance(value, str):
                    row.append(value)
                else:
                    row.append(f"{float(value):.12g}")
            writer.writerow(row)


def safe_name(value: float) -> str:
    return f"{value:.5g}".replace("-", "m").replace(".", "p")


def nice_ticks(lo: float, hi: float, target: int = 5) -> List[float]:
    if not math.isfinite(lo) or not math.isfinite(hi) or hi <= lo:
        return [lo, hi]
    span = hi - lo
    raw_step = span / max(target - 1, 1)
    exponent = math.floor(math.log10(raw_step))
    fraction = raw_step / 10.0**exponent
    if fraction <= 1.5:
        step = 1.0 * 10.0**exponent
    elif fraction <= 3.0:
        step = 2.0 * 10.0**exponent
    elif fraction <= 7.0:
        step = 5.0 * 10.0**exponent
    else:
        step = 10.0 * 10.0**exponent
    start = math.ceil(lo / step) * step
    ticks = []
    value = start
    eps = max(abs(lo), abs(hi), abs(step), 1.0) * 1.0e-10
    while value <= hi + eps:
        if value >= lo - eps:
            ticks.append(0.0 if abs(value) < 1.0e-13 else value)
        value += step
    return ticks if ticks else [lo, hi]


def fmt_tick(value: float) -> str:
    if abs(value) < 1.0e-12:
        return "0"
    if abs(value) >= 1000.0 or abs(value) < 0.001:
        return f"{value:.1e}"
    if abs(value - round(value)) < 1.0e-10:
        return str(int(round(value)))
    return f"{value:.3f}".rstrip("0").rstrip(".")


def dash_pattern(style: str) -> Optional[List[float]]:
    if style == "dash":
        return [7.0, 4.0]
    if style == "dot":
        return [1.5, 3.0]
    if style == "dashdot":
        return [7.0, 3.5, 1.5, 3.5]
    return None


def hex_to_rgb01(value: str) -> Tuple[float, float, float]:
    stripped = value.strip().lstrip("#")
    return tuple(int(stripped[i : i + 2], 16) / 255.0 for i in (0, 2, 4))


GREEK_SVG = {
    "Delta": "&#916;",
    "epsilon": "&#949;",
    "gamma": "&#947;",
    "Gamma": "&#915;",
    "tau": "&#964;",
    "mu": "&#956;",
}
GREEK_PDF_SYMBOL = {
    "Delta": "Δ",
    "epsilon": "ε",
    "gamma": "γ",
    "Gamma": "Γ",
    "tau": "τ",
    "mu": "μ",
}
GREEK_SYMBOLS = tuple(sorted(GREEK_SVG, key=len, reverse=True))


def svg_text(value: str) -> str:
    parts = []
    for kind, text in split_symbol_text(value):
        if kind == "symbol":
            parts.append(GREEK_SVG[text])
        elif kind == "subscript_symbol":
            parts.append(f'<tspan baseline-shift="sub" font-size="70%">{GREEK_SVG[text]}</tspan>')
        elif kind == "subscript_text":
            parts.append(f'<tspan baseline-shift="sub" font-size="70%">{html.escape(text)}</tspan>')
        else:
            parts.append(html.escape(text))
    return "".join(parts)


def split_symbol_text(value: str) -> List[Tuple[str, str]]:
    parts: List[Tuple[str, str]] = []
    index = 0
    while index < len(value):
        if value[index] == "\\":
            matched = False
            for name in GREEK_SYMBOLS:
                token = f"\\{name}"
                if value.startswith(token, index):
                    parts.append(("symbol", name))
                    index += len(token)
                    matched = True
                    break
            if not matched:
                parts.append(("text", value[index]))
                index += 1
            continue
        if value[index] == "_" and index + 1 < len(value):
            next_index = index + 1
            if value[next_index] == "{":
                end_index = value.find("}", next_index + 1)
                if end_index != -1:
                    subscript = value[next_index + 1 : end_index]
                    if subscript.startswith("\\") and subscript[1:] in GREEK_SVG:
                        parts.append(("subscript_symbol", subscript[1:]))
                    else:
                        parts.append(("subscript_text", subscript))
                    index = end_index + 1
                    continue
            matched = False
            for name in GREEK_SYMBOLS:
                token = f"\\{name}"
                if value.startswith(token, next_index):
                    parts.append(("subscript_symbol", name))
                    index = next_index + len(token)
                    matched = True
                    break
            if matched:
                continue
            parts.append(("subscript_text", value[next_index]))
            index = next_index + 1
            continue
        next_positions = [position for position in (value.find("\\", index), value.find("_", index)) if position != -1]
        next_index = min(next_positions) if next_positions else len(value)
        parts.append(("text", value[index:next_index]))
        index = next_index
    return parts


def pdf_symbol_text_width(value: str, font_name: str, font_size: float) -> float:
    width = 0.0
    for kind, text in split_symbol_text(value):
        if kind == "symbol":
            width += pdfmetrics.stringWidth(GREEK_PDF_SYMBOL[text], GREEK_FONT_NAME, font_size)
        elif kind == "subscript_symbol":
            width += pdfmetrics.stringWidth(GREEK_PDF_SYMBOL[text], GREEK_FONT_NAME, 0.72 * font_size)
        elif kind == "subscript_text":
            width += pdfmetrics.stringWidth(text, font_name, 0.72 * font_size)
        else:
            width += pdfmetrics.stringWidth(text, font_name, font_size)
    return width


def draw_pdf_greek_symbol(doc: canvas.Canvas, x: float, y: float, name: str, font_size: float) -> float:
    doc.setStrokeColor(colors.HexColor("#111111"))
    doc.setLineWidth(max(0.45, 0.065 * font_size))
    if name == "Delta":
        width = 0.62 * font_size
        height = 0.78 * font_size
        path = doc.beginPath()
        path.moveTo(x + 0.04 * font_size, y)
        path.lineTo(x + 0.50 * width, y + height)
        path.lineTo(x + width - 0.04 * font_size, y)
        path.close()
        doc.drawPath(path, stroke=1, fill=0)
        return width
    if name == "tau":
        width = 0.66 * font_size
        height = 0.78 * font_size
        doc.setLineWidth(max(0.65, 0.09 * font_size))
        doc.line(x + 0.04 * font_size, y + 0.66 * height, x + width, y + 0.66 * height)
        path = doc.beginPath()
        path.moveTo(x + 0.52 * width, y + 0.66 * height)
        path.curveTo(x + 0.44 * width, y + 0.44 * height, x + 0.43 * width, y + 0.16 * height, x + 0.72 * width, y)
        path.curveTo(x + 0.82 * width, y - 0.02 * height, x + 0.92 * width, y + 0.02 * height, x + 0.98 * width, y + 0.08 * height)
        doc.drawPath(path, stroke=1, fill=0)
        return width
    if name == "mu":
        width = 0.58 * font_size
        height = 0.78 * font_size
        path = doc.beginPath()
        path.moveTo(x + 0.08 * width, y - 0.18 * height)
        path.lineTo(x + 0.08 * width, y + 0.66 * height)
        path.moveTo(x + 0.08 * width, y + 0.15 * height)
        path.curveTo(x + 0.22 * width, y - 0.05 * height, x + 0.42 * width, y - 0.02 * height, x + 0.50 * width, y + 0.18 * height)
        path.lineTo(x + 0.50 * width, y + 0.66 * height)
        path.moveTo(x + 0.50 * width, y + 0.18 * height)
        path.curveTo(x + 0.60 * width, y - 0.02 * height, x + 0.78 * width, y - 0.02 * height, x + 0.88 * width, y + 0.12 * height)
        doc.drawPath(path, stroke=1, fill=0)
        return width
    raise ValueError(f"Unsupported Greek symbol: {name}")


def draw_pdf_symbol_text(
    doc: canvas.Canvas,
    x: float,
    y: float,
    value: str,
    font_name: str,
    font_size: float,
    anchor: str = "left",
) -> None:
    width = pdf_symbol_text_width(value, font_name, font_size)
    cursor = x
    doc.setFillColor(colors.HexColor("#111111"))
    doc.setStrokeColor(colors.HexColor("#111111"))
    if anchor == "center":
        cursor -= 0.5 * width
    elif anchor == "right":
        cursor -= width
    for kind, text in split_symbol_text(value):
        if kind == "symbol":
            doc.setFont(GREEK_FONT_NAME, font_size)
            doc.setFillColor(colors.HexColor("#111111"))
            symbol_text = GREEK_PDF_SYMBOL[text]
            doc.drawString(cursor, y, symbol_text)
            cursor += pdfmetrics.stringWidth(symbol_text, GREEK_FONT_NAME, font_size)
        elif kind == "subscript_symbol":
            subscript_size = 0.72 * font_size
            subscript_y = y - 0.24 * font_size
            doc.setFont(GREEK_FONT_NAME, subscript_size)
            doc.setFillColor(colors.HexColor("#111111"))
            symbol_text = GREEK_PDF_SYMBOL[text]
            doc.drawString(cursor, subscript_y, symbol_text)
            cursor += pdfmetrics.stringWidth(symbol_text, GREEK_FONT_NAME, subscript_size)
        elif kind == "subscript_text":
            subscript_size = 0.72 * font_size
            subscript_y = y - 0.24 * font_size
            doc.setFont(font_name, subscript_size)
            doc.setFillColor(colors.HexColor("#111111"))
            doc.drawString(cursor, subscript_y, text)
            cursor += pdfmetrics.stringWidth(text, font_name, subscript_size)
        else:
            doc.setFont(font_name, font_size)
            doc.setFillColor(colors.HexColor("#111111"))
            doc.drawString(cursor, y, text)
            cursor += pdfmetrics.stringWidth(text, font_name, font_size)


def limits_from_series(
    series: Sequence[Series],
    pad_x: float = 0.04,
    pad_y: float = 0.08,
) -> Tuple[Tuple[float, float], Tuple[float, float]]:
    x_values = np.concatenate([item.x[np.isfinite(item.x)] for item in series])
    y_values = np.concatenate([item.y[np.isfinite(item.y)] for item in series])
    xmin = float(np.min(x_values))
    xmax = float(np.max(x_values))
    ymin = float(np.min(y_values))
    ymax = float(np.max(y_values))
    xrange = max(xmax - xmin, 1.0e-12)
    yrange = max(ymax - ymin, 1.0e-12)
    return (
        (xmin - pad_x * xrange, xmax + pad_x * xrange),
        (ymin - pad_y * yrange, ymax + pad_y * yrange),
    )


def adjusted_limits(panel: Panel, pixel_width: float, pixel_height: float) -> Tuple[Tuple[float, float], Tuple[float, float]]:
    xlim = panel.xlim
    ylim = panel.ylim
    if not panel.equal_aspect:
        return xlim, ylim
    xrange = xlim[1] - xlim[0]
    yrange = ylim[1] - ylim[0]
    if xrange <= 0.0 or yrange <= 0.0:
        return xlim, ylim
    data_ratio = xrange / yrange
    box_ratio = pixel_width / pixel_height
    if data_ratio > box_ratio:
        new_yrange = xrange / box_ratio
        center = 0.5 * (ylim[0] + ylim[1])
        ylim = (center - 0.5 * new_yrange, center + 0.5 * new_yrange)
    else:
        new_xrange = yrange * box_ratio
        center = 0.5 * (xlim[0] + xlim[1])
        xlim = (center - 0.5 * new_xrange, center + 0.5 * new_xrange)
    return xlim, ylim


def svg_points(xs: np.ndarray, ys: np.ndarray) -> str:
    return " ".join(f"{x:.2f},{y:.2f}" for x, y in zip(xs, ys) if math.isfinite(x) and math.isfinite(y))


def write_svg(path: str, panels: Sequence[Panel], size: Tuple[int, int]) -> None:
    width, height = size
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        "<style>",
        "text{font-family:Helvetica,Arial,sans-serif;fill:#111}",
        ".tick{font-size:12px}.label{font-size:14px}.tag{font-size:15px;font-weight:700}",
        "</style>",
        f'<rect x="0" y="0" width="{width}" height="{height}" fill="white"/>',
    ]

    for panel_index, panel in enumerate(panels):
        left = panel.rect[0] * width
        top = panel.rect[1] * height
        panel_width = panel.rect[2] * width
        panel_height = panel.rect[3] * height
        xlim, ylim = adjusted_limits(panel, panel_width, panel_height)
        xticks = list(panel.xticks) if panel.xticks is not None else nice_ticks(*xlim)
        yticks = list(panel.yticks) if panel.yticks is not None else nice_ticks(*ylim)
        clip_id = f"clip_tov_{panel_index}"

        def tx(values: np.ndarray) -> np.ndarray:
            return left + (values - xlim[0]) / (xlim[1] - xlim[0]) * panel_width

        def ty(values: np.ndarray) -> np.ndarray:
            return top + panel_height - (values - ylim[0]) / (ylim[1] - ylim[0]) * panel_height

        parts.append(
            f'<clipPath id="{clip_id}"><rect x="{left:.2f}" y="{top:.2f}" '
            f'width="{panel_width:.2f}" height="{panel_height:.2f}"/></clipPath>'
        )

        for tick in xticks:
            x = float(tx(np.array([tick]))[0])
            parts.append(
                f'<line x1="{x:.2f}" y1="{top:.2f}" x2="{x:.2f}" y2="{top + panel_height:.2f}" '
                f'stroke="#ececec" stroke-width="1"/>'
            )
            parts.append(
                f'<line x1="{x:.2f}" y1="{top + panel_height:.2f}" x2="{x:.2f}" '
                f'y2="{top + panel_height + 4:.2f}" stroke="#111" stroke-width="1"/>'
            )
            parts.append(
                f'<text class="tick" x="{x:.2f}" y="{top + panel_height + 18:.2f}" '
                f'text-anchor="middle">{svg_text(fmt_tick(tick))}</text>'
            )

        for tick in yticks:
            y = float(ty(np.array([tick]))[0])
            parts.append(
                f'<line x1="{left:.2f}" y1="{y:.2f}" x2="{left + panel_width:.2f}" y2="{y:.2f}" '
                f'stroke="#ececec" stroke-width="1"/>'
            )
            parts.append(
                f'<line x1="{left - 4:.2f}" y1="{y:.2f}" x2="{left:.2f}" y2="{y:.2f}" '
                f'stroke="#111" stroke-width="1"/>'
            )
            parts.append(
                f'<text class="tick" x="{left - 8:.2f}" y="{y + 4:.2f}" text-anchor="end">'
                f"{svg_text(fmt_tick(tick))}</text>"
            )

        parts.append(
            f'<rect x="{left:.2f}" y="{top:.2f}" width="{panel_width:.2f}" height="{panel_height:.2f}" '
            f'fill="none" stroke="#111" stroke-width="1.1"/>'
        )

        for item in panel.series:
            xs = tx(item.x)
            ys = ty(item.y)
            dash = dash_pattern(item.style)
            dash_attr = "" if dash is None else f' stroke-dasharray="{",".join(str(d) for d in dash)}"'
            parts.append(
                f'<polyline points="{svg_points(xs, ys)}" fill="none" stroke="{item.color}" '
                f'stroke-width="{item.width:.2f}" stroke-linejoin="round" stroke-linecap="round" '
                f'{dash_attr} clip-path="url(#{clip_id})"/>'
            )

        parts.append(
            f'<text class="label" x="{left + 0.5 * panel_width:.2f}" '
            f'y="{top + panel_height + 41:.2f}" text-anchor="middle">{svg_text(panel.xlabel)}</text>'
        )
        ylabel_x = left - 51.0
        ylabel_y = top + 0.5 * panel_height
        parts.append(
            f'<text class="label" x="{ylabel_x:.2f}" y="{ylabel_y:.2f}" text-anchor="middle" '
            f'transform="rotate(-90 {ylabel_x:.2f} {ylabel_y:.2f})">{svg_text(panel.ylabel)}</text>'
        )

        if panel.legend and panel.series:
            max_label = max(len(item.label) for item in panel.series)
            legend_width = max(112.0, 46.0 + 6.4 * max_label)
            legend_height = 18.0 + 18.0 * len(panel.series)
            if panel.legend_location == "lower left":
                legend_x = left + 9.0
                legend_y = top + panel_height - legend_height - 9.0
            else:
                legend_x = left + panel_width - legend_width - 9.0
                legend_y = top + 9.0
            parts.append(
                f'<rect x="{legend_x:.2f}" y="{legend_y:.2f}" width="{legend_width:.2f}" '
                f'height="{legend_height:.2f}" fill="white" stroke="#cccccc" stroke-width="1" opacity="0.94"/>'
            )
            for row, item in enumerate(panel.series):
                y0 = legend_y + 15.0 + 18.0 * row
                dash = dash_pattern(item.style)
                dash_attr = "" if dash is None else f' stroke-dasharray="{",".join(str(d) for d in dash)}"'
                parts.append(
                    f'<line x1="{legend_x + 9:.2f}" y1="{y0:.2f}" x2="{legend_x + 36:.2f}" y2="{y0:.2f}" '
                    f'stroke="{item.color}" stroke-width="{item.width:.2f}" stroke-linecap="round"{dash_attr}/>'
                )
                parts.append(
                    f'<text class="tick" x="{legend_x + 43:.2f}" y="{y0 + 4:.2f}">'
                    f"{svg_text(item.label)}</text>"
                )

        if panel.inset:
            inset = panel.inset
            inset_left = left + inset.rect[0] * panel_width
            inset_top = top + inset.rect[1] * panel_height
            inset_width = inset.rect[2] * panel_width
            inset_height = inset.rect[3] * panel_height

            def itx(values: np.ndarray) -> np.ndarray:
                return inset_left + (values - inset.xlim[0]) / (inset.xlim[1] - inset.xlim[0]) * inset_width

            def ity(values: np.ndarray) -> np.ndarray:
                return inset_top + inset_height - (values - inset.ylim[0]) / (inset.ylim[1] - inset.ylim[0]) * inset_height

            parts.append(
                f'<rect x="{inset_left:.2f}" y="{inset_top:.2f}" width="{inset_width:.2f}" '
                f'height="{inset_height:.2f}" fill="white" stroke="#bbbbbb" stroke-width="1" opacity="0.94"/>'
            )
            for item in inset.series:
                xs = itx(item.x)
                ys = ity(item.y)
                dash = dash_pattern(item.style)
                dash_attr = "" if dash is None else f' stroke-dasharray="{",".join(str(d) for d in dash)}"'
                parts.append(
                    f'<polyline points="{svg_points(xs, ys)}" fill="none" stroke="{item.color}" '
                    f'stroke-width="{item.width:.2f}" stroke-linejoin="round" stroke-linecap="round"{dash_attr}/>'
                )
            parts.append(
                f'<text class="tick" x="{inset_left + 5:.2f}" y="{inset_top + 13:.2f}">'
                f"{svg_text(inset.label)}</text>"
            )

        if panel.tag:
            tag_width = max(76.0, 9.5 + 7.3 * len(panel.tag))
            tag_y = top + 18.0
            if panel.tag_loc == "upper right":
                tag_x = left + panel_width - tag_width - 3.0
            else:
                tag_x = left + 7.0
            parts.append(
                f'<rect x="{tag_x - 4:.2f}" y="{tag_y - 15:.2f}" width="{tag_width:.2f}" height="20.00" '
                f'fill="white" opacity="0.88"/>'
            )
            parts.append(
                f'<text class="tag" x="{tag_x:.2f}" y="{tag_y:.2f}">{svg_text(panel.tag)}</text>'
            )

    parts.append("</svg>")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(parts))


def draw_pdf_polyline(doc: canvas.Canvas, xs: Sequence[float], ys: Sequence[float]) -> None:
    points = [(float(x), float(y)) for x, y in zip(xs, ys) if math.isfinite(x) and math.isfinite(y)]
    if len(points) < 2:
        return
    path = doc.beginPath()
    path.moveTo(points[0][0], points[0][1])
    for x, y in points[1:]:
        path.lineTo(x, y)
    doc.drawPath(path, stroke=1, fill=0)


_MPL_GREEK: Dict[str, str] = {
    "Delta": "Δ", "Gamma": "Γ", "Lambda": "Λ", "Omega": "Ω", "Phi": "Φ",
    "Pi": "Π", "Psi": "Ψ", "Sigma": "Σ", "Theta": "Θ", "Xi": "Ξ",
    "tau": "τ", "mu": "μ", "gamma": "γ", "alpha": "α",
    "beta": "β", "sigma": "σ", "rho": "ρ", "phi": "φ", "theta": "θ",
    "pi": "π", "lambda": "λ", "epsilon": "ε", "eta": "η", "nu": "ν",
    "omega": "ω", "xi": "ξ", "kappa": "κ", "delta": "δ",
}


def _mpl_label(s: str) -> str:
    import re
    # pass mathtext strings through unchanged
    if s.startswith("$") and s.endswith("$"):
        return s
    result = s
    for name, char in _MPL_GREEK.items():
        result = result.replace("\\" + name, char)
    result = re.sub(r"_\{([^}]*)\}", r"_\1", result)
    return result


def write_png(path: str, panels: Sequence[Panel], size_inches: Tuple[float, float], dpi: int = 300) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    _ls = {"solid": "-", "dash": "--", "dot": ":", "dashdot": "-."}
    fig = plt.figure(figsize=size_inches, facecolor="white", dpi=dpi)
    for panel in panels:
        bot_f = 1.0 - panel.rect[1] - panel.rect[3]
        ax = fig.add_axes([panel.rect[0], bot_f, panel.rect[2], panel.rect[3]])
        ax.set_xlim(panel.xlim)
        ax.set_ylim(panel.ylim)
        if panel.equal_aspect:
            ax.set_aspect("equal", adjustable="box")
        if panel.xticks is not None:
            ax.set_xticks(list(panel.xticks))
        if panel.yticks is not None:
            ax.set_yticks(list(panel.yticks))
        ax.tick_params(axis="both", labelsize=8.0)
        ax.set_xlabel(_mpl_label(panel.xlabel), fontsize=9.0)
        ax.set_ylabel(_mpl_label(panel.ylabel), fontsize=9.0)
        ax.grid(True, color="#ececec", linewidth=0.5, zorder=0)
        for item in panel.series:
            ax.plot(item.x, item.y, color=item.color, linewidth=item.width * 0.65,
                    linestyle=_ls.get(item.style, "-"), label=_mpl_label(item.label))
        if panel.legend and panel.series:
            _loc_map = {"lower left": "lower left", "upper left": "upper left",
                        "lower right": "lower right"}
            loc = _loc_map.get(panel.legend_location, "upper right")
            ax.legend(loc=loc, fontsize=8.0, framealpha=0.94)
        if panel.tag:
            _tx = 0.97 if panel.tag_loc == "upper right" else 0.03
            _ha = "right" if panel.tag_loc == "upper right" else "left"
            ax.text(_tx, 0.97, _mpl_label(panel.tag), transform=ax.transAxes, fontsize=9.0,
                    fontweight="bold", va="top", ha=_ha,
                    bbox=dict(facecolor="white", edgecolor="none", alpha=0.88, pad=2))
        if panel.inset is not None:
            ins = panel.inset
            # Convert axes-normalised inset rect → figure coordinates
            pan_l = panel.rect[0]
            pan_b = 1.0 - panel.rect[1] - panel.rect[3]
            pan_w = panel.rect[2]
            pan_h = panel.rect[3]
            ins_ax_l = ins.rect[0]
            ins_ax_b = 1.0 - ins.rect[1] - ins.rect[3]
            fig_l = pan_l + ins_ax_l * pan_w
            fig_b = pan_b + ins_ax_b * pan_h
            fig_w = ins.rect[2] * pan_w
            fig_h = ins.rect[3] * pan_h
            iax = fig.add_axes([fig_l, fig_b, fig_w, fig_h])
            iax.set_facecolor("white")
            for sp in iax.spines.values():
                sp.set_linewidth(0.8)
                sp.set_edgecolor("0.40")
            iax.set_xlim(ins.xlim)
            iax.set_ylim(ins.ylim)
            iax.tick_params(labelsize=6.5)
            if not ins.show_yticks:
                iax.set_yticklabels([])
            for item in ins.series:
                iax.plot(item.x, item.y, color=item.color, linewidth=item.width * 0.65,
                         linestyle=_ls.get(item.style, "-"))
            iax.text(0.05, 0.93, _mpl_label(ins.label), transform=iax.transAxes,
                     fontsize=7.0, va="top")
            if ins.connect_zoom:
                from matplotlib.patches import Rectangle
                ax.add_patch(Rectangle(
                    (ins.xlim[0], ins.ylim[0]),
                    ins.xlim[1] - ins.xlim[0],
                    ins.ylim[1] - ins.ylim[0],
                    linewidth=0.75, edgecolor="0.50", facecolor="none", alpha=0.85,
                    zorder=4,
                ))
    fig.savefig(path, dpi=dpi, facecolor="white", bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)


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
        top = panel.rect[1] * height
        panel_width = panel.rect[2] * width
        panel_height = panel.rect[3] * height
        bottom = height - top - panel_height
        xlim, ylim = adjusted_limits(panel, panel_width, panel_height)
        xticks = list(panel.xticks) if panel.xticks is not None else nice_ticks(*xlim)
        yticks = list(panel.yticks) if panel.yticks is not None else nice_ticks(*ylim)

        def tx(values: np.ndarray) -> np.ndarray:
            return left + (values - xlim[0]) / (xlim[1] - xlim[0]) * panel_width

        def ty(values: np.ndarray) -> np.ndarray:
            return bottom + (values - ylim[0]) / (ylim[1] - ylim[0]) * panel_height

        doc.setStrokeColor(colors.HexColor("#ececec"))
        doc.setLineWidth(0.45)
        for tick in xticks:
            x = float(tx(np.array([tick]))[0])
            doc.line(x, bottom, x, bottom + panel_height)
        for tick in yticks:
            y = float(ty(np.array([tick]))[0])
            doc.line(left, y, left + panel_width, y)

        doc.setStrokeColor(colors.HexColor("#111111"))
        doc.setLineWidth(0.75)
        doc.rect(left, bottom, panel_width, panel_height, stroke=1, fill=0)
        doc.setFont("Helvetica", 8.0)
        for tick in xticks:
            x = float(tx(np.array([tick]))[0])
            doc.line(x, bottom, x, bottom - 2.5)
            draw_pdf_symbol_text(doc, x, bottom - 13.0, fmt_tick(tick), "Helvetica", 8.0, "center")
        for tick in yticks:
            y = float(ty(np.array([tick]))[0])
            doc.line(left - 2.5, y, left, y)
            draw_pdf_symbol_text(doc, left - 5.5, y - 3.0, fmt_tick(tick), "Helvetica", 8.0, "right")

        doc.setFont("Helvetica", 9.0)
        draw_pdf_symbol_text(doc, left + 0.5 * panel_width, bottom - 29.0, panel.xlabel, "Helvetica", 9.0, "center")
        doc.saveState()
        doc.translate(left - 40.0, bottom + 0.5 * panel_height)
        doc.rotate(90)
        draw_pdf_symbol_text(doc, 0.0, 0.0, panel.ylabel, "Helvetica", 9.0, "center")
        doc.restoreState()
        for item in panel.series:
            doc.setStrokeColor(colors.Color(*hex_to_rgb01(item.color)))
            doc.setLineWidth(item.width * 0.62)
            pattern = dash_pattern(item.style)
            if pattern:
                doc.setDash(pattern)
            else:
                doc.setDash()
            draw_pdf_polyline(doc, tx(item.x), ty(item.y))
        doc.setDash()

        if panel.legend and panel.series:
            max_label = max(len(item.label) for item in panel.series)
            legend_width = max(90.0, 32.0 + 4.8 * max_label)
            legend_height = 14.0 + 12.5 * len(panel.series)
            if panel.legend_location == "lower left":
                legend_x = left + 6.0
                legend_y = bottom + 6.0
            else:
                legend_x = left + panel_width - legend_width - 6.0
                legend_y = bottom + panel_height - legend_height - 6.0
            doc.setFillColor(colors.white)
            doc.setStrokeColor(colors.HexColor("#cccccc"))
            doc.setLineWidth(0.5)
            doc.rect(legend_x, legend_y, legend_width, legend_height, stroke=1, fill=1)
            doc.setFont("Helvetica", 8.0)
            doc.setFillColor(colors.HexColor("#111111"))
            for row, item in enumerate(panel.series):
                y0 = legend_y + legend_height - 10.0 - 12.5 * row
                doc.setStrokeColor(colors.Color(*hex_to_rgb01(item.color)))
                doc.setLineWidth(item.width * 0.62)
                pattern = dash_pattern(item.style)
                if pattern:
                    doc.setDash(pattern)
                else:
                    doc.setDash()
                doc.line(legend_x + 6.0, y0, legend_x + 24.0, y0)
                doc.setDash()
                doc.setFillColor(colors.HexColor("#111111"))
                draw_pdf_symbol_text(doc, legend_x + 29.0, y0 - 3.2, item.label, "Helvetica", 8.0)

        if panel.inset:
            inset = panel.inset
            inset_left = left + inset.rect[0] * panel_width
            inset_top = bottom + panel_height - (inset.rect[1] + inset.rect[3]) * panel_height
            inset_width = inset.rect[2] * panel_width
            inset_height = inset.rect[3] * panel_height

            def itx(values: np.ndarray) -> np.ndarray:
                return inset_left + (values - inset.xlim[0]) / (inset.xlim[1] - inset.xlim[0]) * inset_width

            def ity(values: np.ndarray) -> np.ndarray:
                return inset_top + (values - inset.ylim[0]) / (inset.ylim[1] - inset.ylim[0]) * inset_height

            doc.setFillColor(colors.white)
            doc.setStrokeColor(colors.HexColor("#bbbbbb"))
            doc.setLineWidth(0.45)
            doc.rect(inset_left, inset_top, inset_width, inset_height, stroke=1, fill=1)
            for item in inset.series:
                doc.setStrokeColor(colors.Color(*hex_to_rgb01(item.color)))
                doc.setLineWidth(item.width * 0.56)
                pattern = dash_pattern(item.style)
                if pattern:
                    doc.setDash(pattern)
                else:
                    doc.setDash()
                draw_pdf_polyline(doc, itx(item.x), ity(item.y))
            doc.setDash()
            doc.setFillColor(colors.HexColor("#111111"))
            doc.setFont("Helvetica", 6.4)
            draw_pdf_symbol_text(doc, inset_left + 3.0, inset_top + inset_height - 8.0, inset.label, "Helvetica", 6.4)

        if panel.tag:
            doc.setFont("Helvetica-Bold", 9.2)
            tag_x = left + 6.0
            tag_y = bottom + panel_height - 13.0
            tag_width = max(58.0, pdf_symbol_text_width(panel.tag, "Helvetica-Bold", 9.2) + 8.0)
            doc.setFillColor(colors.white)
            doc.rect(tag_x - 3.0, tag_y - 3.0, tag_width, 12.0, stroke=0, fill=1)
            doc.setFillColor(colors.HexColor("#111111"))
            draw_pdf_symbol_text(doc, tag_x, tag_y, panel.tag, "Helvetica-Bold", 9.2)

    doc.showPage()
    doc.save()


def save_figure(
    outdir: str,
    basename: str,
    panels: Sequence[Panel],
    svg_size: Tuple[int, int],
    pdf_size_inches: Tuple[float, float],
) -> None:
    ensure_dir(outdir)
    write_svg(os.path.join(outdir, f"{basename}.svg"), panels, svg_size)
    write_png(os.path.join(outdir, f"{basename}.png"), panels, pdf_size_inches)
    if HAS_REPORTLAB:
        write_pdf(os.path.join(outdir, f"{basename}.pdf"), panels, (pdf_size_inches[0] * 72.0, pdf_size_inches[1] * 72.0))


def panel_grid_2x2() -> List[Tuple[float, float, float, float]]:
    return [
        (0.095, 0.095, 0.375, 0.34),
        (0.585, 0.095, 0.355, 0.34),
        (0.095, 0.575, 0.375, 0.32),
        (0.585, 0.575, 0.355, 0.32),
    ]


def make_profiles(sol: TOVSolution, outdir: str, datadir: str) -> None:
    radius = sol.sample(1400)
    x = radius / sol.radius
    pressure = np.asarray(sol.interp(sol.pressure, radius), dtype=float)
    epsilon = np.asarray(sol.interp(sol.epsilon, radius), dtype=float)
    mass = np.asarray(sol.interp(sol.mass, radius), dtype=float)
    lapse = np.asarray(sol.interp(sol.lapse, radius), dtype=float)
    radial_metric = np.asarray(sol.interp(sol.radial_metric, radius), dtype=float)
    speed, gamma = local_speed(sol, radius)

    rects = panel_grid_2x2()
    panels = [
        Panel(
            rects[0],
            (0.0, 1.0),
            (-0.03, 1.08),
            "r / R",
            "normalized profile",
            [
                Series(x, pressure / sol.central_pressure, "p / p_c", PALETTE["pressure"], 1.55),
                Series(x, epsilon / epsilon[0], "\\epsilon / \\epsilon_c", PALETTE["density"], 1.55, "dash"),
            ],
            "(a)",
        ),
        Panel(
            rects[1],
            (0.0, 1.0),
            (-0.03, 1.08),
            "r / R",
            "m(r) / M",
            [Series(x, mass / sol.total_mass, "mass", PALETTE["mass"], 1.65)],
            "(b)",
            legend=False,
        ),
        Panel(
            rects[2],
            (0.0, 1.0),
            (min(float(np.min(lapse / sol.surface_lapse)), 0.75) * 0.96, max(float(np.max(radial_metric)), 1.05) * 1.04),
            "r / R",
            "metric functions",
            [
                Series(x, lapse / sol.surface_lapse, "A / A(R)", PALETTE["lapse"], 1.55),
                Series(x, radial_metric, "B", PALETTE["metric"], 1.55, "dash"),
            ],
            "(c)",
        ),
        Panel(
            rects[3],
            (0.0, 1.0),
            (-0.02, max(float(np.max(speed)), float(np.max(gamma - 1.0))) * 1.13),
            "r / R",
            "release kinematics",
            [
                Series(x, speed, "v / c", PALETTE["coordinate"], 1.55),
                Series(x, gamma - 1.0, "\\gamma - 1", PALETTE["proper"], 1.55, "dash"),
            ],
            "(d)",
        ),
    ]
    base = "fig_01_tov_profiles"
    save_figure(outdir, base, panels, (1250, 820), (7.0, 4.6))
    save_csv(
        os.path.join(datadir, f"{base}.csv"),
        [
            ("r_over_R", x),
            ("p_over_pc", pressure / sol.central_pressure),
            ("eps_over_epsc", epsilon / epsilon[0]),
            ("m_over_M", mass / sol.total_mass),
            ("A_over_AR", lapse / sol.surface_lapse),
            ("B", radial_metric),
            ("v_over_c", speed),
            ("gamma_minus_1", gamma - 1.0),
        ],
    )


def make_indices(sol: TOVSolution, outdir: str, datadir: str) -> None:
    radius = sol.sample(1300, stop=0.995)
    x = radius / sol.radius
    n_coordinate = np.asarray(effective_index(sol, radius, "coordinate"), dtype=float)
    n_proper = np.asarray(effective_index(sol, radius, "proper"), dtype=float)
    ratio = n_coordinate / n_proper
    lapse = np.asarray(sol.interp(sol.lapse, radius), dtype=float)

    left_series = [
        Series(x, np.log10(n_coordinate), r"$n_t$", PALETTE["coordinate"], 1.6),
        Series(x, np.log10(n_proper), r"$n_\tau$", PALETTE["proper"], 1.6, "dash"),
    ]
    right_series = [
        Series(x, ratio, r"$n_t / n_\tau$", PALETTE["black"], 1.6),
        Series(x, 1.0 / lapse, r"$1/A$", PALETTE["accent"], 1.3, "dash"),
    ]
    left_xlim, left_ylim = limits_from_series(left_series, pad_x=0.0, pad_y=0.08)
    right_xlim, right_ylim = limits_from_series(right_series, pad_x=0.0, pad_y=0.08)
    # Pannelli impilati verticalmente — singola colonna (3.375 in)
    # Convenzione rect: y=0 in ALTO (pixel) → pannello superiore ha y0 piccolo
    panels = [
        Panel(
            (0.15, 0.06, 0.79, 0.40),   # pannello superiore (a)
            (0.0, 1.0),
            left_ylim,
            r"$r/R$",
            r"$\log_{10}\,n$",
            left_series,
            "(a)",
        ),
        Panel(
            (0.15, 0.52, 0.79, 0.40),   # pannello inferiore (b)
            (0.0, 1.0),
            right_ylim,
            r"$r/R$",
            r"$n_t / n_\tau$",
            right_series,
            "(b)",
        ),
    ]
    base = "fig_02_effective_indices"
    save_figure(outdir, base, panels, (1013, 1080), (3.375, 3.6))
    save_csv(
        os.path.join(datadir, f"{base}.csv"),
        [
            ("r_over_R", x),
            ("n_coordinate", n_coordinate),
            ("n_proper", n_proper),
            ("log10_n_coordinate", np.log10(n_coordinate)),
            ("log10_n_proper", np.log10(n_proper)),
            ("n_coordinate_over_n_proper", ratio),
            ("one_over_A", 1.0 / lapse),
        ],
    )


def make_curves(sol: TOVSolution, outdir: str, datadir: str, deltas_deg: Sequence[float]) -> None:
    theta = np.linspace(0.0, 2.0 * math.pi, 900)
    surface_x = np.cos(theta)
    surface_y = np.sin(theta)
    panels = []
    all_columns: List[Tuple[str, np.ndarray]] = []

    for index, delta_deg in enumerate(deltas_deg):
        delta = math.radians(float(delta_deg))
        x_coordinate, y_coordinate, q_coordinate, delta_coordinate = full_curve_xy(sol, "coordinate", delta)
        x_proper, y_proper, q_proper, delta_proper = full_curve_xy(sol, "proper", delta)
        chord_y = np.linspace(-math.sin(0.5 * delta), math.sin(0.5 * delta), 220)
        chord_x = np.full_like(chord_y, math.cos(0.5 * delta))
        panel_series = [
            Series(surface_x, surface_y, "surface", PALETTE["surface"], 1.0, "dot"),
        ]
        if index == 1:
            panel_series.append(Series(chord_x, chord_y, "chord", PALETTE["light_gray"], 1.2, "dashdot"))
        panel_series.extend(
            [
                Series(x_coordinate, y_coordinate, "coordinate", PALETTE["coordinate"], 1.55),
                Series(x_proper, y_proper, "proper", PALETTE["proper"], 1.55, "dash"),
            ]
        )
        inset = None
        if index == 1:
            inset_x_values = np.concatenate([x_coordinate, x_proper, chord_x])
            inset_y_values = np.concatenate([y_coordinate, y_proper, chord_y])
            inset_xmin = float(np.min(inset_x_values))
            inset_xmax = float(np.max(inset_x_values))
            inset_ymin = float(np.min(inset_y_values))
            inset_ymax = float(np.max(inset_y_values))
            inset_xpad = 0.08 * (inset_xmax - inset_xmin)
            inset_ypad = 0.04 * (inset_ymax - inset_ymin)
            inset = Inset(
                rect=(0.08, 0.56, 0.48, 0.34),
                xlim=(inset_xmin - inset_xpad, inset_xmax + inset_xpad),
                ylim=(inset_ymin - inset_ypad, inset_ymax + inset_ypad),
                label="zoom",
                series=[
                    Series(chord_x, chord_y, "chord", PALETTE["light_gray"], 1.0, "dashdot"),
                    Series(x_coordinate, y_coordinate, "coordinate", PALETTE["coordinate"], 1.25),
                    Series(x_proper, y_proper, "proper", PALETTE["proper"], 1.25, "dash"),
                ],
            )
        left = 0.065 + index * 0.315
        panels.append(
            Panel(
                (left, 0.14, 0.27, 0.72),
                (-1.08, 1.08),
                (-1.08, 1.08),
                "x / R",
                "y / R",
                panel_series,
                f"({chr(ord('a') + index)}) \\Delta={delta_deg:g} deg",
                legend=index == 0,
                legend_location="lower left" if index == 0 else "upper right",
                inset=inset,
                equal_aspect=True,
                xticks=[-1.0, -0.5, 0.0, 0.5, 1.0],
                yticks=[-1.0, -0.5, 0.0, 0.5, 1.0],
            )
        )
        n_rows = max(len(x_coordinate), len(x_proper))
        pad_coordinate = np.full(n_rows, np.nan)
        pad_proper = np.full(n_rows, np.nan)
        pad_coordinate_y = np.full(n_rows, np.nan)
        pad_proper_y = np.full(n_rows, np.nan)
        pad_coordinate[: len(x_coordinate)] = x_coordinate
        pad_coordinate_y[: len(y_coordinate)] = y_coordinate
        pad_proper[: len(x_proper)] = x_proper
        pad_proper_y[: len(y_proper)] = y_proper
        save_csv(
            os.path.join(datadir, f"fig_03_curves_delta_{safe_name(delta_deg)}.csv"),
            [
                ("x_coordinate_over_R", pad_coordinate),
                ("y_coordinate_over_R", pad_coordinate_y),
                ("x_proper_over_R", pad_proper),
                ("y_proper_over_R", pad_proper_y),
                ("q_coordinate", np.full(n_rows, q_coordinate)),
                ("q_proper", np.full(n_rows, q_proper)),
                ("delta_coordinate_deg", np.full(n_rows, math.degrees(delta_coordinate))),
                ("delta_proper_deg", np.full(n_rows, math.degrees(delta_proper))),
            ],
        )
        all_columns.append((f"q_coordinate_delta_{safe_name(delta_deg)}", np.array([q_coordinate])))
        all_columns.append((f"q_proper_delta_{safe_name(delta_deg)}", np.array([q_proper])))

    base = "fig_03_brachistochrone_curves"
    save_figure(outdir, base, panels, (1450, 520), (7.0, 2.55))

    # Figura singola colonna: solo l'ultimo pannello (Δ più grande), senza tag
    last = panels[-1]
    last_delta = deltas_deg[-1]
    single_panel = Panel(
        (0.13, 0.08, 0.80, 0.80),
        last.xlim,
        last.ylim,
        last.xlabel,
        last.ylabel,
        last.series,
        f"$\\Delta = {last_delta:g}^\\circ$",
        legend=True,
        legend_location="lower left",
        inset=last.inset,
        equal_aspect=True,
        xticks=last.xticks,
        yticks=last.yticks,
    )
    single_base = f"fig_03_brachistochrone_delta{safe_name(last_delta)}"
    save_figure(outdir, single_base, [single_panel], (1013, 1013), (3.375, 3.375))


def make_turning(sol: TOVSolution, outdir: str, datadir: str, deltas_deg: np.ndarray) -> Dict[str, np.ndarray]:
    scan = aperture_scan(sol, deltas_deg)
    series = [
        Series(scan["delta_deg"], scan["q_coordinate"], "coordinate", PALETTE["coordinate"], 1.65),
        Series(scan["delta_deg"], scan["q_proper"], "proper", PALETTE["proper"], 1.65, "dash"),
    ]
    panel = Panel(
        (0.18, 0.08, 0.74, 0.66),
        (float(np.min(deltas_deg)), float(np.max(deltas_deg))),
        (0.0, 1.02),
        "\\Delta [deg]",
        "$r_*/R$",
        series,
        "",
        xticks=[20, 50, 80, 110, 140, 170],
        yticks=[0.0, 0.25, 0.5, 0.75, 1.0],
    )
    base = "fig_04_turning_radius_vs_delta"
    save_figure(outdir, base, [panel], (1050, 560), (3.37, 2.05))
    save_csv(
        os.path.join(datadir, f"{base}.csv"),
        [
            ("delta_deg", scan["delta_deg"]),
            ("q_coordinate", scan["q_coordinate"]),
            ("q_proper", scan["q_proper"]),
            ("q_proper_minus_q_coordinate", scan["q_proper"] - scan["q_coordinate"]),
        ],
    )
    return scan


def make_times(scan: Dict[str, np.ndarray], outdir: str, datadir: str) -> None:
    left_series = [
        Series(scan["delta_deg"], scan["time_coordinate"], "t c / R", PALETTE["coordinate"], 1.65),
        Series(scan["delta_deg"], scan["time_proper"], "\\tau c / R", PALETTE["proper"], 1.65, "dash"),
    ]
    ratio = scan["time_coordinate"] / scan["time_proper"]
    right_series = [Series(scan["delta_deg"], ratio, "t / \\tau", PALETTE["black"], 1.65)]
    left_xlim, left_ylim = limits_from_series(left_series, pad_x=0.0, pad_y=0.24)
    right_xlim, right_ylim = limits_from_series(right_series, pad_x=0.0, pad_y=0.24)
    panels = [
        Panel(
            (0.09, 0.17, 0.39, 0.66),
            left_xlim,
            left_ylim,
            "\\Delta [deg]",
            "travel time",
            left_series,
            "(a)",
            xticks=[20, 60, 100, 140, 170],
        ),
        Panel(
            (0.60, 0.17, 0.34, 0.66),
            right_xlim,
            right_ylim,
            "\\Delta [deg]",
            "time ratio",
            right_series,
            "(b)",
            legend=False,
            xticks=[20, 60, 100, 140, 170],
        ),
    ]
    base = "fig_05_travel_times_vs_delta"
    save_figure("{}".format(outdir), base, panels, (1250, 520), (7.0, 2.9))
    save_csv(
        os.path.join(datadir, f"{base}.csv"),
        [
            ("delta_deg", scan["delta_deg"]),
            ("t_coordinate_c_over_R", scan["time_coordinate"]),
            ("tau_proper_c_over_R", scan["time_proper"]),
            ("t_over_tau", ratio),
        ],
    )


def make_compactness_scan(
    eos: PolytropicEOS | TabulatedEOS,
    central_pressures: Sequence[float],
    dr: float,
    delta_deg: float,
    outdir: str,
    datadir: str,
) -> None:
    solutions = compactness_sequence(central_pressures, eos=eos, dr=dr)
    compactness = np.array([sol.compactness for sol in solutions])
    radii = np.array([sol.radius for sol in solutions])
    masses = np.array([sol.total_mass for sol in solutions])
    central_densities = np.array([float(eos.rho_from_pressure(pressure)) for pressure in central_pressures])
    q_coordinate = np.full(len(solutions), np.nan)
    q_proper = np.full(len(solutions), np.nan)
    t_coordinate = np.full(len(solutions), np.nan)
    t_proper = np.full(len(solutions), np.nan)
    delta = math.radians(delta_deg)
    for index, sol in enumerate(solutions):
        scan = aperture_scan(sol, [delta_deg])
        q_coordinate[index] = scan["q_coordinate"][0]
        q_proper[index] = scan["q_proper"][0]
        t_coordinate[index] = scan["time_coordinate"][0]
        t_proper[index] = scan["time_proper"][0]

    left_series = [
        Series(compactness, q_coordinate, "coordinate", PALETTE["coordinate"], 1.65),
        Series(compactness, q_proper, "proper", PALETTE["proper"], 1.65, "dash"),
    ]
    right_series = [
        Series(compactness, t_coordinate, "t c / R", PALETTE["coordinate"], 1.65),
        Series(compactness, t_proper, "\\tau c / R", PALETTE["proper"], 1.65, "dash"),
    ]
    _, left_ylim = limits_from_series(left_series, pad_x=0.04, pad_y=0.24)
    _, right_ylim = limits_from_series(right_series, pad_x=0.04, pad_y=0.24)
    xlim = (float(np.min(compactness)) * 0.96, float(np.max(compactness)) * 1.04)
    panels = [
        Panel(
            (0.09, 0.17, 0.39, 0.66),
            xlim,
            left_ylim,
            "\\mu = M / R [1]",
            "r* / R",
            left_series,
            f"(a) \\Delta={delta_deg:g} deg",
        ),
        Panel(
            (0.60, 0.17, 0.34, 0.66),
            xlim,
            right_ylim,
            "\\mu = M / R [1]",
            "travel time",
            right_series,
            "(b)",
        ),
    ]
    base = f"fig_06_compactness_scan_delta_{safe_name(delta_deg)}"
    save_figure(outdir, base, panels, (1250, 520), (7.0, 2.9))
    save_csv(
        os.path.join(datadir, f"{base}.csv"),
        [
            ("central_pressure", np.asarray(central_pressures, dtype=float)),
            ("central_density", central_densities),
            ("mu_M_over_R", compactness),
            ("radius", radii),
            ("mass", masses),
            ("mass_Msun", masses / M_SUN_KM),
            ("q_coordinate", q_coordinate),
            ("q_proper", q_proper),
            ("t_coordinate_c_over_R", t_coordinate),
            ("tau_proper_c_over_R", t_proper),
            ("q_proper_minus_q_coordinate", q_proper - q_coordinate),
            ("t_over_tau", t_coordinate / t_proper),
        ],
    )


def find_pc_for_compactness(
    target_mu: float,
    eos: PolytropicEOS | TabulatedEOS,
    dr: float,
    pc_lo: float = 5.0e-5,
    pc_hi: float = 8.0e-3,
    n_iter: int = 52,
) -> float:
    def _mu(pc: float) -> float:
        try:
            return integrate_tov(pc, eos=eos, dr=dr).compactness
        except Exception:
            return float("nan")

    mu_lo, mu_hi = _mu(pc_lo), _mu(pc_hi)
    if not (math.isfinite(mu_lo) and math.isfinite(mu_hi)):
        raise ValueError("Boundary TOV integrations failed for alt EOS.")
    if not (min(mu_lo, mu_hi) < target_mu < max(mu_lo, mu_hi)):
        raise ValueError(
            "Cannot bracket target_mu=%.4f; alt EOS gives mu in [%.4f, %.4f]."
            % (target_mu, min(mu_lo, mu_hi), max(mu_lo, mu_hi))
        )
    for _ in range(n_iter):
        pc_mid = math.sqrt(pc_lo * pc_hi)
        mu_mid = _mu(pc_mid)
        if not math.isfinite(mu_mid):
            pc_hi = pc_mid
            continue
        if (mu_mid - target_mu) * (mu_lo - target_mu) <= 0.0:
            pc_hi, mu_hi = pc_mid, mu_mid
        else:
            pc_lo, mu_lo = pc_mid, mu_mid
    return math.sqrt(pc_lo * pc_hi)


def make_eos_comparison(
    sol_ref: TOVSolution,
    eos_alt: PolytropicEOS | TabulatedEOS,
    dr: float,
    outdir: str,
    datadir: str,
    delta_deg: float = 135.0,
) -> None:
    """fig_07: single-panel square EOS comparison with turning-point inset."""
    try:
        pc_alt = find_pc_for_compactness(sol_ref.compactness, eos_alt, dr)
    except ValueError as exc:
        print("EOS comparison skipped: %s" % exc)
        return
    sol_alt = integrate_tov(pc_alt, eos=eos_alt, dr=dr)
    alt_name = (getattr(eos_alt, "name", "alt") or "alt")[:7]

    delta = math.radians(delta_deg)
    x_rt, y_rt, q_rt, _ = full_curve_xy(sol_ref, "coordinate", delta)
    x_rp, y_rp, q_rp, _ = full_curve_xy(sol_ref, "proper", delta)
    x_at, y_at, q_at, _ = full_curve_xy(sol_alt, "coordinate", delta)
    x_ap, y_ap, q_ap, _ = full_curve_xy(sol_alt, "proper", delta)

    theta = np.linspace(0.0, 2.0 * math.pi, 900)
    surf_x, surf_y = np.cos(theta), np.sin(theta)

    series: List[Series] = [
        Series(surf_x, surf_y, "surface", PALETTE["surface"], 1.0, "dot"),
        Series(x_rt, y_rt, "poly t", PALETTE["coordinate"], 1.9),
        Series(x_rp, y_rp, "poly \\tau", PALETTE["proper"], 1.9, "dash"),
        Series(x_at, y_at, "%s t" % alt_name, "#009E73", 1.9, "dashdot"),
        Series(x_ap, y_ap, "%s \\tau" % alt_name, "#E69F00", 1.9, "dot"),
    ]

    q_min = min(float(q_rp), float(q_ap))
    q_max = max(float(q_rt), float(q_at))
    gap = max(q_max - q_min, 0.04)
    inset_xlim = (q_min - gap * 0.6, q_max + gap * 0.6)
    inset_ylim = (-0.28, 0.28)
    inset_series = [
        Series(x_rt, y_rt, "", PALETTE["coordinate"], 1.9),
        Series(x_rp, y_rp, "", PALETTE["proper"], 1.9, "dash"),
        Series(x_at, y_at, "", "#009E73", 1.9, "dashdot"),
        Series(x_ap, y_ap, "", "#E69F00", 1.9, "dot"),
    ]
    inset = Inset(
        rect=(0.04, 0.04, 0.42, 0.40),
        xlim=inset_xlim,
        ylim=inset_ylim,
        label="zoom: turning pt.",
        series=inset_series,
    )

    panel = Panel(
        (0.14, 0.07, 0.82, 0.82),
        (-1.08, 1.08), (-1.08, 1.08),
        "x / R", "y / R",
        series,
        "",
        legend=True,
        legend_location="lower left",
        inset=inset,
        equal_aspect=True,
        xticks=[-1.0, -0.5, 0.0, 0.5, 1.0],
        yticks=[-1.0, -0.5, 0.0, 0.5, 1.0],
    )
    base = "fig_07_eos_comparison"
    save_figure(outdir, base, [panel], (700, 700), (3.375, 3.375))
    save_csv(
        os.path.join(datadir, "%s.csv" % base),
        [
            ("eos_ref", np.array([sol_ref.eos.name])),
            ("eos_alt", np.array([sol_alt.eos.name])),
            ("mu_ref", np.array([sol_ref.compactness])),
            ("mu_alt", np.array([sol_alt.compactness])),
            ("q_poly_t", np.array([float(q_rt)])),
            ("q_poly_tau", np.array([float(q_rp)])),
            ("q_alt_t", np.array([float(q_at)])),
            ("q_alt_tau", np.array([float(q_ap)])),
        ],
    )
    print(
        "EOS comparison: %s mu=%.4f  vs  %s mu=%.4f"
        % (sol_ref.eos.name, sol_ref.compactness, sol_alt.eos.name, sol_alt.compactness)
    )


def make_depth_difference(
    eos: PolytropicEOS | TabulatedEOS,
    central_pressures: Sequence[float],
    dr: float,
    deltas_deg: Sequence[float],
    outdir: str,
    datadir: str,
) -> None:
    """fig_08: q_t - q_tau vs compactness for multiple Delta values."""
    solutions = compactness_sequence(list(central_pressures), eos=eos, dr=dr)
    mu = np.array([2.0 * sol.compactness for sol in solutions])   # μ = r_s/R = 2M/R
    xlim = (float(np.min(mu)) * 0.96, float(np.max(mu)) * 1.04)

    _pal = ["#111111", "#0072B2", "#D55E00", "#009E73", "#CC79A7"]
    series_list: List[Series] = []
    csv_cols: List[Tuple[str, np.ndarray]] = [("mu_rs_over_R", mu)]

    for i, delta_deg in enumerate(deltas_deg):
        diff = np.full(len(solutions), np.nan)
        for j, sol in enumerate(solutions):
            scan = aperture_scan(sol, [float(delta_deg)])
            diff[j] = float(scan["q_coordinate"][0]) - float(scan["q_proper"][0])
        series_list.append(Series(
            mu, diff, "\\Delta=%d" % int(delta_deg),
            _pal[i % len(_pal)], 1.65,
        ))
        csv_cols.append(("diff_q_delta_%s" % safe_name(float(delta_deg)), diff))

    _, ylim = limits_from_series(series_list, pad_x=0.0, pad_y=0.12)
    ylim = (max(0.0, ylim[0]), ylim[1])

    panel = Panel(
        (0.18, 0.14, 0.74, 0.70),
        xlim, ylim,
        "\\mu = r_s / R",
        "q_{t} - q_{\\tau}",
        series_list,
        "",
    )
    base = "fig_08_depth_difference"
    save_figure(outdir, base, [panel], (1050, 560), (3.37, 2.4))
    save_csv(os.path.join(datadir, "%s.csv" % base), csv_cols)


def make_compactness_scan_multidelta(
    eos: PolytropicEOS | TabulatedEOS,
    central_pressures: Sequence[float],
    dr: float,
    deltas_deg: Sequence[float],
    outdir: str,
    datadir: str,
) -> None:
    """fig_09: r*/R and t/tau vs compactness for multiple Delta values."""
    solutions = compactness_sequence(list(central_pressures), eos=eos, dr=dr)
    mu = np.array([2.0 * sol.compactness for sol in solutions])   # μ = r_s/R = 2M/R
    xlim = (float(np.min(mu)) * 0.96, float(np.max(mu)) * 1.04)

    _pal = [PALETTE["coordinate"], PALETTE["accent"], PALETTE["proper"]]
    left_series: List[Series] = []
    right_series: List[Series] = []
    csv_cols: List[Tuple[str, np.ndarray]] = [("mu_rs_over_R", mu)]

    for i, delta_deg in enumerate(deltas_deg):
        col = _pal[i % len(_pal)]
        q_t = np.full(len(solutions), np.nan)
        q_p = np.full(len(solutions), np.nan)
        ratio = np.full(len(solutions), np.nan)
        for j, sol in enumerate(solutions):
            scan = aperture_scan(sol, [float(delta_deg)])
            q_t[j] = scan["q_coordinate"][0]
            q_p[j] = scan["q_proper"][0]
            tc = float(scan["time_coordinate"][0])
            tp = float(scan["time_proper"][0])
            if math.isfinite(tc) and math.isfinite(tp) and tp > 0.0:
                ratio[j] = tc / tp
        lbl = "\\Delta=%d" % int(delta_deg)
        left_series.append(Series(mu, q_t, "t %s" % lbl, col, 1.55))
        left_series.append(Series(mu, q_p, "\\tau %s" % lbl, col, 1.55, "dash"))
        right_series.append(Series(mu, ratio, lbl, col, 1.65))
        csv_cols += [
            ("q_t_delta_%s" % safe_name(float(delta_deg)), q_t),
            ("q_tau_delta_%s" % safe_name(float(delta_deg)), q_p),
            ("t_over_tau_delta_%s" % safe_name(float(delta_deg)), ratio),
        ]

    _, left_ylim = limits_from_series(left_series, pad_x=0.0, pad_y=0.10)
    _, right_ylim = limits_from_series(right_series, pad_x=0.0, pad_y=0.12)

    panels = [
        Panel(
            (0.09, 0.17, 0.39, 0.66),
            xlim, left_ylim,
            "\\mu = r_s / R", "r* / R",
            left_series, "(a)",
        ),
        Panel(
            (0.60, 0.17, 0.34, 0.66),
            xlim, right_ylim,
            "\\mu = r_s / R", "t / \\tau",
            right_series, "(b)",
        ),
    ]
    base = "fig_09_compactness_multidelta"
    save_figure(outdir, base, panels, (1250, 520), (7.0, 2.9))
    save_csv(os.path.join(datadir, "%s.csv" % base), csv_cols)


_SEQ_PAL = ["#313695", "#4575b4", "#56B4E9", "#E69F00", "#D55E00", "#CC79A7"]


def make_gamma_scan(
    central_pressure: float,
    dr: float,
    delta_deg: float,
    gamma_values: Sequence[float],
    mode: str,
    outdir: str,
    datadir: str,
    kappa: float = 100.0,
) -> None:
    """Square figure: brachistochrones at fixed central pressure for multiple Gamma values."""
    delta = math.radians(delta_deg)
    mode_tag = "t" if mode == "coordinate" else "tau"
    mode_sym = "t" if mode == "coordinate" else "\\tau"
    theta = np.linspace(0.0, 2.0 * math.pi, 900)
    surf_x, surf_y = np.cos(theta), np.sin(theta)

    n_gamma = len(gamma_values)
    curve_data: List[tuple] = []
    for i, gamma in enumerate(gamma_values):
        eos_g = PolytropicEOS(kappa=kappa, gamma=gamma)
        try:
            sol_g = integrate_tov(central_pressure, eos=eos_g, dr=dr)
            x, y, q, _ = full_curve_xy(sol_g, mode, delta)
            label = "\\Gamma=%.1f  (\\mu=%.3f)" % (gamma, sol_g.compactness)
            col_idx = round(i * (len(_SEQ_PAL) - 1) / max(n_gamma - 1, 1))
            curve_data.append((x, y, float(q), label, _SEQ_PAL[col_idx]))
        except Exception as exc:
            print("gamma_scan: Gamma=%.1f skipped: %s" % (gamma, exc))

    series_list: List[Series] = [
        Series(surf_x, surf_y, "surface", PALETTE["surface"], 1.0, "dot"),
    ]
    inset_series: List[Series] = []
    for x, y, q, label, col in curve_data:
        series_list.append(Series(x, y, label, col, 1.75))
        inset_series.append(Series(x, y, "", col, 1.75))

    inset = None
    if curve_data:
        q_vals = [d[2] for d in curve_data]
        q_min, q_max = min(q_vals), max(q_vals)
        gap = max(q_max - q_min, 0.04)
        inset = Inset(
            rect=(0.14, 0.04, 0.44, 0.42),
            xlim=(q_min - gap * 0.8, q_max + gap * 0.8),
            ylim=(-0.28, 0.28),
            label="zoom: turning pt.",
            series=inset_series,
            show_yticks=True,
        )

    panel = Panel(
        (0.14, 0.07, 0.82, 0.82),
        (-1.08, 1.08), (-1.08, 1.08),
        "x / R", "y / R",
        series_list,
        "",
        legend=True,
        legend_location="lower left",
        inset=inset,
        equal_aspect=True,
        xticks=[-1.0, -0.5, 0.0, 0.5, 1.0],
        yticks=[-1.0, -0.5, 0.0, 0.5, 1.0],
    )
    base = "fig_10_gamma_scan_%s" % mode_tag
    if mode == "proper":
        base = "fig_11_gamma_scan_%s" % mode_tag
    save_figure(outdir, base, [panel], (700, 700), (3.375, 3.375))


def make_compose_curves(
    eos_compose: TabulatedEOS,
    dr: float,
    target_mu: float,
    delta_values: Sequence[float],
    mode: str,
    outdir: str,
    datadir: str,
) -> None:
    """Square figure: brachistochrones for CompOSE EOS at multiple Delta values, fixed mu."""
    try:
        pc = find_pc_for_compactness(target_mu, eos_compose, dr)
        sol_c = integrate_tov(pc, eos=eos_compose, dr=dr)
    except Exception as exc:
        print("compose_curves: could not build solution (%s)." % exc)
        return

    mode_tag = "t" if mode == "coordinate" else "tau"
    mode_sym = "t" if mode == "coordinate" else "\\tau"
    theta = np.linspace(0.0, 2.0 * math.pi, 900)
    surf_x, surf_y = np.cos(theta), np.sin(theta)

    series_list: List[Series] = [
        Series(surf_x, surf_y, "surface", PALETTE["surface"], 1.0, "dot"),
    ]
    for i, delta_deg in enumerate(delta_values):
        x, y, _, _ = full_curve_xy(sol_c, mode, math.radians(float(delta_deg)))
        series_list.append(
            Series(x, y, "\\Delta=%d" % int(delta_deg), _SEQ_PAL[i % len(_SEQ_PAL)], 1.75)
        )

    base_idx = "12" if mode == "coordinate" else "13"
    panel = Panel(
        (0.14, 0.07, 0.82, 0.82),
        (-1.08, 1.08), (-1.08, 1.08),
        "x / R", "y / R",
        series_list,
        "%s-brachistochrone, \\mu=%.3f" % (mode_sym, sol_c.compactness),
        legend=True,
        legend_location="lower left",
        equal_aspect=True,
        xticks=[-1.0, -0.5, 0.0, 0.5, 1.0],
        yticks=[-1.0, -0.5, 0.0, 0.5, 1.0],
    )
    base = "fig_%s_compose_curves_%s" % (base_idx, mode_tag)
    save_figure(outdir, base, [panel], (700, 700), (3.375, 3.375))


# =============================================================================
# Schwarzschild validation: incompressible EOS + polytrope convergence
# =============================================================================

@dataclass
class IncompressibleEOS:
    """Constant energy density: exact Schwarzschild interior metric solution."""
    epsilon0: float = 1.0e-4
    name: str = "Schwarzschild (\\rho=const)"

    def epsilon_from_pressure(self, pressure: float | np.ndarray) -> float | np.ndarray:
        if np.isscalar(pressure):
            return self.epsilon0
        return np.full_like(np.asarray(pressure, dtype=float), self.epsilon0)

    def rho_from_pressure(self, pressure: float | np.ndarray) -> float | np.ndarray:
        return self.epsilon_from_pressure(pressure)


def _incompressible_scan(
    epsilon0: float,
    dr: float,
    delta: float,
    mu_min: float = 0.02,
    mu_max: float = 0.38,
    n_pc: int = 60,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Scan pc for IncompressibleEOS; return (mu, q_t, q_tau) sorted by mu."""
    eos = IncompressibleEOS(epsilon0=epsilon0)
    # Analytical p_c for constant-density star: mu here is r_s/R = 2M/R
    # p_c/eps0 = (1 - sqrt(1 - mu)) / (3*sqrt(1 - mu) - 1)
    def pc_from_mu(mu: float) -> float:
        s = math.sqrt(max(1.0 - mu, 0.0))
        return epsilon0 * (1.0 - s) / max(3.0 * s - 1.0, 1e-15)

    pc_lo = pc_from_mu(mu_min) * 0.5
    pc_hi = pc_from_mu(mu_max) * 2.0
    pc_values = np.geomspace(pc_lo, pc_hi, n_pc)

    mus, q_ts, q_taus = [], [], []
    for pc in pc_values:
        try:
            sol = integrate_tov(pc, eos=eos, dr=dr)
            mu_paper = 2.0 * sol.compactness   # μ = r_s/R
            if not (mu_min * 0.5 < mu_paper < mu_max * 1.1):
                continue
            _, _, qt, _ = full_curve_xy(sol, "coordinate", delta)
            _, _, qtau, _ = full_curve_xy(sol, "proper", delta)
            mus.append(mu_paper)
            q_ts.append(float(qt))
            q_taus.append(float(qtau))
        except Exception:
            pass

    order = np.argsort(mus)
    return np.array(mus)[order], np.array(q_ts)[order], np.array(q_taus)[order]


def _polytrope_scan(
    gamma: float,
    dr: float,
    delta: float,
    n_pc: int = 50,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Scan pc for PolytropicEOS on the STABLE branch only (μ monotone in pc).

    Past the maximum-mass star μ starts decreasing — we stop there to avoid
    integrating on the unstable branch where the brachistochrone is ill-defined.
    """
    eos = PolytropicEOS(kappa=100.0, gamma=gamma)
    pc_values = np.geomspace(1e-6, 5.0, n_pc)
    mus, q_ts, q_taus = [], [], []
    mu_max_seen = 0.0
    for pc in pc_values:
        try:
            sol = integrate_tov(pc, eos=eos, dr=dr)
            mu_paper = 2.0 * sol.compactness   # μ = r_s/R
            if not (0.01 < mu_paper < 0.88):
                continue
            # Stop if we have moved backward by >2% — past the mass maximum
            if mus and mu_paper < mu_max_seen * 0.98:
                break
            mu_max_seen = max(mu_max_seen, mu_paper)
            _, _, qt, _ = full_curve_xy(sol, "coordinate", delta)
            _, _, qtau, _ = full_curve_xy(sol, "proper", delta)
            qt, qtau = float(qt), float(qtau)
            # Sanity: both turning points must be valid fractions
            if not (0.0 < qtau < qt < 1.0):
                continue
            # Sanity: smooth variation — reject wild jumps
            if q_ts and abs(qt - q_ts[-1]) > 0.15:
                break
            mus.append(mu_paper)
            q_ts.append(qt)
            q_taus.append(qtau)
        except Exception:
            pass
    if not mus:
        return np.array([]), np.array([]), np.array([])
    order = np.argsort(mus)
    return np.array(mus)[order], np.array(q_ts)[order], np.array(q_taus)[order]


def make_schwarzschild_validation(
    gamma_values: Sequence[float],
    dr: float,
    delta_deg: float,
    outdir: str,
    datadir: str,
    mu_targets: Optional[Sequence[float]] = None,
) -> None:
    """fig_15: q_t and q_tau vs mu — polytropes vs Schwarzschild (rho=const)."""
    delta = math.radians(delta_deg)
    q0_newt = 1.0 - delta_deg / 180.0  # Newtonian limit for constant density

    PAL = ["#0072B2", "#D55E00", "#009E73", "#CC79A7"]

    series_t: List[Series] = []
    series_tau: List[Series] = []

    # Polytropic EOS curves (full pc scan to avoid missing high-mu region)
    mu_max_all = 0.0
    for idx, gamma in enumerate(gamma_values):
        col = PAL[idx % len(PAL)]
        label = "\\Gamma=%g" % gamma
        mus, qts, qtaus = _polytrope_scan(gamma, dr, delta, n_pc=50)
        if len(mus):
            mu_max_all = max(mu_max_all, float(np.max(mus)))
            series_t.append(Series(mus, qts, label, col, 1.6))
            series_tau.append(Series(mus, qtaus, label, col, 1.6))

    # Schwarzschild reference (incompressible)
    mus_s, qts_s, qtaus_s = _incompressible_scan(
        epsilon0=1.0e-4, dr=dr, delta=delta,
        mu_min=0.02, mu_max=0.84,   # in μ=r_s/R units; Buchdahl = 8/9 ≈ 0.889
    )
    if len(mus_s):
        mu_max_all = max(mu_max_all, float(np.max(mus_s)))
        schw_label = "Schwarzschild (\\rho=const)"
        series_t.append(Series(mus_s, qts_s, schw_label, "#111111", 2.4))
        series_tau.append(Series(mus_s, qtaus_s, schw_label, "#111111", 2.4))

    # Newtonian limit for constant density (horizontal dashed)
    mu_line = np.array([0.0, mu_max_all * 1.08])
    for sl in (series_t, series_tau):
        sl.append(Series(mu_line, np.full(2, q0_newt),
                         "Newtonian (\\rho=const)", "#888888", 0.9, "dash"))

    all_qt  = np.concatenate([s.y[np.isfinite(s.y)] for s in series_t])
    all_qtau = np.concatenate([s.y[np.isfinite(s.y)] for s in series_tau])
    ylo_t  = float(np.nanmin(all_qt))  - 0.02
    yhi_t  = float(np.nanmax(all_qt))  + 0.02
    ylo_tau = float(np.nanmin(all_qtau)) - 0.02
    yhi_tau = float(np.nanmax(all_qtau)) + 0.02
    xlim = (0.0, mu_max_all * 1.06)

    top_rect = (0.16, 0.04, 0.78, 0.44)
    bot_rect = (0.16, 0.52, 0.78, 0.44)

    panel_t = Panel(
        top_rect, xlim, (ylo_t, yhi_t),
        "", "$q^{(t)}$",
        series_t,
        "\\Delta = %g deg" % delta_deg,
        tag_loc="upper right",
        legend=True, legend_location="upper left",
    )
    panel_tau = Panel(
        bot_rect, xlim, (ylo_tau, yhi_tau),
        "$\\mu = r_s/R$", "$q^{(\\tau)}$",
        series_tau,
        "",
        legend=False,
    )

    base = "fig_15_schwarzschild_validation"
    save_figure(outdir, base, [panel_t, panel_tau], (700, 950), (3.375, 4.6))
    print("schwarzschild validation: %s/%s.png" % (outdir, base))


# =============================================================================
# Weak-field t–τ splitting formula  (Sec. weak_field)
# Convention: μ_paper = r_s/R = 2M/R = 2 * sol.compactness
# W(x) = ∫_x^1 m̄(x')/x'² dx',  P₀²(x) = x²/W(x)
# =============================================================================

def _build_W_interp(sol: TOVSolution):
    """Interpolant for W(x) = ∫_x^1 m̄(x')/x'² dx' (right-to-left trapz)."""
    from scipy.interpolate import interp1d
    x    = sol.r / sol.radius
    mbar = sol.mass / sol.total_mass
    with np.errstate(divide="ignore", invalid="ignore"):
        intgd = np.where(x > 1e-10, mbar / x**2, 0.0)
    W = np.zeros_like(x)
    for i in range(len(x) - 2, -1, -1):
        W[i] = W[i + 1] + 0.5 * (intgd[i] + intgd[i + 1]) * (x[i + 1] - x[i])
    W = np.clip(W, 0.0, None)
    return interp1d(x, W, kind="cubic", bounds_error=False,
                    fill_value=(float(W[0]), 0.0))


def _P0sq(x: float, W_func) -> float:
    W_x = float(W_func(x))
    return x**2 / W_x if W_x > 1e-30 else 1e30


def _F0_integral(q0: float, W_func) -> float:
    """F⁽⁰⁾[q₀] = ∫_{q₀}^1 √P₀²(q₀) / (x √(P₀²(x)−P₀²(q₀))) dx."""
    import warnings
    from scipy import integrate as sci_int
    P0sq_q = _P0sq(q0, W_func)
    def intgd(x: float) -> float:
        d = _P0sq(x, W_func) - P0sq_q
        return math.sqrt(P0sq_q) / (x * math.sqrt(d)) if d > 1e-30 else 0.0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        val, _ = sci_int.quad(intgd, q0 * (1.0 + 1e-7), 1.0 - 1e-6,
                              limit=200, epsabs=1e-6, epsrel=1e-6)
    return val


def _find_q0(delta: float, W_func) -> float:
    """Solve F⁽⁰⁾[q₀] = Δ/2 for the Newtonian turning point."""
    from scipy.optimize import brentq
    return brentq(lambda q: _F0_integral(q, W_func) - delta / 2.0,
                  1e-4, 1.0 - 1e-4, xtol=1e-7)


def _G_functional(q0: float, W_func) -> float:
    """G[q₀] = ∫_{q₀}^1 I(x;q₀)·P₀²(x)·(W(q₀)−W(x))/(P₀²(x)−P₀²(q₀)) dx."""
    import warnings
    from scipy import integrate as sci_int
    P0sq_q = _P0sq(q0, W_func)
    W_q0   = float(W_func(q0))
    def intgd(x: float) -> float:
        W_x    = float(W_func(x))
        P0sq_x = _P0sq(x, W_func)
        d = P0sq_x - P0sq_q
        if d <= 1e-30:
            return 0.0
        return (math.sqrt(P0sq_q) / (x * math.sqrt(d))) * P0sq_x * (W_q0 - W_x) / d
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        val, _ = sci_int.quad(intgd, q0 * (1.0 + 1e-7), 1.0 - 1e-6,
                              limit=200, epsabs=1e-6, epsrel=1e-6)
    return val


def _dF0dq(q0: float, W_func, h: float = 3e-4) -> float:
    """∂_q F⁽⁰⁾(q₀) by central finite difference."""
    return (_F0_integral(q0 + h, W_func) - _F0_integral(q0 - h, W_func)) / (2.0 * h)


def _build_mbar_interp(sol: TOVSolution):
    """Interpolant for m̄(x) = m(r)/M (normalised mass profile)."""
    from scipy.interpolate import interp1d
    x    = sol.r / sol.radius
    mbar = sol.mass / sol.total_mass
    return interp1d(x, mbar, kind="cubic", bounds_error=False,
                    fill_value=(0.0, 1.0))


def _M_integral(q0: float, W_func, mbar_func) -> float:
    """M[q₀] = ∫_{q₀}^1 I₀(x;q₀)·m̄(x)/(2x) dx.

    Arises from the O(μ) correction due to the radial-metric factor √b ≈ 1+μm̄/(2x).
    This term is identical for the t and τ Fermat integrals, so it cancels in G
    but enters individually in Φ₁^t = 3G/4 + M and Φ₁^τ = −G/4 + M.
    """
    import warnings
    from scipy import integrate as sci_int
    P0sq_q = _P0sq(q0, W_func)
    def intgd(x: float) -> float:
        d = _P0sq(x, W_func) - P0sq_q
        if d <= 1e-30:
            return 0.0
        return math.sqrt(P0sq_q) / (x * math.sqrt(d)) * float(mbar_func(x)) / (2.0 * x)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        val, _ = sci_int.quad(intgd, q0 * (1.0 + 1e-7), 1.0 - 1e-6,
                              limit=200, epsabs=1e-6, epsrel=1e-6)
    return val


def _schw_dq_o2_coeff(q0: float) -> Tuple[float, float]:
    """C₁, C₂ for Schwarzschild (uniform density): δq = C₁μ + C₂μ².

    Uses analytical W(x)=(1−x²)/2 and m̄(x)=x³ (no interpolation).
    """
    import warnings
    # Analytical uniform-density profile: W(x) = (1-x^2)/2, mbar(x) = x^3
    def W_func(x):
        x = np.asarray(x, float)
        return np.where(x <= 1.0, (1.0 - x**2) / 2.0, 0.0)
    def mbar_func(x):
        x = np.asarray(x, float)
        return np.minimum(x**3, 1.0)

    # For the Schwarzschild profile P₀²(x)→∞ near x=1 causes spurious
    # convergence warnings in the Romberg table; the results are accurate.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        G   = _G_functional(q0, W_func)
        M   = _M_integral(q0, W_func, mbar_func)
        dF  = _dF0dq(q0, W_func)
        if abs(dF) < 1e-15:
            return (0.0, 0.0)
        h   = 3e-4
        G_p = (_G_functional(q0 + h, W_func) - _G_functional(q0 - h, W_func)) / (2 * h)
        M_p = (_M_integral(q0 + h, W_func, mbar_func) - _M_integral(q0 - h, W_func, mbar_func)) / (2 * h)
        dF2 = (_dF0dq(q0 + h, W_func) - _dF0dq(q0 - h, W_func)) / (2 * h)
    C1 = -G / dF
    C2 = (G * G_p / 2.0 + G * M_p + M * G_p) / dF**2 \
         - G * (G / 2.0 + 2.0 * M) * dF2 / (2.0 * dF**3)
    return C1, C2


def dq_weakfield(sol: TOVSolution, delta: float) -> float:
    """O(μ) formula: q^(t)−q^(τ) = −μ·G[q₀]/∂_q F⁽⁰⁾(q₀)."""
    W_func = _build_W_interp(sol)
    q0     = _find_q0(delta, W_func)
    G      = _G_functional(q0, W_func)
    dF     = _dF0dq(q0, W_func)
    if abs(dF) < 1e-15:
        return float("nan")
    return -2.0 * sol.compactness * G / dF


def dq_weakfield_o2(sol: TOVSolution, delta: float) -> float:
    """O(μ²) formula for q^(t)−q^(τ).

    Extends the O(μ) result by including the kinematic second-order terms that
    arise when q₁^t and q₁^τ feed back through the perturbed Fermat conditions:

        δq = −μ G/F′₀  +  μ² [(G G′/2 + G M′ + M G′)/F′₀²
                                − G(G/2 + 2M) F″₀ / (2 F′₀³)]

    where  G  = _G_functional (first-order optical-index splitting),
           M  = _M_integral   (radial-metric √b correction, same for t and τ),
           and primes denote d/dq₀.

    Derivation: expand n_t and n_τ to O(μ) in the lapse A ≈ A₀ exp(−μW),
    insert into the Fermat angle condition, collect at each order.
    Individual first-order functionals: Φ₁^t = 3G/4 + M, Φ₁^τ = −G/4 + M.
    The intrinsic O(μ²) correction to the lapse (from second-order TOV) is
    dropped (it requires EOS-specific pressure integrals and is numerically
    small for μ ≲ 0.3).
    """
    W_func    = _build_W_interp(sol)
    mbar_func = _build_mbar_interp(sol)
    q0        = _find_q0(delta, W_func)
    mu        = 2.0 * sol.compactness

    G  = _G_functional(q0, W_func)
    M  = _M_integral(q0, W_func, mbar_func)
    dF = _dF0dq(q0, W_func)
    if abs(dF) < 1e-15:
        return float("nan")

    h    = 3e-4
    G_p  = (_G_functional(q0 + h, W_func) - _G_functional(q0 - h, W_func)) / (2 * h)
    M_p  = (_M_integral(q0 + h, W_func, mbar_func) - _M_integral(q0 - h, W_func, mbar_func)) / (2 * h)
    dF2  = (_dF0dq(q0 + h, W_func) - _dF0dq(q0 - h, W_func)) / (2 * h)

    dq1 = -mu * G / dF
    dq2 = mu**2 * (
        (G * G_p / 2.0 + G * M_p + M * G_p) / dF**2
        - G * (G / 2.0 + 2.0 * M) * dF2 / (2.0 * dF**3)
    )
    return dq1 + dq2


def make_t_tau_splitting(
    eos_pairs: Sequence[Tuple[str, Any]],
    dr: float,
    delta_deg: float,
    outdir: str,
    datadir: str,
    mu_targets: Optional[Sequence[float]] = None,
    compute_weakfield: bool = True,
) -> None:
    """fig_14a/b: q^(t)−q^(τ) vs μ — two versions with O(μ) and O(μ²) weak-field.

    Both versions share the same numerical data; only the analytic comparison
    curves differ.  Outputs:
        fig_14_t_tau_splitting_o1.png  (O(μ) weak-field)
        fig_14_t_tau_splitting_o2.png  (O(μ²) weak-field)
    """
    if mu_targets is None:
        mu_targets = list(np.linspace(0.09, 0.22, 14))
    delta   = math.radians(delta_deg)
    q0_schw = 1.0 - delta_deg / 180.0          # Schwarzschild: q₀ = 1 − Δ/π

    PAL = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#56B4E9"]
    mu_max = 0.0

    # Collect per-EOS data in one pass (expensive TOV/brachistochrone integrals)
    eos_data: List[Tuple] = []   # (name, col, mu_arr, dq_num, dq_o1, dq_o2)
    for idx, (eos_name, eos) in enumerate(eos_pairs):
        col = PAL[idx % len(PAL)]
        mu_num: List[float] = []
        dq_num: List[float] = []
        dq_o1:  List[float] = []
        dq_o2:  List[float] = []
        for mu_t in mu_targets:
            try:
                pc  = find_pc_for_compactness(mu_t, eos, dr)
                sol = integrate_tov(pc, eos=eos, dr=dr)
                _, _, q_t,   _ = full_curve_xy(sol, "coordinate", delta)
                _, _, q_tau, _ = full_curve_xy(sol, "proper",     delta)
                mu_num.append(2.0 * sol.compactness)
                dq_num.append(float(q_t) - float(q_tau))
                if compute_weakfield:
                    dq_o1.append(dq_weakfield(sol, delta))
                    dq_o2.append(dq_weakfield_o2(sol, delta))
            except Exception as exc:
                print("splitting %s mu=%.3f: %s" % (eos_name, mu_t, exc))
        if not mu_num:
            continue
        mu_arr = np.array(mu_num)
        mu_max = max(mu_max, float(np.max(mu_arr)))
        eos_data.append((eos_name, col, mu_arr, np.array(dq_num),
                         np.array(dq_o1), np.array(dq_o2)))

    # Schwarzschild analytical curves
    mu_line  = np.linspace(0.0, (mu_max or 0.48) * 1.08, 300)
    # O(μ): classical formula
    dq_schw_o1 = mu_line * q0_schw * (1.0 - q0_schw**2)**2 / 4.0
    # O(μ²): extended formula using uniform-density coefficients
    C1_schw, C2_schw = _schw_dq_o2_coeff(q0_schw)
    dq_schw_o2 = C1_schw * mu_line + C2_schw * mu_line**2

    def _build_panel_and_save(order: int) -> None:
        label_wf = "O(\\mu)"   if order == 1 else "O(\\mu^2)"
        schw_y   = dq_schw_o1  if order == 1 else dq_schw_o2
        base     = "fig_14_t_tau_splitting_o%d" % order

        series_list: List[Series] = []
        for (eos_name, col, mu_arr, dq_num_arr, dq_o1_arr, dq_o2_arr) in eos_data:
            series_list.append(Series(mu_arr, dq_num_arr, eos_name, col, 1.8))
            wf_arr = dq_o1_arr if order == 1 else dq_o2_arr
            if compute_weakfield and len(wf_arr):
                series_list.append(
                    Series(mu_arr, wf_arr,
                           "%s  (%s)" % (eos_name, label_wf), col, 1.1, "dash"))
        series_list.append(
            Series(mu_line, schw_y,
                   "Schw.  (analytic, %s)" % label_wf, "#333333", 1.5, "dashdot"))

        all_y = np.concatenate([s.y[np.isfinite(s.y)] for s in series_list if len(s.y)])
        ylim  = (0.0, float(np.max(all_y)) * 1.12) if len(all_y) else (0.0, 0.1)
        xlim  = (0.0, (mu_max or 0.24) * 1.08)

        panel = Panel(
            (0.14, 0.08, 0.82, 0.84),
            xlim, ylim,
            "$\\mu = r_s/R$",
            "$q^{(t)} - q^{(\\tau)}$",
            series_list,
            "\\Delta = %g deg" % delta_deg,
            tag_loc="upper right",
            legend=True,
            legend_location="upper left",
        )
        save_figure(outdir, base, [panel], (1013, 920), (3.375, 3.1))
        print("t-tau splitting (%s): %s/%s.png" % (label_wf, outdir, base))

    _build_panel_and_save(1)   # O(μ)  version
    _build_panel_and_save(2)   # O(μ²) version


def write_readme(
    outdir: str,
    datadir: str,
    sol: TOVSolution,
    args: argparse.Namespace,
) -> None:
    path = os.path.join(outdir, "README.md")
    lines = [
        "# TOV brachistochrone publication figures",
        "",
        "Generated by `TOVInternalMetric/genera_grafici_tov.py`.",
        "",
        "Conventions:",
        "- Geometrized units `G = c = 1`.",
        "- Travel times are plotted as `t c / R` and `\\tau c / R`.",
        "- Vector `PDF` and `SVG` are generated for APS-style workflows.",
        "",
        "EOS:",
        f"- Type: `{args.eos}`",
        f"- Name: `{sol.eos.name}`",
    ]
    if args.eos == "polytrope":
        lines.extend(
            [
                f"- `K = {args.kappa:g}`",
                f"- `Gamma = {args.gamma:g}`",
            ]
        )
    else:
        lines.extend(
            [
                f"- Table: `{args.eos_table or default_table_path()}`" if args.eos == "table" else f"- CompOSE dir: `{args.compose_dir or default_compose_dir()}`",
                "- Expected table relation: `pressure = pressure(density)`." if args.eos == "table" else "- CompOSE input: `eos.nb.ns` + `eos.thermo.ns`.",
                "- If no epsilon column is present, the density column is interpreted as total energy density." if args.eos == "table" else "- CompOSE units are converted from MeV/fm^3 to km^-2 for TOV integration.",
            ]
        )
    lines.extend(
        [
            "",
        "Main model:",
        f"- `pc = {args.central_pressure:g}`",
        f"- central density input = `{args.central_density}`",
        f"- `R = {sol.radius:.8g}`",
        f"- `M = {sol.total_mass:.8g} km = {sol.total_mass / M_SUN_KM:.8g} M_sun`",
        f"- `\\mu = M/R = {sol.compactness:.8g}`",
        f"- `A(R) = {sol.surface_lapse:.8g}`",
        "",
        "Figures:",
        "- `fig_01_tov_profiles`: pressure, energy density, mass, metric functions and release kinematics.",
        "- `fig_02_effective_indices`: coordinate-time and proper-time effective indices.",
        "- `fig_03_brachistochrone_curves`: spatial curves at selected surface apertures.",
        "- `fig_04_turning_radius_vs_delta`: turning radius as a function of aperture.",
        "- `fig_05_travel_times_vs_delta`: coordinate and proper travel times.",
        "- `fig_06_compactness_scan_*`: compactness dependence along an EOS sequence.",
        "",
        f"CSV source data are in `{os.path.relpath(datadir, outdir)}`.",
        ]
    )
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")


def build_parser() -> argparse.ArgumentParser:
    script_dir = os.path.dirname(os.path.abspath(__file__))
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--figure",
        choices=[
            "all", "profiles", "indices", "curves", "turning", "times",
            "compactness", "comparison", "depth_diff", "multidelta",
            "gamma_t", "gamma_tau", "compose_curves", "splitting", "schwarzschild_val",
        ],
        default="all",
    )
    parser.add_argument("--outdir", default=os.path.join(script_dir, "figures"))
    parser.add_argument("--datadir", default=None)
    parser.add_argument("--eos", choices=["polytrope", "table", "compose"], default="polytrope")
    parser.add_argument("--eos-table", default=None)
    parser.add_argument("--eos-name", default=None)
    parser.add_argument("--compose-dir", default=None)
    parser.add_argument("--compose-nb", default=None)
    parser.add_argument("--compose-thermo", default=None)
    parser.add_argument("--density-column", default=None)
    parser.add_argument("--pressure-column", default=None)
    parser.add_argument("--epsilon-column", default=None)
    parser.add_argument("--central-pressure", type=float, default=4.0e-4)
    parser.add_argument("--central-density", type=float, default=None)
    parser.add_argument("--kappa", type=float, default=100.0)
    parser.add_argument("--gamma", type=float, default=2.0)
    parser.add_argument("--dr", type=float, default=4.0e-3)
    parser.add_argument("--curve-deltas", default="60,75,135")
    parser.add_argument("--delta-min", type=float, default=10.0)
    parser.add_argument("--delta-max", type=float, default=170.0)
    parser.add_argument("--delta-count", type=int, default=65)
    parser.add_argument("--compactness-delta", type=float, default=90.0)
    parser.add_argument("--sequence-pc", default="5e-5,1e-4,2e-4,4e-4,8e-4,1.2e-3,2e-3")
    parser.add_argument("--sequence-density", default=None)
    parser.add_argument("--multidelta-deltas", default="45,90,135",
                        help="comma-separated Delta values for multi-delta scan (fig_09)")
    parser.add_argument("--depth-diff-deltas", default="30,60,90,120,150",
                        help="comma-separated Delta values for depth difference (fig_08)")
    parser.add_argument("--gamma-values", default="2.0,2.5,3.0",
                        help="comma-separated Gamma values for polytropic EOS scan")
    parser.add_argument("--gamma-delta", type=float, default=135.0,
                        help="aperture angle in degrees for gamma scan figures")
    parser.add_argument("--compose-delta-values", default="30,60,90,120,150",
                        help="comma-separated Delta values for CompOSE curve figures")
    parser.add_argument("--splitting-delta", type=float, default=135.0,
                        help="aperture angle for t-tau splitting figure (fig_14)")
    parser.add_argument("--splitting-gammas", default="2.0,3.0,10.0",
                        help="comma-separated Gamma values for splitting figure EOS comparison")
    parser.add_argument("--splitting-mu",
                        default="0.02,0.025,0.03,0.035,0.06,0.08,0.10,0.12,0.14,0.16,0.18,0.20",
                        help="comma-separated target mu values for splitting scan")
    parser.add_argument("--no-weakfield", action="store_true",
                        help="skip the slow O(mu) weak-field formula computation")
    parser.add_argument("--schw-gammas", default="5.0,10.0,20.0",
                        help="Gamma values for Schwarzschild validation figure")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    outdir = os.path.abspath(args.outdir)
    datadir = os.path.abspath(args.datadir) if args.datadir else os.path.join(outdir, "data")
    ensure_dir(outdir)
    ensure_dir(datadir)

    eos = build_eos(args)
    central_pressure = central_pressure_from_args(args, eos)
    args.central_pressure = central_pressure
    sol = integrate_tov(central_pressure, eos=eos, dr=args.dr)
    deltas_grid = np.linspace(args.delta_min, args.delta_max, args.delta_count)
    curve_deltas = parse_float_list(args.curve_deltas)

    if args.figure in ("all", "profiles"):
        make_profiles(sol, outdir, datadir)
    if args.figure in ("all", "indices"):
        make_indices(sol, outdir, datadir)
    if args.figure in ("all", "curves"):
        make_curves(sol, outdir, datadir, curve_deltas)
    scan = None
    if args.figure in ("all", "turning", "times"):
        scan = make_turning(sol, outdir, datadir, deltas_grid)
    if args.figure in ("all", "times"):
        if scan is None:
            scan = aperture_scan(sol, deltas_grid)
        make_times(scan, outdir, datadir)
    if args.figure in ("all", "compactness"):
        make_compactness_scan(
            eos,
            sequence_pressures_from_args(args, eos),
            args.dr,
            args.compactness_delta,
            outdir,
            datadir,
        )

    multidelta_deltas = parse_float_list(args.multidelta_deltas)
    depth_diff_deltas = parse_float_list(args.depth_diff_deltas)

    gamma_values = parse_float_list(args.gamma_values)
    compose_delta_values = parse_float_list(args.compose_delta_values)

    def _load_compose_eos():
        nb = args.compose_nb or os.path.join(default_compose_dir(), "eos.nb.ns")
        th = args.compose_thermo or os.path.join(default_compose_dir(), "eos.thermo.ns")
        return TabulatedEOS.from_compose_ns(nb, th, name=args.eos_name or "qmc_rmf_1")

    if args.figure in ("all", "comparison"):
        eos_alt = None
        if args.eos == "polytrope":
            try:
                eos_alt = _load_compose_eos()
            except Exception as exc:
                print("EOS comparison: could not load CompOSE EOS (%s)." % exc)
        elif args.eos in ("compose", "table"):
            eos_alt = PolytropicEOS(kappa=args.kappa, gamma=args.gamma)
        if eos_alt is not None:
            make_eos_comparison(sol, eos_alt, args.dr, outdir, datadir)
        elif args.figure == "comparison":
            print("EOS comparison requires a paired EOS. Use --eos polytrope with CompOSE data present.")

    if args.figure in ("all", "depth_diff"):
        make_depth_difference(
            eos, sequence_pressures_from_args(args, eos),
            args.dr, depth_diff_deltas, outdir, datadir,
        )

    if args.figure in ("all", "multidelta"):
        make_compactness_scan_multidelta(
            eos, sequence_pressures_from_args(args, eos),
            args.dr, multidelta_deltas, outdir, datadir,
        )

    if args.figure in ("all", "gamma_t"):
        make_gamma_scan(central_pressure, args.dr, args.gamma_delta, gamma_values,
                        "coordinate", outdir, datadir, kappa=args.kappa)

    if args.figure in ("all", "gamma_tau"):
        make_gamma_scan(central_pressure, args.dr, args.gamma_delta, gamma_values,
                        "proper", outdir, datadir, kappa=args.kappa)

    if args.figure in ("all", "compose_curves"):
        eos_compose = None
        if args.eos == "polytrope":
            try:
                eos_compose = _load_compose_eos()
            except Exception as exc:
                print("compose_curves: could not load CompOSE EOS (%s)." % exc)
        elif args.eos in ("compose", "table"):
            eos_compose = eos
        if eos_compose is not None:
            make_compose_curves(eos_compose, args.dr, sol.compactness,
                                compose_delta_values, "coordinate", outdir, datadir)
            make_compose_curves(eos_compose, args.dr, sol.compactness,
                                compose_delta_values, "proper", outdir, datadir)

    if args.figure in ("all", "splitting"):
        splitting_gammas = parse_float_list(args.splitting_gammas)
        splitting_mu     = parse_float_list(args.splitting_mu)
        eos_pairs: List[Tuple[str, Any]] = [
            ("\\Gamma=%.1f" % g, PolytropicEOS(kappa=args.kappa, gamma=g))
            for g in splitting_gammas
        ]
        make_t_tau_splitting(
            eos_pairs, args.dr, args.splitting_delta,
            outdir, datadir,
            mu_targets=splitting_mu,
            compute_weakfield=not args.no_weakfield,
        )

    if args.figure in ("all", "schwarzschild_val"):
        schw_gammas = parse_float_list(args.schw_gammas)
        schw_mu = list(np.linspace(0.04, 0.35, 20))
        make_schwarzschild_validation(
            schw_gammas, args.dr, args.splitting_delta,
            outdir, datadir, mu_targets=schw_mu,
        )

    write_readme(outdir, datadir, sol, args)

    print(f"Wrote figures to: {outdir}")
    print(f"Wrote CSV data to: {datadir}")
    print(f"EOS: {sol.eos.name}")
    print(f"Model: R={sol.radius:.6g}, M={sol.total_mass:.6g}, mu={sol.compactness:.6g}, A(R)={sol.surface_lapse:.6g}")
    if not HAS_REPORTLAB:
        print("reportlab is not available; PDF export was skipped.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
