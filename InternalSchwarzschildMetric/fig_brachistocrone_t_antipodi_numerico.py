#!/usr/bin/env python3
"""
fig_brachistocrone_t_antipodi_numerico.py

Figura: famiglia di brachistocrone t (tempo coordinato) nella metrica interna
di Schwarzschild, con estremi ai poli Nord/Sud (Delta = pi), al variare della
compattezza mu da 0.1 verso 8/9.

Importante:
  - usa integrazione numerica completa (nessuna espansione post-newtoniana);
  - per ogni mu sceglie la soluzione q con tempo coordinato minimo tra le radici
    di Delta(q)=pi trovate numericamente;
  - quando il root-finding non intercetta un cambio di segno (caso limite q~0),
    usa la soluzione numerica con Delta massimo nel campionamento (q piccolissimo).
"""

import argparse
import colorsys
import math
import os
from typing import List, Tuple

import numpy as np
from PIL import Image, ImageChops, ImageDraw, ImageFont

from simula_brachistocrona_schwarzschild_interna import (
    U_MAX,
    delta_da_q,
    profilo_phi_mezza_curva,
    tempo_coord_adim_da_q,
)


def _safe_delta_minus_target(q: float, u: float, delta_target: float, n_int: int) -> float:
    try:
        return float(delta_da_q(q, u, n=n_int) - delta_target)
    except Exception:
        return float("nan")


def trova_radici_delta_pi(u: float, n_int_scan: int = 1800) -> Tuple[List[float], float, float]:
    """
    Restituisce:
      - lista radici q di Delta(q)=pi
      - q_peak dove Delta(q) e' massimo nel campionamento
      - Delta_peak (radianti)
    """
    delta_target = math.pi

    q_small = np.geomspace(1e-7, 5e-3, 280, endpoint=False)
    q_large = np.linspace(5e-3, 0.999, 520)
    q_grid = np.concatenate([q_small, q_large])

    vals = np.array([_safe_delta_minus_target(float(q), u, delta_target, n_int_scan) for q in q_grid])
    delta_vals = vals + delta_target

    valid = np.isfinite(delta_vals)
    if not np.any(valid):
        return [], 1e-6, float("nan")

    i_peak = int(np.nanargmax(delta_vals))
    q_peak = float(q_grid[i_peak])
    delta_peak = float(delta_vals[i_peak])

    roots = []
    for i in range(len(q_grid) - 1):
        f1 = vals[i]
        f2 = vals[i + 1]
        if not (np.isfinite(f1) and np.isfinite(f2)):
            continue

        if abs(f1) < 1e-8:
            roots.append(float(q_grid[i]))
            continue

        if f1 * f2 < 0.0:
            a = float(q_grid[i])
            b = float(q_grid[i + 1])
            fa = f1
            for _ in range(62):
                m = 0.5 * (a + b)
                fm = _safe_delta_minus_target(m, u, delta_target, n_int_scan + 800)
                if not np.isfinite(fm):
                    break
                if fa * fm <= 0.0:
                    b = m
                else:
                    a = m
                    fa = fm
            roots.append(0.5 * (a + b))

    roots_sorted = []
    for q in sorted(roots):
        if (not roots_sorted) or abs(q - roots_sorted[-1]) > 2e-5:
            roots_sorted.append(float(q))

    return roots_sorted, q_peak, delta_peak


def seleziona_q_brachistocrona_t(u: float) -> Tuple[float, str, float]:
    """
    Seleziona q della brachistocrona t per Delta=pi:
      - root minimo-tempo se ci sono radici;
      - fallback al q di Delta massimo (in pratica q~0 per mu piccoli).

    Restituisce: q, mode, delta_peak_deg
      mode in {"root-min-time", "peak-fallback"}
    """
    if abs(u) < 1e-14:
        # Newtonian antipodal limit: straight diameter.
        return 0.0, "newtonian-diameter", 180.0

    roots, q_peak, delta_peak = trova_radici_delta_pi(u)

    if roots:
        best_q = None
        best_t = None
        for q in roots:
            try:
                tbar = tempo_coord_adim_da_q(q, u, n=5500)
            except Exception:
                continue
            if (best_t is None) or (tbar < best_t):
                best_t = tbar
                best_q = q
        if best_q is not None:
            return float(best_q), "root-min-time", math.degrees(delta_peak)

    return float(q_peak), "peak-fallback", math.degrees(delta_peak)


def traiettoria_da_q_antipodi(u: float, q: float, n: int = 2600) -> Tuple[np.ndarray, np.ndarray]:
    """
    Costruisce la traiettoria per Delta=pi a partire da q usando il profilo
    numerico phi(x) e una piccola rinormalizzazione angolare ai poli.
    """
    if abs(u) < 1e-14 or q <= 1e-12:
        y = np.linspace(1.0, -1.0, max(1200, n))
        x = np.zeros_like(y)
        return x, y

    delta_target = math.pi
    x_half, phi_half = profilo_phi_mezza_curva(q, u, n=n)

    if abs(phi_half[-1]) > 1e-14:
        scale = (0.5 * delta_target) / phi_half[-1]
        phi_half = phi_half * scale

    x_left = x_half[::-1]
    phi_left = -phi_half[::-1]
    x_right = x_half[1:]
    phi_right = phi_half[1:]

    x_rad = np.concatenate([x_left, x_right])
    phi = np.concatenate([phi_left, phi_right])

    x = x_rad * np.cos(phi)
    y = x_rad * np.sin(phi)
    return x, y


def _make_color(i: int, n: int) -> Tuple[int, int, int]:
    # blu -> rosso al crescere di mu
    t = i / max(1, n - 1)
    hue = 0.66 * (1.0 - t)  # 0.66~blu, 0.0~rosso
    sat = 0.82
    val = 0.92
    r, g, b = colorsys.hsv_to_rgb(hue, sat, val)
    return int(255 * r), int(255 * g), int(255 * b)


def _load_font(size: int) -> ImageFont.ImageFont:
    for name in ["Arial.ttf", "Helvetica.ttc", "DejaVuSans.ttf"]:
        try:
            return ImageFont.truetype(name, size=size)
        except Exception:
            continue
    return ImageFont.load_default()


def _to_px(x: float, y: float, cx: float, cy: float, rpx: float) -> Tuple[float, float]:
    return cx + rpx * x, cy - rpx * y


def draw_figure(
    out_file: str,
    u_values: List[float],
    curve_points: List[Tuple[np.ndarray, np.ndarray]],
    q_values: List[float],
    modes: List[str],
) -> None:
    ss = 2  # supersampling per linee piu' pulite
    W, H = 2200 * ss, 1600 * ss
    img = Image.new("RGB", (W, H), (255, 255, 255))
    d = ImageDraw.Draw(img)

    cx, cy = 860 * ss, 800 * ss
    rpx = 660 * ss

    # Sphere + N-S comparison line
    d.ellipse((cx - rpx, cy - rpx, cx + rpx, cy + rpx), outline=(0, 0, 0), width=4 * ss)
    d.line((cx, cy - rpx, cx, cy + rpx), fill=(110, 110, 110), width=3 * ss)

    # Curves
    n = len(u_values)
    for i, (u, (x, y), q) in enumerate(zip(u_values, curve_points, q_values)):
        col = _make_color(i, n)
        pts = [_to_px(float(xx), float(yy), cx, cy, rpx) for xx, yy in zip(x, y)]
        d.line(pts, fill=col, width=4 * ss)

    # Poli
    np_px = _to_px(0.0, 1.0, cx, cy, rpx)
    sp_px = _to_px(0.0, -1.0, cx, cy, rpx)
    rr = 8 * ss
    d.ellipse((np_px[0] - rr, np_px[1] - rr, np_px[0] + rr, np_px[1] + rr), fill=(0, 0, 0))
    d.ellipse((sp_px[0] - rr, sp_px[1] - rr, sp_px[0] + rr, sp_px[1] + rr), fill=(0, 0, 0))

    f_sub = _load_font(26 * ss)
    f_leg = _load_font(25 * ss)

    # Pole labels: N slightly above the circle, S below the circle
    d.text((int(np_px[0] + 14 * ss), int(np_px[1] - 58 * ss)), "N", fill=(20, 20, 20), font=f_sub)
    d.text((int(sp_px[0] + 14 * ss), int(sp_px[1] + 12 * ss)), "S", fill=(20, 20, 20), font=f_sub)

    # Legend (mu only)
    lx = 1530 * ss
    ly = 210 * ss
    d.text((lx, ly - 70 * ss), "μ and corresponding q (Buchdahl: μ < 8/9):", fill=(20, 20, 20), font=f_sub)
    row_h = 44 * ss
    for i, (u, q) in enumerate(zip(u_values, q_values)):
        y0 = ly + i * row_h
        col = _make_color(i, n)
        d.line((lx, y0 + 12 * ss, lx + 72 * ss, y0 + 12 * ss), fill=col, width=5 * ss)
        txt = f"μ = {u:.4f}   q = {q:.4f}"
        d.text((lx + 88 * ss, y0), txt, fill=(20, 20, 20), font=f_leg)

    # downsample finale
    lanczos = getattr(Image, "Resampling", Image).LANCZOS
    img_out = img.resize((W // ss, H // ss), lanczos)

    # Trim external white margins to keep the figure compact for publication.
    bg = Image.new("RGB", img_out.size, (255, 255, 255))
    diff = ImageChops.difference(img_out, bg)
    bbox = diff.getbbox()
    if bbox is not None:
        pad_left = 16
        pad_top = 16
        pad_bottom = 16
        pad_right = 2
        left = max(0, bbox[0] - pad_left)
        top = max(0, bbox[1] - pad_top)
        right = min(img_out.size[0], bbox[2] + pad_right)
        bottom = min(img_out.size[1], bbox[3] + pad_bottom)
        img_out = img_out.crop((left, top, right, bottom))
    os.makedirs(os.path.dirname(out_file), exist_ok=True)
    img_out.save(out_file, format="PNG")


def parse_u_list(text: str) -> List[float]:
    vals = []
    for raw in text.split(","):
        s = raw.strip()
        if not s:
            continue
        vals.append(float(s))
    if not vals:
        raise ValueError("u list is empty.")
    return vals


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--u-list",
        type=str,
        default="0.0,0.78,0.82,0.85,0.87,0.88,0.885,0.888,0.8885,0.8888,0.8888878889",
        help="comma-separated compactness values mu in (0, 8/9)",
    )
    parser.add_argument(
        "--out-file",
        type=str,
        default="figures/fig_brachistocrone_t_antipodi_compattezza.png",
        help="output PNG path",
    )
    parser.add_argument("--n-curve", type=int, default=2600, help="samples per trajectory")
    args = parser.parse_args()

    u_values = parse_u_list(args.u_list)
    for u in u_values:
        if not (0.0 <= u < U_MAX):
            raise ValueError(f"u={u} is outside physical range [0, 8/9).")
    if args.n_curve < 1200:
        raise ValueError("n-curve must be >= 1200.")

    curve_points = []
    q_values = []
    modes = []

    print()
    print("=== Antipodal t-brachistochrone family (numerical integrals) ===")
    for u in u_values:
        q, mode, _dpeak_deg = seleziona_q_brachistocrona_t(u)
        x, y = traiettoria_da_q_antipodi(u=u, q=q, n=args.n_curve)
        curve_points.append((x, y))
        q_values.append(q)
        modes.append(mode)
        print(f"mu={u:.6f}  q={q:.7f}  mode={mode}")

    out_file = args.out_file
    if not os.path.isabs(out_file):
        out_file = os.path.abspath(out_file)

    draw_figure(
        out_file=out_file,
        u_values=u_values,
        curve_points=curve_points,
        q_values=q_values,
        modes=modes,
    )

    print()
    print(f"Saved figure: {out_file}")


if __name__ == "__main__":
    main()
