#!/usr/bin/env python3
"""
fig_indice_rifrazione_nr_vs_mu.py

Grafico di n(r) per la brachistocrona-t (metrica interna di Schwarzschild),
valutato a r = 0.1 R, in funzione della compattezza mu in [0, 8/9).

Definizione usata (coerente con il funzionale del tempo coordinato):
    n_t(x, mu) = 1 / ( A(x,mu) * sqrt(1 - A(x,mu)^2 / A(1,mu)^2) )
con:
    A(x,mu) = 1/2 * (3 sqrt(1-mu) - sqrt(1-mu x^2))
    x = r/R
"""

import argparse
import math
import os
from typing import Tuple

import numpy as np
from PIL import Image, ImageChops, ImageDraw, ImageFont


MU_MAX = 8.0 / 9.0


def A_lapse(x, mu):
    return 0.5 * (3.0 * np.sqrt(1.0 - mu) - np.sqrt(1.0 - mu * x * x))


def n_t_brachistocrona(x: float, mu: np.ndarray) -> np.ndarray:
    """
    Indice effettivo del problema di minimo tempo coordinato.
    """
    e = np.sqrt(1.0 - mu)  # A(1,mu)
    a = A_lapse(x, mu)
    rad = np.maximum(1.0 - (a * a) / (e * e), 1e-300)
    return 1.0 / (a * np.sqrt(rad))


def campiona_mu(mu_eps: float) -> np.ndarray:
    """
    Campionamento denso vicino ai bordi mu->0 e mu->8/9.
    """
    left = np.geomspace(mu_eps, 3e-2, 260, endpoint=False)
    mid = np.linspace(3e-2, 8.2e-1, 360, endpoint=False)
    right_span = MU_MAX - 8.2e-1
    right = MU_MAX - np.geomspace(mu_eps, right_span, 280)
    mu = np.concatenate([left, mid, right])
    mu = np.unique(np.clip(mu, mu_eps, MU_MAX - mu_eps))
    mu.sort()
    return mu


def _load_font(size: int):
    for name in ["Arial.ttf", "Helvetica.ttc", "DejaVuSans.ttf"]:
        try:
            return ImageFont.truetype(name, size=size)
        except Exception:
            continue
    return ImageFont.load_default()


def _ticks_log(ymin: float, ymax: float) -> list[float]:
    out = []
    dec_min = int(math.floor(math.log10(ymin)))
    dec_max = int(math.ceil(math.log10(ymax)))
    for d in range(dec_min, dec_max + 1):
        base = 10.0**d
        for m in [1, 2, 5]:
            v = m * base
            if ymin <= v <= ymax:
                out.append(v)
    return out


def draw_plot(out_file: str, mu: np.ndarray, nvals: np.ndarray, x_fixed: float) -> None:
    ss = 2
    W, H = 1800 * ss, 1200 * ss
    img = Image.new("RGB", (W, H), (255, 255, 255))
    d = ImageDraw.Draw(img)

    left = 150 * ss
    right = 120 * ss
    top = 85 * ss
    bottom = 140 * ss
    x0, y0 = left, top
    x1, y1 = W - right, H - bottom

    # axes
    d.rectangle((x0, y0, x1, y1), outline=(0, 0, 0), width=2 * ss)

    f_tick = _load_font(20 * ss)
    f_lab = _load_font(25 * ss)
    f_title = _load_font(30 * ss)

    mu_min = float(mu[0])
    mu_max = float(mu[-1])
    y_min = float(np.min(nvals))
    y_max = float(np.max(nvals))
    ly_min = math.log10(y_min)
    ly_max = math.log10(y_max)

    def map_x(u: float) -> float:
        return x0 + (u - 0.0) / (MU_MAX - 0.0) * (x1 - x0)

    def map_y(v: float) -> float:
        lv = math.log10(v)
        t = (lv - ly_min) / (ly_max - ly_min)
        return y1 - t * (y1 - y0)

    # grid + ticks x
    xticks = [0.0, 0.2, 0.4, 0.6, 0.8, MU_MAX]
    xlabels = ["0", "0.2", "0.4", "0.6", "0.8", "8/9"]
    for xv, lbl in zip(xticks, xlabels):
        xx = map_x(xv)
        d.line((xx, y1, xx, y1 + 10 * ss), fill=(0, 0, 0), width=2 * ss)
        d.line((xx, y0, xx, y1), fill=(228, 228, 228), width=1 * ss)
        tb = d.textbbox((0, 0), lbl, font=f_tick)
        tw = tb[2] - tb[0]
        d.text((xx - tw / 2, y1 + 18 * ss), lbl, fill=(0, 0, 0), font=f_tick)

    # grid + ticks y (log scale)
    for yv in _ticks_log(y_min, y_max):
        yy = map_y(yv)
        d.line((x0 - 10 * ss, yy, x0, yy), fill=(0, 0, 0), width=2 * ss)
        d.line((x0, yy, x1, yy), fill=(232, 232, 232), width=1 * ss)
        lbl = f"{yv:g}"
        tb = d.textbbox((0, 0), lbl, font=f_tick)
        tw = tb[2] - tb[0]
        th = tb[3] - tb[1]
        d.text((x0 - 18 * ss - tw, yy - th / 2), lbl, fill=(0, 0, 0), font=f_tick)

    # curve
    pts = [(map_x(float(u)), map_y(float(v))) for u, v in zip(mu, nvals)]
    d.line(pts, fill=(33, 98, 194), width=3 * ss)

    # labels
    xlabel = "compactness μ"
    ylabel = f"n_t(r = {x_fixed:.1f} R)"
    title = "Effective refractive index for t-brachistochrone"
    subtitle = "(log scale on y axis)"

    tbx = d.textbbox((0, 0), xlabel, font=f_lab)
    d.text(((x0 + x1) / 2 - (tbx[2] - tbx[0]) / 2, H - 70 * ss), xlabel, fill=(0, 0, 0), font=f_lab)

    # y label (simple horizontal near top-left for portability)
    d.text((22 * ss, y0 - 6 * ss), ylabel, fill=(0, 0, 0), font=f_lab)

    d.text((x0, 18 * ss), title, fill=(0, 0, 0), font=f_title)
    d.text((x0, 52 * ss), subtitle, fill=(60, 60, 60), font=f_tick)

    # mark asymptotic boundaries
    x_left = map_x(mu_min)
    x_right = map_x(mu_max)
    d.line((x_left, y0, x_left, y1), fill=(160, 160, 160), width=1 * ss)
    d.line((x_right, y0, x_right, y1), fill=(160, 160, 160), width=1 * ss)

    # downsample + trim
    lanczos = getattr(Image, "Resampling", Image).LANCZOS
    out = img.resize((W // ss, H // ss), lanczos)
    bg = Image.new("RGB", out.size, (255, 255, 255))
    diff = ImageChops.difference(out, bg)
    bbox = diff.getbbox()
    if bbox:
        pad = 10
        l = max(0, bbox[0] - pad)
        t = max(0, bbox[1] - pad)
        r = min(out.size[0], bbox[2] + pad)
        b = min(out.size[1], bbox[3] + pad)
        out = out.crop((l, t, r, b))

    os.makedirs(os.path.dirname(out_file), exist_ok=True)
    out.save(out_file, format="PNG")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--x-r-over-R", type=float, default=0.1, help="fixed x = r/R value")
    parser.add_argument("--mu-eps", type=float, default=1e-6,
                        help="distance from endpoints 0 and 8/9")
    parser.add_argument(
        "--out-file",
        type=str,
        default="figures/fig_n_of_mu_at_r01R.png",
        help="output PNG path",
    )
    args = parser.parse_args()

    x = args.x_r_over_R
    if not (0.0 < x < 1.0):
        raise ValueError("x-r-over-R must be in (0,1).")
    if not (0.0 < args.mu_eps < 1e-2):
        raise ValueError("mu-eps should be in (0, 1e-2).")

    mu = campiona_mu(args.mu_eps)
    nvals = n_t_brachistocrona(x=x, mu=mu)

    out_file = args.out_file
    if not os.path.isabs(out_file):
        out_file = os.path.abspath(out_file)

    draw_plot(out_file=out_file, mu=mu, nvals=nvals, x_fixed=x)

    i_min = int(np.argmin(nvals))
    print()
    print("=== n_t(r) vs mu ===")
    print(f"x = r/R                 : {x:.6f}")
    print(f"mu interval             : [{mu[0]:.3e}, {mu[-1]:.10f}]")
    print(f"n_t min                 : {nvals[i_min]:.6f} at mu={mu[i_min]:.6f}")
    print(f"n_t(mu->0+) approx      : {nvals[0]:.6f}")
    print(f"n_t(mu->8/9-) approx    : {nvals[-1]:.6f}")
    print(f"Saved figure            : {out_file}")


if __name__ == "__main__":
    main()
