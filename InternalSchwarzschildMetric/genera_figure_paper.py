#!/usr/bin/env python3
"""
genera_figure_paper.py

Genera le illustrazioni principali per il paper sul tunnel brachistocrono:

1) Curve q(Delta) newtoniana / NLO / NNLO (espansione in mu) per una lista di mu.
2) Confronto della serie (Newton, NLO, NNLO) con la soluzione numerica "esatta"
   della brachistocrona relativistica interna (Schwarzschild interna).
3) Grafico dell'errore relativo NLO/NNLO rispetto alla soluzione numerica.

Output:
    Cartella figure (default: ./figures), con PNG + PDF + CSV di supporto.
"""

import argparse
import csv
import math
import os
from typing import List, Tuple

import numpy as np

from simula_brachistocrona_schwarzschild_interna import scegli_q_ottimo


def parse_mu_list(text: str) -> List[float]:
    vals = []
    for raw in text.split(","):
        s = raw.strip()
        if not s:
            continue
        vals.append(float(s))
    if not vals:
        raise ValueError("mu list is empty.")
    return vals


def mu_tag(mu: float) -> str:
    return f"{mu:.4f}".replace(".", "p")


def q_series(delta_rad: np.ndarray, mu: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Serie perturbativa:
        q = q0 + mu*q1 + mu^2*q2 + O(mu^3)
    dove q1 e q2 sono qui inclusi direttamente in q_NLO e q_NNLO.
    """
    q0 = 1.0 - delta_rad / np.pi
    q_nlo = q0 + mu * q0 * (15.0 - 7.0 * q0**2) / 32.0
    q_nnlo = (
        q_nlo
        + mu**2
        * q0
        * (1827.0 - 3489.0 * q0**2 + 2457.0 * q0**4 - 539.0 * q0**6)
        / 4096.0
    )
    return q0, q_nlo, q_nnlo


def curva_da_q_newton_riancorata(delta_rad: float, q: float, n: int = 3000) -> Tuple[np.ndarray, np.ndarray]:
    """
    Curva tipo brachistocrona newtoniana parametrizzata da q, riancorata per
    avere estremi separati esattamente da Delta target.

    Nota:
      La formula chiusa newtoniana impone naturalmente Delta=(1-q)pi.
      Se q arriva da espansione NLO/NNLO a Delta fissato, puo' non rispettare
      esattamente questa relazione; quindi si applica una piccola riscalatura
      angolare per imporre gli estremi desiderati.
    """
    q = float(q)
    q = max(0.0, min(0.999999999, q))

    if n < 1000:
        n = 1000

    nh = n // 2
    u = np.linspace(0.0, math.pi / 2.0, nh)

    rho = np.sqrt(q * q + (1.0 - q * q) * np.sin(u) ** 2)
    phi = np.arctan2(np.sin(u), q * np.cos(u)) - q * u

    rho_left = rho[::-1]
    phi_left = -phi[::-1]
    rho_right = rho[1:]
    phi_right = phi[1:]

    rho_full = np.concatenate([rho_left, rho_right])
    phi_full = np.concatenate([phi_left, phi_right])

    # Riancora sugli estremi target: phi_end = Delta/2
    phi_end = abs(phi_full[0])
    if phi_end > 1e-15:
        scale = (0.5 * delta_rad) / phi_end
        phi_full = phi_full * scale

    x = rho_full * np.cos(phi_full)
    y = rho_full * np.sin(phi_full)
    return x, y


def maybe_clip(q: np.ndarray, clip_physical: bool) -> np.ndarray:
    if clip_physical:
        return np.clip(q, 0.0, 1.0)
    return q


def build_exact_q(delta_deg: np.ndarray, mu: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Soluzione numerica "esatta" via risoluzione di Delta(q)=Delta_target.
    Restituisce:
      q_exact, err_delta_deg
    """
    q_vals = np.full_like(delta_deg, np.nan, dtype=float)
    err_vals = np.full_like(delta_deg, np.nan, dtype=float)

    for i, ddeg in enumerate(delta_deg):
        delta = math.radians(float(ddeg))
        try:
            q, _tbar, _delta_num, err, _cands = scegli_q_ottimo(delta_target=delta, u=mu)
            q_vals[i] = q
            err_vals[i] = math.degrees(err)
        except Exception:
            q_vals[i] = np.nan
            err_vals[i] = np.nan

    return q_vals, err_vals


def save_q_series_figure(
    delta_deg: np.ndarray,
    q0: np.ndarray,
    q_nlo: np.ndarray,
    q_nnlo: np.ndarray,
    mu: float,
    outdir: str,
    dpi: int,
    clip_physical: bool,
) -> None:
    import matplotlib.pyplot as plt

    q0p = maybe_clip(q0, clip_physical)
    q1p = maybe_clip(q_nlo, clip_physical)
    q2p = maybe_clip(q_nnlo, clip_physical)

    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    ax.plot(delta_deg, q0p, color="#1f77b4", linewidth=2.2, label="Newtonian")
    ax.plot(delta_deg, q1p, color="#d62728", linestyle="--", linewidth=2.0, label=f"NLO (mu={mu:g})")
    ax.plot(delta_deg, q2p, color="#2ca02c", linestyle="-.", linewidth=2.0, label=f"NNLO (mu={mu:g})")

    ax.set_xlabel("Delta [degrees]")
    ax.set_ylabel("q = r_min / R")
    ax.set_title(f"q(Delta) comparison: Newtonian, NLO, NNLO (mu={mu:g})")
    ax.grid(True, alpha=0.28)
    ax.legend()
    ax.set_xlim(float(delta_deg[0]), float(delta_deg[-1]))
    ax.set_ylim(-0.02, 1.02)
    fig.tight_layout()

    tag = mu_tag(mu)
    png = os.path.join(outdir, f"fig_q_vs_delta_series_mu_{tag}.png")
    pdf = os.path.join(outdir, f"fig_q_vs_delta_series_mu_{tag}.pdf")
    fig.savefig(png, dpi=dpi)
    fig.savefig(pdf)
    plt.close(fig)


def save_validation_figures(
    delta_dense_deg: np.ndarray,
    mu: float,
    q0_dense: np.ndarray,
    q_nlo_dense: np.ndarray,
    q_nnlo_dense: np.ndarray,
    delta_exact_deg: np.ndarray,
    q_exact: np.ndarray,
    outdir: str,
    dpi: int,
    clip_physical: bool,
) -> None:
    import matplotlib.pyplot as plt

    q0e, q1e, q2e = q_series(np.deg2rad(delta_exact_deg), mu)
    if clip_physical:
        q0e = np.clip(q0e, 0.0, 1.0)
        q1e = np.clip(q1e, 0.0, 1.0)
        q2e = np.clip(q2e, 0.0, 1.0)
        q0_dense = np.clip(q0_dense, 0.0, 1.0)
        q_nlo_dense = np.clip(q_nlo_dense, 0.0, 1.0)
        q_nnlo_dense = np.clip(q_nnlo_dense, 0.0, 1.0)

    # Figura 1: curve + punti numerici esatti
    fig1, ax1 = plt.subplots(figsize=(8.2, 5.2))
    ax1.plot(delta_dense_deg, q0_dense, color="#1f77b4", linewidth=2.0, label="Newtonian")
    ax1.plot(delta_dense_deg, q_nlo_dense, color="#d62728", linestyle="--", linewidth=2.0, label="NLO")
    ax1.plot(delta_dense_deg, q_nnlo_dense, color="#2ca02c", linestyle="-.", linewidth=2.1, label="NNLO")
    ax1.scatter(delta_exact_deg, q_exact, color="#111111", s=18, alpha=0.9, label="Exact numerical")

    ax1.set_xlabel("Delta [degrees]")
    ax1.set_ylabel("q = r_min / R")
    ax1.set_title(f"Series validation vs numerical solution (mu={mu:g})")
    ax1.grid(True, alpha=0.28)
    ax1.legend()
    ax1.set_xlim(float(delta_dense_deg[0]), float(delta_dense_deg[-1]))
    ax1.set_ylim(-0.02, 1.02)
    fig1.tight_layout()

    tag = mu_tag(mu)
    fig1_png = os.path.join(outdir, f"fig_validazione_q_vs_delta_mu_{tag}.png")
    fig1_pdf = os.path.join(outdir, f"fig_validazione_q_vs_delta_mu_{tag}.pdf")
    fig1.savefig(fig1_png, dpi=dpi)
    fig1.savefig(fig1_pdf)
    plt.close(fig1)

    # Figura 2: errore relativo NLO e NNLO
    eps = 1e-12
    mask = np.isfinite(q_exact)
    rel_nlo = np.full_like(q_exact, np.nan)
    rel_nnlo = np.full_like(q_exact, np.nan)
    rel_nlo[mask] = np.abs(q1e[mask] - q_exact[mask]) / np.maximum(np.abs(q_exact[mask]), eps)
    rel_nnlo[mask] = np.abs(q2e[mask] - q_exact[mask]) / np.maximum(np.abs(q_exact[mask]), eps)

    fig2, ax2 = plt.subplots(figsize=(8.2, 5.2))
    ax2.plot(delta_exact_deg, 100.0 * rel_nlo, color="#d62728", linestyle="--", linewidth=1.9, label="NLO")
    ax2.plot(delta_exact_deg, 100.0 * rel_nnlo, color="#2ca02c", linestyle="-.", linewidth=2.1, label="NNLO")
    ax2.set_xlabel("Delta [degrees]")
    ax2.set_ylabel("Relative error on q [%]")
    ax2.set_title(f"Relative error: series vs numerical (mu={mu:g})")
    ax2.grid(True, alpha=0.28)
    ax2.legend()
    ax2.set_xlim(float(delta_exact_deg[0]), float(delta_exact_deg[-1]))
    fig2.tight_layout()

    fig2_png = os.path.join(outdir, f"fig_errore_relativo_mu_{tag}.png")
    fig2_pdf = os.path.join(outdir, f"fig_errore_relativo_mu_{tag}.pdf")
    fig2.savefig(fig2_png, dpi=dpi)
    fig2.savefig(fig2_pdf)
    plt.close(fig2)

    # CSV di supporto per il paper
    csv_path = os.path.join(outdir, f"tab_validazione_mu_{tag}.csv")
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "Delta_deg",
                "q_exact_numeric",
                "q0_newton",
                "q_nlo",
                "q_nnlo",
                "relerr_nlo",
                "relerr_nnlo",
            ]
        )
        for i in range(len(delta_exact_deg)):
            w.writerow(
                [
                    float(delta_exact_deg[i]),
                    float(q_exact[i]) if np.isfinite(q_exact[i]) else "",
                    float(q0e[i]),
                    float(q1e[i]),
                    float(q2e[i]),
                    float(rel_nlo[i]) if np.isfinite(rel_nlo[i]) else "",
                    float(rel_nnlo[i]) if np.isfinite(rel_nnlo[i]) else "",
                ]
            )


def save_trajectory_figures(
    delta_deg: float,
    mu_ref: float,
    mu_min: float,
    mu_max: float,
    mu_count: int,
    outdir: str,
    dpi: int,
    clip_physical: bool,
    n_curve: int = 3200,
) -> None:
    import matplotlib.pyplot as plt

    if mu_count < 2:
        raise ValueError("mu_count must be at least 2.")

    delta_rad = math.radians(delta_deg)
    q0, _, _ = q_series(np.array([delta_rad]), 0.0)
    q0 = float(q0[0])

    # Figura A: Newton vs NLO a mu_ref
    _q0r, qnlo_ref_arr, _qnnlo_ref_arr = q_series(np.array([delta_rad]), mu_ref)
    qnlo_ref = float(qnlo_ref_arr[0])
    if clip_physical:
        qnlo_ref = float(np.clip(qnlo_ref, 0.0, 1.0))
        q0_plot = float(np.clip(q0, 0.0, 1.0))
    else:
        q0_plot = q0

    x_newt, y_newt = curva_da_q_newton_riancorata(delta_rad, q0_plot, n=n_curve)
    x_ref, y_ref = curva_da_q_newton_riancorata(delta_rad, qnlo_ref, n=n_curve)

    tt = np.linspace(0.0, 2.0 * math.pi, 1000)
    fig_a, ax_a = plt.subplots(figsize=(7.2, 7.2))
    ax_a.plot(np.cos(tt), np.sin(tt), "k-", linewidth=1.15, label="surface")
    ax_a.plot(x_newt, y_newt, color="#1f77b4", linewidth=2.5, label=f"Newtonian (q={q0_plot:.4f})")
    ax_a.plot(x_ref, y_ref, color="#d62728", linestyle="--", linewidth=2.4,
              label=f"NLO mu={mu_ref:g} (q={qnlo_ref:.4f})")
    ax_a.set_aspect("equal", adjustable="box")
    ax_a.set_xlim(-1.08, 1.08)
    ax_a.set_ylim(-1.08, 1.08)
    ax_a.set_xlabel("x/R")
    ax_a.set_ylabel("y/R")
    ax_a.set_title(f"Trajectory inside the sphere: Delta={delta_deg:.1f} deg")
    ax_a.grid(True, alpha=0.28)
    ax_a.legend(loc="upper right")
    fig_a.tight_layout()

    tag_ref = mu_tag(mu_ref)
    fa_png = os.path.join(outdir, f"fig_traiettoria_newton_vs_nlo_mu_{tag_ref}_delta_{delta_deg:.0f}deg.png")
    fa_pdf = os.path.join(outdir, f"fig_traiettoria_newton_vs_nlo_mu_{tag_ref}_delta_{delta_deg:.0f}deg.pdf")
    fig_a.savefig(fa_png, dpi=dpi)
    fig_a.savefig(fa_pdf)
    plt.close(fig_a)

    # Figura B: famiglia NLO per mu in [mu_min, mu_max]
    mu_vals = np.linspace(mu_min, mu_max, mu_count)
    q_rows = []

    fig_b, ax_b = plt.subplots(figsize=(7.5, 7.5))
    ax_b.plot(np.cos(tt), np.sin(tt), "k-", linewidth=1.15, label="surface")
    ax_b.plot(x_newt, y_newt, color="#1f77b4", linewidth=2.6, label=f"Newtonian (q={q0_plot:.4f})")

    cmap = plt.get_cmap("viridis")
    for i, mu in enumerate(mu_vals):
        _q0m, qnlo_m_arr, _qnnlo_m_arr = q_series(np.array([delta_rad]), float(mu))
        qnlo_m = float(qnlo_m_arr[0])
        if clip_physical:
            qnlo_m = float(np.clip(qnlo_m, 0.0, 1.0))
        x_m, y_m = curva_da_q_newton_riancorata(delta_rad, qnlo_m, n=n_curve)
        col = cmap(i / max(1, mu_count - 1))
        ax_b.plot(x_m, y_m, color=col, linewidth=1.8, alpha=0.95)
        q_rows.append((float(mu), qnlo_m))

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=mu_min, vmax=mu_max))
    sm.set_array([])
    cbar = fig_b.colorbar(sm, ax=ax_b, fraction=0.045, pad=0.03)
    cbar.set_label("mu (NLO)")

    ax_b.set_aspect("equal", adjustable="box")
    ax_b.set_xlim(-1.08, 1.08)
    ax_b.set_ylim(-1.08, 1.08)
    ax_b.set_xlabel("x/R")
    ax_b.set_ylabel("y/R")
    ax_b.set_title(
        f"NLO trajectory family inside the sphere, Delta={delta_deg:.1f} deg, mu in [{mu_min:.2f}, {mu_max:.2f}]"
    )
    ax_b.grid(True, alpha=0.28)
    ax_b.legend(loc="upper right")
    fig_b.tight_layout()

    ftag = f"{mu_min:.2f}_{mu_max:.2f}".replace(".", "p")
    fb_png = os.path.join(outdir, f"fig_famiglia_traiettorie_nlo_delta_{delta_deg:.0f}deg_mu_{ftag}.png")
    fb_pdf = os.path.join(outdir, f"fig_famiglia_traiettorie_nlo_delta_{delta_deg:.0f}deg_mu_{ftag}.pdf")
    fig_b.savefig(fb_png, dpi=dpi)
    fig_b.savefig(fb_pdf)
    plt.close(fig_b)

    # CSV di supporto per q_NLO(Delta fisso, mu variabile)
    csv_path = os.path.join(outdir, f"tab_qnlo_delta_{delta_deg:.0f}deg_mu_{ftag}.csv")
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["Delta_deg", "mu", "q_newton", "q_nlo"])
        for mu, qn in q_rows:
            w.writerow([delta_deg, mu, q0_plot, qn])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mu-list", type=str, default="0.05,0.10,0.30",
                        help="comma-separated mu values for series figures")
    parser.add_argument("--mu-validate", type=float, default=0.10,
                        help="mu value used for validation against numerical solution")
    parser.add_argument("--delta-min-deg", type=float, default=2.0,
                        help="minimum Delta in degrees (endpoint avoided)")
    parser.add_argument("--delta-max-deg", type=float, default=178.0,
                        help="maximum Delta in degrees (endpoint avoided)")
    parser.add_argument("--n-delta", type=int, default=500,
                        help="number of points for continuous curves")
    parser.add_argument("--n-exact", type=int, default=61,
                        help="number of exact numerical points for validation")
    parser.add_argument("--clip-physical", action="store_true",
                        help="clip q to physical range [0,1] in plots")
    parser.add_argument("--traj-delta-deg", type=float, default=90.0,
                        help="Delta [deg] for trajectory geometry figures")
    parser.add_argument("--traj-mu-ref", type=float, default=0.10,
                        help="reference mu for Newtonian vs NLO trajectory figure")
    parser.add_argument("--traj-mu-min", type=float, default=0.01,
                        help="minimum mu for NLO trajectory family")
    parser.add_argument("--traj-mu-max", type=float, default=0.30,
                        help="maximum mu for NLO trajectory family")
    parser.add_argument("--traj-mu-count", type=int, default=12,
                        help="number of NLO curves in family figure")
    parser.add_argument("--traj-n-points", type=int, default=3200,
                        help="number of points per curve in trajectory figures")
    parser.add_argument("--outdir", type=str, default="figures",
                        help="output folder")
    parser.add_argument("--dpi", type=int, default=220,
                        help="PNG resolution")
    args = parser.parse_args()

    if args.n_delta < 20:
        raise ValueError("n-delta must be >= 20.")
    if args.n_exact < 5:
        raise ValueError("n-exact must be >= 5.")
    if not (0.0 < args.delta_min_deg < args.delta_max_deg < 180.0):
        raise ValueError("Require 0 < delta-min-deg < delta-max-deg < 180.")
    if args.dpi < 100:
        raise ValueError("dpi too low: use at least 100.")
    if not (0.0 < args.traj_delta_deg < 180.0):
        raise ValueError("traj-delta-deg must lie in (0, 180).")
    if args.traj_mu_count < 2:
        raise ValueError("traj-mu-count must be >= 2.")
    if args.traj_n_points < 1000:
        raise ValueError("traj-n-points must be >= 1000.")
    if not (0.0 < args.traj_mu_min <= args.traj_mu_ref <= args.traj_mu_max < 8.0 / 9.0):
        raise ValueError("Require 0 < traj-mu-min <= traj-mu-ref <= traj-mu-max < 8/9.")

    mu_list = parse_mu_list(args.mu_list)
    for mu in mu_list + [args.mu_validate]:
        if not (0.0 < mu < 8.0 / 9.0):
            raise ValueError(f"mu={mu} is outside physical range (0, 8/9).")

    os.makedirs(args.outdir, exist_ok=True)

    delta_dense_deg = np.linspace(args.delta_min_deg, args.delta_max_deg, args.n_delta)
    delta_dense_rad = np.deg2rad(delta_dense_deg)

    print()
    print("=== Paper figure generation ===")
    print(f"outdir            : {args.outdir}")
    print(f"mu-list           : {mu_list}")
    print(f"mu-validate       : {args.mu_validate}")
    print(f"Delta range [deg] : [{args.delta_min_deg}, {args.delta_max_deg}]")
    print(f"n-delta           : {args.n_delta}")
    print(f"n-exact           : {args.n_exact}")
    print(f"clip-physical     : {args.clip_physical}")
    print(f"traj Delta [deg]  : {args.traj_delta_deg}")
    print(f"traj mu-ref       : {args.traj_mu_ref}")
    print(f"traj mu range     : [{args.traj_mu_min}, {args.traj_mu_max}] x {args.traj_mu_count}")

    # Figure serie per ogni mu richiesto
    for mu in mu_list:
        q0, q_nlo, q_nnlo = q_series(delta_dense_rad, mu)
        save_q_series_figure(
            delta_deg=delta_dense_deg,
            q0=q0,
            q_nlo=q_nlo,
            q_nnlo=q_nnlo,
            mu=mu,
            outdir=args.outdir,
            dpi=args.dpi,
            clip_physical=args.clip_physical,
        )
        print(f"Saved series figure for mu={mu:g}")

    # Validazione numerica su mu specificato
    mu_v = args.mu_validate
    q0_v, q1_v, q2_v = q_series(delta_dense_rad, mu_v)

    delta_exact_deg = np.linspace(args.delta_min_deg, args.delta_max_deg, args.n_exact)
    q_exact, err_delta_deg = build_exact_q(delta_exact_deg, mu_v)
    n_ok = int(np.sum(np.isfinite(q_exact)))
    if n_ok == 0:
        raise RuntimeError("Validation failed: no converged exact numerical points.")

    save_validation_figures(
        delta_dense_deg=delta_dense_deg,
        mu=mu_v,
        q0_dense=q0_v,
        q_nlo_dense=q1_v,
        q_nnlo_dense=q2_v,
        delta_exact_deg=delta_exact_deg,
        q_exact=q_exact,
        outdir=args.outdir,
        dpi=args.dpi,
        clip_physical=args.clip_physical,
    )

    err_abs = np.abs(err_delta_deg[np.isfinite(err_delta_deg)])
    if len(err_abs) > 0:
        print(f"Validation: {n_ok}/{args.n_exact} converged points")
        print(f"Mean Delta error [deg]: {np.mean(err_abs):.3e}")
        print(f"Max Delta error [deg] : {np.max(err_abs):.3e}")

    # Figure geometriche richieste per il paper
    save_trajectory_figures(
        delta_deg=args.traj_delta_deg,
        mu_ref=args.traj_mu_ref,
        mu_min=args.traj_mu_min,
        mu_max=args.traj_mu_max,
        mu_count=args.traj_mu_count,
        outdir=args.outdir,
        dpi=args.dpi,
        clip_physical=args.clip_physical,
        n_curve=args.traj_n_points,
    )
    print("Trajectory geometry figures: completed")

    print("Completato.")


if __name__ == "__main__":
    main()
