#!/usr/bin/env python3
"""
animazione_confronto_brachistocrone.py

Confronta nella stessa animazione:
  1) brachistocrona Newtoniana (sfera uniforme)
  2) brachistocrona relativistica (metrica di Schwarzschild interna)

I due moti partono dagli stessi punti di superficie (separazione Delta) e
vengono animati nello stesso asse temporale fisico T.
"""

import argparse
import math
import numpy as np

from simula_caduta_tunnel import curva_brachistocrona, profilo_temporale_curva, posizione_al_tempo
from simula_brachistocrona_schwarzschild_interna import (
    A_lapse,
    B_radiale,
    C,
    G,
    costruisci_traiettoria,
    compattezza_u,
    M_TERRA,
    R_TERRA,
    U_MAX,
)


def massa_da_u_R(u, R):
    return u * C * C * R / (2.0 * G)


def omega_uniforme_generica(M, R):
    return math.sqrt(G * M / (R ** 3))


def profilo_temporale_relativistico(x_rad, phi, u, R):
    """
    Restituisce T_cum (s) lungo la traiettoria relativistica gia' campionata.

    Formula:
      dt_bar = ds_bar / (A * v_loc)
      ds_bar^2 = B(x)^2 dx^2 + x^2 dphi^2
      v_loc = sqrt(1 - A(x)^2 / A(1)^2)
      T = (R/c) * t_bar
    """
    x_rad = np.asarray(x_rad)
    phi = np.asarray(phi)

    dx = np.diff(x_rad)
    dphi = np.diff(phi)

    xm = 0.5 * (x_rad[:-1] + x_rad[1:])
    Bm = B_radiale(xm, u)
    Am = A_lapse(xm, u)
    e = math.sqrt(1.0 - u)  # A(1)

    ds_bar = np.sqrt((Bm * dx) ** 2 + (xm * dphi) ** 2)
    v_loc = np.sqrt(np.maximum(1.0 - (Am * Am) / (e * e), 1e-18))
    dt_bar = ds_bar / np.maximum(Am * v_loc, 1e-18)

    tbar_cum = np.concatenate(([0.0], np.cumsum(dt_bar)))
    T_cum = (R / C) * tbar_cum
    return T_cum


def salva_animazione_confronto(
    delta_deg,
    u,
    R,
    out_file,
    fps=18,
    anim_dur=10.0,
    N=6000,
):
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation

    delta = math.radians(delta_deg)

    # Newton: stessa stella (u, R) -> M coerente
    M = massa_da_u_R(u, R)
    omega = omega_uniforme_generica(M, R)
    xN, yN = curva_brachistocrona(delta, N=N)
    tauN, TN, _tau_cum_N, T_cum_N = profilo_temporale_curva(xN, yN, omega)

    # Relativistica interna
    risR = costruisci_traiettoria(delta_target=delta, u=u, n=max(2200, N // 2))
    xR = risR["x"]
    yR = risR["y"]
    xradR = risR["x_rad"]
    phiR = risR["phi"]
    T_cum_R = profilo_temporale_relativistico(xradR, phiR, u=u, R=R)
    TR = T_cum_R[-1]

    T_max = max(TN, TR)
    T_anim = 1.12 * T_max
    n_frames = max(2, int(fps * anim_dur))
    t_frames = np.linspace(0.0, T_anim, n_frames)

    fig, ax = plt.subplots(figsize=(8, 8))
    t = np.linspace(0.0, 2.0 * math.pi, 1000)
    ax.plot(np.cos(t), np.sin(t), "k-", linewidth=1.1)

    colN = "#c9a227"
    colR = "#8f63d8"
    ax.plot(xN, yN, color=colN, linewidth=2.8, label=f"Newtoniana ({TN/60.0:.2f} min)")
    ax.plot(xR, yR, color=colR, linewidth=2.8, label=f"Relativistica ({TR/60.0:.2f} min)")

    pN, = ax.plot([xN[0]], [yN[0]], "o", color=colN, markersize=9,
                  markeredgecolor="black", markeredgewidth=0.6)
    pR, = ax.plot([xR[0]], [yR[0]], "o", color=colR, markersize=9,
                  markeredgecolor="black", markeredgewidth=0.6)

    tempo_text = ax.text(
        0.02, 0.98, "", transform=ax.transAxes, va="top", ha="left",
        fontsize=11, bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.7", alpha=0.9)
    )

    if TN < TR:
        ordine = "Arrivo: Newtoniana prima"
    elif TR < TN:
        ordine = "Arrivo: Relativistica prima"
    else:
        ordine = "Arrivo: pari (entro tolleranza numerica)"

    ax.text(0.02, 0.02, ordine, transform=ax.transAxes, va="bottom", ha="left",
            fontsize=10, color="0.15")

    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(-1.08, 1.08)
    ax.set_ylim(-1.08, 1.08)
    ax.set_xlabel("x/R")
    ax.set_ylabel("y/R")
    ax.set_title(f"Brachistocrona Newton vs Schwarzschild interna, Delta={delta_deg:.1f} deg, u={u:.4f}")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right")

    def update(frame_idx):
        tcur = t_frames[frame_idx]
        tempo_text.set_text(f"t = {tcur:7.2f} s   ({tcur/60.0:.3f} min)")

        xn, yn = posizione_al_tempo(xN, yN, T_cum_N, tcur)
        xr, yr = posizione_al_tempo(xR, yR, T_cum_R, tcur)
        pN.set_data([xn], [yn])
        pR.set_data([xr], [yr])

        return [tempo_text, pN, pR]

    ani = FuncAnimation(
        fig, update, frames=n_frames, interval=1000.0 / fps, blit=True, repeat=False
    )

    try:
        if out_file.lower().endswith(".gif"):
            from matplotlib.animation import PillowWriter
            ani.save(out_file, writer=PillowWriter(fps=fps), dpi=150)
        else:
            from matplotlib.animation import FFMpegWriter
            ani.save(out_file, writer=FFMpegWriter(fps=fps), dpi=150)
    except Exception:
        # fallback robusto: se mp4 fallisce, salva gif
        fallback = out_file.rsplit(".", 1)[0] + ".gif"
        from matplotlib.animation import PillowWriter
        ani.save(fallback, writer=PillowWriter(fps=fps), dpi=150)
        out_file = fallback

    plt.close(fig)

    return {
        "out_file": out_file,
        "TN": TN,
        "TR": TR,
        "q_rel": risR["q"],
        "tau_newton": tauN,
        "delta_err_deg_rel": math.degrees(risR["err_delta"]),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--delta-deg", type=float, default=90.0,
                        help="angolo centrale tra i due punti di superficie")
    parser.add_argument("--u", type=float, default=0.2,
                        help="compattezza u=r_s/R, deve stare in (0, 8/9)")
    parser.add_argument("--earth-u", action="store_true",
                        help="usa la compattezza terrestre reale")
    parser.add_argument("--R", type=float, default=R_TERRA,
                        help="raggio fisico in metri (scala i tempi in secondi)")
    parser.add_argument("--fps", type=int, default=18,
                        help="frame per secondo")
    parser.add_argument("--anim-dur", type=float, default=10.0,
                        help="durata video in secondi")
    parser.add_argument("--out-file", type=str,
                        default="confronto_brachistocrone_newton_rel.mp4",
                        help="file output (.mp4 o .gif)")
    parser.add_argument("--N", type=int, default=6000,
                        help="numero punti campionati per curva Newtoniana")
    args = parser.parse_args()

    if not (0.0 < args.delta_deg <= 180.0):
        raise ValueError("delta-deg deve stare in (0, 180].")
    if args.fps <= 0:
        raise ValueError("fps deve essere > 0.")
    if args.anim_dur <= 0.0:
        raise ValueError("anim-dur deve essere > 0.")
    if args.R <= 0.0:
        raise ValueError("R deve essere > 0.")
    if args.N < 1000:
        raise ValueError("N deve essere >= 1000 per buona stabilita' numerica.")

    if args.earth_u:
        u = compattezza_u(M_TERRA, R_TERRA)
    else:
        u = args.u

    if not (0.0 < u < U_MAX):
        raise ValueError("u deve stare in (0, 8/9).")

    out = salva_animazione_confronto(
        delta_deg=args.delta_deg,
        u=u,
        R=args.R,
        out_file=args.out_file,
        fps=args.fps,
        anim_dur=args.anim_dur,
        N=args.N,
    )

    print()
    print("=== Confronto brachistocrone (Newton vs Rel) ===")
    print(f"Delta                 = {args.delta_deg:.6f} deg")
    print(f"u                     = {u:.10f}")
    print(f"Tempo Newtoniano      = {out['TN']:.6f} s = {out['TN']/60.0:.6f} min")
    print(f"Tempo Relativistico   = {out['TR']:.6f} s = {out['TR']/60.0:.6f} min")
    print(f"Rapporto TR/TN        = {out['TR']/out['TN']:.9f}")
    print(f"q relativistico       = {out['q_rel']:.9f}")
    print(f"tau Newton (omega*T)  = {out['tau_newton']:.9f}")
    print(f"Err Delta rel         = {out['delta_err_deg_rel']:.3e} deg")
    print(f"Animazione salvata in: {out['out_file']}")


if __name__ == "__main__":
    main()
