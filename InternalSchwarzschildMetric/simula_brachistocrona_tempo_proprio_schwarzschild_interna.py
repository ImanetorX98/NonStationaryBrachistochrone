#!/usr/bin/env python3
"""
simula_brachistocrona_tempo_proprio_schwarzschild_interna.py

Scopo:
    Visualizzare la brachistocrona che minimizza il tempo proprio all'interno
    di una sfera uniforme descritta dalla metrica interna di Schwarzschild.

Nota fisica importante:
    Una minimizzazione "assoluta" del tempo proprio fra due eventi non e'
    ben posta senza vincoli (si puo' tendere a traiettorie quasi nulle con
    tempo proprio arbitrariamente piccolo). Qui usiamo il vincolo standard:
    particella libera rilasciata da ferma alla superficie (energia specifica
    fissata e = A(1)).

Metrica equatoriale (x = r/R):
    ds^2 = -A(x)^2 c^2 dt^2 + B(x)^2 R^2 dx^2 + R^2 x^2 dphi^2

con:
    A(x) = 1/2 * (3*sqrt(1-u) - sqrt(1-u*x^2))
    B(x) = 1/sqrt(1-u*x^2)
    u    = r_s/R = 2GM/(c^2 R), con u < 8/9.

Con e = A(1), la velocita' locale soddisfa:
    v(x)/c = sqrt(1 - A(x)^2/e^2)

Minimizzazione:
    - tempo coordinato: n_t(x)  = 1 / (A * sqrt(1-A^2/e^2))
    - tempo proprio:    n_tau(x)= A / sqrt(e^2-A^2)

Entrambi i problemi hanno prima integrale:
    dphi/dx = [B(x)/x] / sqrt(W(x)/W(q) - 1)
con:
    W(x) = x^2 * n(x)^2
    q    = r_min/R.

Uso:
    python simula_brachistocrona_tempo_proprio_schwarzschild_interna.py --delta-deg 90 --plot
    python simula_brachistocrona_tempo_proprio_schwarzschild_interna.py --delta-deg 120 --u 0.2 --plot --compare-coordinate
"""

import argparse
import math
import numpy as np


G = 6.67430e-11
C = 299792458.0
M_TERRA = 5.9722e24
R_TERRA = 6.371e6
U_MAX = 8.0 / 9.0


def compattezza_u(M, R):
    return 2.0 * G * M / (C * C * R)


def A_lapse(x, u):
    return 0.5 * (3.0 * math.sqrt(1.0 - u) - np.sqrt(1.0 - u * x * x))


def B_radiale(x, u):
    return 1.0 / np.sqrt(1.0 - u * x * x)


def energia_superficie(u):
    return math.sqrt(1.0 - u)


def W_coord(x, u):
    """
    W per minimizzazione del tempo coordinato.
    """
    e = energia_superficie(u)
    Ax = A_lapse(x, u)
    den = Ax * Ax * (1.0 - (Ax * Ax) / (e * e))
    den = np.maximum(den, 1e-20)
    return x * x / den


def W_proprio(x, u):
    """
    W per minimizzazione del tempo proprio (con e fissata).
    """
    e = energia_superficie(u)
    Ax = A_lapse(x, u)
    den = np.maximum(e * e - Ax * Ax, 1e-20)
    return x * x * Ax * Ax / den


def W_jacobi(x, u, kind):
    if kind == "coord":
        return W_coord(x, u)
    if kind == "proprio":
        return W_proprio(x, u)
    raise ValueError("kind deve essere 'coord' o 'proprio'.")


def dphi_dx(x, q, u, kind):
    Wq = W_jacobi(q, u, kind)
    ratio = W_jacobi(x, u, kind) / Wq - 1.0
    ratio = np.maximum(ratio, 1e-18)
    return (B_radiale(x, u) / x) / np.sqrt(ratio)


def phi_half_midpoint(q, u, kind, x_top=1.0, n=4000):
    """
    Integrale phi(q->x_top) con sostituzione x = q + (x_top-q) s^2, s in [0,1].
    """
    if not (0.0 < q < x_top <= 1.0):
        raise ValueError("Richiesto 0 < q < x_top <= 1.")

    i = np.arange(n)
    s = (i + 0.5) / n
    x = q + (x_top - q) * s * s
    dx_ds = 2.0 * (x_top - q) * s
    integrando = dphi_dx(x, q, u, kind) * dx_ds
    return integrando.sum() / n


def delta_da_q(q, u, kind, n=4000):
    return 2.0 * phi_half_midpoint(q=q, u=u, kind=kind, n=n)


def tempi_adim_da_q(q, u, kind, n=6000):
    """
    Restituisce:
      tbar  = c*T/R
      taubar= c*Tau/R
    lungo la curva associata al dato q e al dato funzionale (kind).
    """
    i = np.arange(n)
    s = (i + 0.5) / n
    x = q + (1.0 - q) * s * s
    dx_ds = 2.0 * (1.0 - q) * s

    phip = dphi_dx(x, q, u, kind)
    ds_dx = np.sqrt(B_radiale(x, u) ** 2 + (x * phip) ** 2)

    e = energia_superficie(u)
    Ax = A_lapse(x, u)
    vel_loc = np.sqrt(np.maximum(1.0 - (Ax * Ax) / (e * e), 1e-18))

    dtbar_dx = ds_dx / np.maximum(Ax * vel_loc, 1e-18)
    dtaubar_dx = ds_dx * Ax / np.sqrt(np.maximum(e * e - Ax * Ax, 1e-18))

    dtbar_ds = dtbar_dx * dx_ds
    dtaubar_ds = dtaubar_dx * dx_ds

    tbar_half = dtbar_ds.sum() / n
    taubar_half = dtaubar_ds.sum() / n
    return 2.0 * tbar_half, 2.0 * taubar_half


def trova_q_candidati(delta_target, u, kind, n_scan=260, n_int=2200):
    q_min = 1e-4
    q_max = 1.0 - 1e-6
    qs = np.linspace(q_min, q_max, n_scan)
    vals = np.array([delta_da_q(q, u, kind, n=n_int) - delta_target for q in qs])

    candidati = []
    for i in range(len(qs) - 1):
        f1 = vals[i]
        f2 = vals[i + 1]
        if np.isnan(f1) or np.isnan(f2):
            continue
        if f1 == 0.0:
            candidati.append(qs[i])
            continue
        if f1 * f2 < 0.0:
            a = qs[i]
            b = qs[i + 1]
            fa = f1
            fb = f2
            for _ in range(60):
                m = 0.5 * (a + b)
                fm = delta_da_q(m, u, kind, n=n_int + 600) - delta_target
                if fa * fm <= 0.0:
                    b = m
                    fb = fm
                else:
                    a = m
                    fa = fm
            _ = fb
            candidati.append(0.5 * (a + b))

    candidati_ordinati = []
    for q in sorted(candidati):
        if (not candidati_ordinati) or abs(q - candidati_ordinati[-1]) > 1e-5:
            candidati_ordinati.append(q)
    return candidati_ordinati


def scegli_q_ottimo(delta_target, u, kind):
    """
    Sceglie q con delta(q)=delta_target e minimo del funzionale indicato da kind.
    """
    candidati = trova_q_candidati(delta_target, u, kind)
    if not candidati:
        raise RuntimeError("Nessuna soluzione q trovata: prova a cambiare Delta o u.")

    best = None
    for q in candidati:
        tbar, taubar = tempi_adim_da_q(q, u, kind=kind, n=5000)
        costo = taubar if kind == "proprio" else tbar
        if (best is None) or (costo < best["costo"]):
            best = {"q": q, "tbar": tbar, "taubar": taubar, "costo": costo}

    q_best = best["q"]
    delta_num = delta_da_q(q_best, u, kind, n=7000)
    err = delta_num - delta_target
    return q_best, best["tbar"], best["taubar"], delta_num, err, candidati


def profilo_phi_mezza_curva(q, u, kind, n=2200):
    s_edges = np.linspace(0.0, 1.0, n + 1)
    i = np.arange(n)
    s_mid = (i + 0.5) / n

    x_edges = q + (1.0 - q) * s_edges * s_edges
    x_mid = q + (1.0 - q) * s_mid * s_mid
    dx_ds_mid = 2.0 * (1.0 - q) * s_mid
    dphi_ds_mid = dphi_dx(x_mid, q, u, kind) * dx_ds_mid

    ds = 1.0 / n
    dphi = dphi_ds_mid * ds
    phi_edges = np.concatenate(([0.0], np.cumsum(dphi)))
    return x_edges, phi_edges


def costruisci_traiettoria(delta_target, u, kind, n=2200):
    q, tbar, taubar, delta_num, err, candidati = scegli_q_ottimo(delta_target, u, kind)

    x_half, phi_half = profilo_phi_mezza_curva(q, u, kind=kind, n=n)
    scala = (0.5 * delta_target) / phi_half[-1]
    phi_half = phi_half * scala

    x_left = x_half[::-1]
    phi_left = -phi_half[::-1]
    x_right = x_half[1:]
    phi_right = phi_half[1:]

    x_rad = np.concatenate([x_left, x_right])
    phi = np.concatenate([phi_left, phi_right])

    x_cart = x_rad * np.cos(phi)
    y_cart = x_rad * np.sin(phi)

    return {
        "kind": kind,
        "q": q,
        "tbar": tbar,
        "taubar": taubar,
        "delta_num": delta_num,
        "err_delta": err,
        "candidati_q": candidati,
        "x_rad": x_rad,
        "phi": phi,
        "x": x_cart,
        "y": y_cart,
    }


def visualizza_traiettoria(
    ris_tau,
    delta_deg,
    u,
    out_file,
    ris_coord=None,
):
    import matplotlib.pyplot as plt

    t = np.linspace(0.0, 2.0 * math.pi, 1000)
    plt.figure(figsize=(8, 8))
    plt.plot(np.cos(t), np.sin(t), "k-", linewidth=1.2, label="superficie x=1")

    plt.plot(
        ris_tau["x"],
        ris_tau["y"],
        color="#c9a227",
        linewidth=2.8,
        label="ottima per tempo proprio",
    )

    if ris_coord is not None:
        plt.plot(
            ris_coord["x"],
            ris_coord["y"],
            color="#4e79a7",
            linewidth=2.0,
            linestyle="--",
            label="ottima per tempo coordinato",
        )

    xa, ya = ris_tau["x"][0], ris_tau["y"][0]
    xb, yb = ris_tau["x"][-1], ris_tau["y"][-1]
    s = np.linspace(0.0, 1.0, 500)
    xr = (1.0 - s) * xa + s * xb
    yr = (1.0 - s) * ya + s * yb
    plt.plot(xr, yr, ":", color="0.45", linewidth=1.4, alpha=0.9, label="corda riferimento")

    plt.gca().set_aspect("equal", adjustable="box")
    plt.xlim(-1.08, 1.08)
    plt.ylim(-1.08, 1.08)
    plt.xlabel("x = r cos(phi) / R")
    plt.ylabel("y = r sin(phi) / R")
    plt.title(f"Brachistocrona per tempo proprio, Delta={delta_deg:.1f} deg, u={u:.4f}")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_file, dpi=180)
    plt.close()
    return out_file


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--delta-deg", type=float, default=90.0,
                        help="angolo centrale tra i punti estremi in gradi (0, 180]")
    parser.add_argument("--u", type=float, default=0.2,
                        help="compattezza u = r_s/R = 2GM/(c^2 R), in (0, 8/9)")
    parser.add_argument("--R", type=float, default=R_TERRA,
                        help="raggio fisico R in metri per conversione tempi")
    parser.add_argument("--earth-u", action="store_true",
                        help="ignora --u e usa la compattezza terrestre")
    parser.add_argument("--plot", action="store_true",
                        help="salva il grafico della traiettoria ottima")
    parser.add_argument("--out", type=str,
                        default="brachistocrona_tempo_proprio_schwarzschild_interna.png",
                        help="nome file grafico PNG")
    parser.add_argument("--compare-coordinate", action="store_true",
                        help="disegna anche la soluzione ottima per tempo coordinato")
    args = parser.parse_args()

    if not (0.0 < args.delta_deg <= 180.0):
        raise ValueError("delta-deg deve stare in (0, 180].")
    if args.R <= 0.0:
        raise ValueError("R deve essere > 0.")

    if args.earth_u:
        u = compattezza_u(M_TERRA, R_TERRA)
    else:
        u = args.u

    if not (0.0 < u < U_MAX):
        raise ValueError("u deve stare in (0, 8/9).")

    delta = math.radians(args.delta_deg)

    ris_tau = costruisci_traiettoria(delta_target=delta, u=u, kind="proprio")
    T_tau_coord = ris_tau["tbar"] * args.R / C
    Tau_prop = ris_tau["taubar"] * args.R / C

    print()
    print("=== Brachistocrona (minimo tempo proprio) ===")
    print(f"Delta target         = {args.delta_deg:.6f} deg")
    print(f"u = r_s/R            = {u:.10f}")
    print(f"q ottimo proprio     = {ris_tau['q']:.10f}")
    print(f"Delta numerico       = {math.degrees(ris_tau['delta_num']):.6f} deg")
    print(f"Errore Delta         = {math.degrees(ris_tau['err_delta']):.3e} deg")
    print(f"q candidati          = {len(ris_tau['candidati_q'])}")
    print(f"c*T/R su curva       = {ris_tau['tbar']:.10f}")
    print(f"c*Tau/R minimo       = {ris_tau['taubar']:.10f}")
    print(f"T coordinato         = {T_tau_coord:.6f} s = {T_tau_coord/60.0:.6f} min")
    print(f"Tau proprio          = {Tau_prop:.6f} s = {Tau_prop/60.0:.6f} min")

    ris_coord = None
    if args.compare_coordinate:
        ris_coord = costruisci_traiettoria(delta_target=delta, u=u, kind="coord")
        T_coord = ris_coord["tbar"] * args.R / C
        Tau_coord_path = ris_coord["taubar"] * args.R / C

        print()
        print("--- Confronto con soluzione ottima per tempo coordinato ---")
        print(f"q ottimo coordinato  = {ris_coord['q']:.10f}")
        print(f"T coordinato minimo  = {T_coord:.6f} s = {T_coord/60.0:.6f} min")
        print(f"Tau su curva coord   = {Tau_coord_path:.6f} s = {Tau_coord_path/60.0:.6f} min")
        print(f"Rapporto Tau_proprio/Tau_coord_path = {Tau_prop/max(Tau_coord_path,1e-18):.9f}")
        print(f"Rapporto T_proprio_path/T_coord_min = {T_tau_coord/max(T_coord,1e-18):.9f}")

    if args.plot:
        path = visualizza_traiettoria(
            ris_tau=ris_tau,
            delta_deg=args.delta_deg,
            u=u,
            out_file=args.out,
            ris_coord=ris_coord,
        )
        print()
        print(f"Grafico salvato in: {path}")


if __name__ == "__main__":
    main()
