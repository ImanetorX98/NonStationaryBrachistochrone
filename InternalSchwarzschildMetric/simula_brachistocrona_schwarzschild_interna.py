#!/usr/bin/env python3
"""
simula_brachistocrona_schwarzschild_interna.py

Scopo:
    Calcolare numericamente la traiettoria di caduta a tempo coordinato minimo
    (brachistocrona relativistica) dentro una sfera uniforme descritta dalla
    metrica di Schwarzschild interna, e visualizzarla.

Idea fisica (equatoriale, metrica statica sferica):
    ds^2 = -A(r)^2 c^2 dt^2 + B(r)^2 dr^2 + r^2 dphi^2

    con metrica interna (densita' costante):
      A(x) = 1/2 * (3*sqrt(1-u) - sqrt(1-u*x^2))
      B(x) = 1/sqrt(1-u*x^2)

    dove:
      x = r/R in [0,1]
      u = r_s / R = 2GM/(c^2 R), con u < 8/9 (limite di Buchdahl).

    Per una particella rilasciata da ferma alla superficie x=1:
      e = A(1) = sqrt(1-u)

    La minimizzazione del tempo coordinato t equivale a una geodesica di una
    metrica di Jacobi (conforme), da cui segue una prima integrale:

      phi(x; q) = integral_q^x [ B(xi)/xi ] / sqrt(W(xi)/W(q)-1) dxi

      W(x) = x^2 / ( A(x)^2 * (1 - A(x)^2/e^2) )

    La separazione angolare totale e':
      Delta(q) = 2 * phi(1; q)

    Dato Delta, risolviamo numericamente Delta(q)=Delta_target.

Uso:
    python simula_brachistocrona_schwarzschild_interna.py --delta-deg 90
    python simula_brachistocrona_schwarzschild_interna.py --delta-deg 120 --u 0.2 --plot
"""

import argparse
import math
import numpy as np


# Costanti fisiche (solo per conversione tempi in secondi)
G = 6.67430e-11
C = 299792458.0
M_TERRA = 5.9722e24
R_TERRA = 6.371e6

# Limite compattezza per soluzione interna a densita' costante
U_MAX = 8.0 / 9.0


def compattezza_u(M, R):
    return 2.0 * G * M / (C * C * R)


def A_lapse(x, u):
    return 0.5 * (3.0 * math.sqrt(1.0 - u) - np.sqrt(1.0 - u * x * x))


def B_radiale(x, u):
    return 1.0 / np.sqrt(1.0 - u * x * x)


def W_jacobi(x, u):
    """
    Fattore che compare nella prima integrale per la traiettoria ottima.
    """
    e = math.sqrt(1.0 - u)  # A(1), particella rilasciata da ferma alla superficie
    Ax = A_lapse(x, u)
    den = Ax * Ax * (1.0 - (Ax * Ax) / (e * e))
    den = np.maximum(den, 1e-20)
    return x * x / den


def _phi_half_midpoint(q, x_top, u, n=4000):
    """
    Integrale phi(q->x_top) con regola del punto medio e sostituzione:
        x = q + (x_top-q) s^2, s in [0,1]
    che regolarizza la singolarita' integrabile in x=q.
    """
    if not (0.0 < q < x_top <= 1.0):
        raise ValueError("Richiesto 0 < q < x_top <= 1.")

    i = np.arange(n)
    s = (i + 0.5) / n
    x = q + (x_top - q) * s * s
    dx_ds = 2.0 * (x_top - q) * s

    Wq = W_jacobi(q, u)
    ratio = W_jacobi(x, u) / Wq - 1.0
    # ratio > 0 su (q, x_top]; uso clip difensivo per rumore numerico
    ratio = np.maximum(ratio, 1e-18)

    dphi_dx = (B_radiale(x, u) / x) / np.sqrt(ratio)
    integrando = dphi_dx * dx_ds
    return integrando.sum() / n


def delta_da_q(q, u, n=4000):
    return 2.0 * _phi_half_midpoint(q=q, x_top=1.0, u=u, n=n)


def dphi_dx(x, q, u):
    Wq = W_jacobi(q, u)
    ratio = W_jacobi(x, u) / Wq - 1.0
    ratio = np.maximum(ratio, 1e-18)
    return (B_radiale(x, u) / x) / np.sqrt(ratio)


def tempo_coord_adim_da_q(q, u, n=6000):
    """
    Tempo coordinato totale adimensionale:
        t_bar = c*T/R
    ottenuto integrando da q a 1 e raddoppiando.
    """
    i = np.arange(n)
    s = (i + 0.5) / n
    x = q + (1.0 - q) * s * s
    dx_ds = 2.0 * (1.0 - q) * s

    Ax = A_lapse(x, u)
    e = math.sqrt(1.0 - u)
    vel_loc = np.sqrt(np.maximum(1.0 - (Ax * Ax) / (e * e), 1e-18))

    phip = dphi_dx(x, q, u)
    ds_dx = np.sqrt(B_radiale(x, u) ** 2 + (x * phip) ** 2)
    dtbar_dx = ds_dx / (Ax * vel_loc)
    dtbar_ds = dtbar_dx * dx_ds

    tbar_half = dtbar_ds.sum() / n
    return 2.0 * tbar_half


def trova_q_candidati(delta_target, u, n_scan=260, n_int=2200):
    """
    Cerca tutti i candidati q in (0,1) tali che Delta(q)=Delta_target.
    """
    q_min = 1e-4
    q_max = 1.0 - 1e-6
    qs = np.linspace(q_min, q_max, n_scan)
    vals = np.array([delta_da_q(q, u, n=n_int) - delta_target for q in qs])

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
                fm = delta_da_q(m, u, n=n_int + 600) - delta_target
                if fa * fm <= 0.0:
                    b = m
                    fb = fm
                else:
                    a = m
                    fa = fm
            _ = fb  # solo per simmetria con fa; evita warning in alcuni linter
            candidati.append(0.5 * (a + b))

    # Rimuove duplicati numerici
    candidati_ordinati = []
    for q in sorted(candidati):
        if (not candidati_ordinati) or abs(q - candidati_ordinati[-1]) > 1e-5:
            candidati_ordinati.append(q)

    return candidati_ordinati


def scegli_q_ottimo(delta_target, u):
    """
    Se esistono piu' soluzioni q della stessa Delta, sceglie quella con
    tempo coordinato minore.
    """
    candidati = trova_q_candidati(delta_target, u)
    if not candidati:
        raise RuntimeError(
            "Nessuna soluzione q trovata: prova a cambiare Delta o compattezza u."
        )

    best = None
    for q in candidati:
        tbar = tempo_coord_adim_da_q(q, u, n=5000)
        if (best is None) or (tbar < best["tbar"]):
            best = {"q": q, "tbar": tbar}

    q_best = best["q"]
    delta_num = delta_da_q(q_best, u, n=7000)
    err = delta_num - delta_target
    return q_best, best["tbar"], delta_num, err, candidati


def profilo_phi_mezza_curva(q, u, n=2200):
    """
    Restituisce campioni x(s), phi(s) per mezza curva (da q a 1), con
    integrazione cumulativa al punto medio.
    """
    s_edges = np.linspace(0.0, 1.0, n + 1)
    i = np.arange(n)
    s_mid = (i + 0.5) / n

    x_edges = q + (1.0 - q) * s_edges * s_edges

    x_mid = q + (1.0 - q) * s_mid * s_mid
    dx_ds_mid = 2.0 * (1.0 - q) * s_mid
    dphi_ds_mid = dphi_dx(x_mid, q, u) * dx_ds_mid

    ds = 1.0 / n
    dphi = dphi_ds_mid * ds
    phi_edges = np.concatenate(([0.0], np.cumsum(dphi)))

    return x_edges, phi_edges


def costruisci_traiettoria(delta_target, u, n=2200):
    q, tbar, delta_num, err, candidati = scegli_q_ottimo(delta_target, u)

    x_half, phi_half = profilo_phi_mezza_curva(q, u, n=n)
    # Normalizza all'angolo richiesto (micro-correzione numerica)
    scala = (0.5 * delta_target) / phi_half[-1]
    phi_half = phi_half * scala

    # Ramo sinistro (superficie->fondo) e destro (fondo->superficie)
    x_left = x_half[::-1]
    phi_left = -phi_half[::-1]
    x_right = x_half[1:]
    phi_right = phi_half[1:]

    x_rad = np.concatenate([x_left, x_right])
    phi = np.concatenate([phi_left, phi_right])

    x_cart = x_rad * np.cos(phi)
    y_cart = x_rad * np.sin(phi)

    return {
        "q": q,
        "tbar": tbar,
        "delta_num": delta_num,
        "err_delta": err,
        "candidati_q": candidati,
        "x_rad": x_rad,
        "phi": phi,
        "x": x_cart,
        "y": y_cart,
    }


def visualizza_traiettoria(ris, delta_deg, u, out_file="brachistocrona_schwarzschild_interna.png"):
    import matplotlib.pyplot as plt

    x = ris["x"]
    y = ris["y"]

    t = np.linspace(0.0, 2.0 * math.pi, 1000)
    plt.figure(figsize=(8, 8))
    plt.plot(np.cos(t), np.sin(t), "k-", linewidth=1.2, label="superficie x=1")
    plt.plot(x, y, color="#c9a227", linewidth=2.8, label="traiettoria ottima")

    # Segmento rettilineo di riferimento sugli stessi estremi
    xa, ya = x[0], y[0]
    xb, yb = x[-1], y[-1]
    s = np.linspace(0.0, 1.0, 500)
    xr = (1.0 - s) * xa + s * xb
    yr = (1.0 - s) * ya + s * yb
    plt.plot(xr, yr, "--", color="#4e79a7", linewidth=1.4, alpha=0.9, label="corda (riferimento)")

    plt.gca().set_aspect("equal", adjustable="box")
    plt.xlim(-1.08, 1.08)
    plt.ylim(-1.08, 1.08)
    plt.xlabel("x = r cos(phi) / R")
    plt.ylabel("y = r sin(phi) / R")
    plt.title(f"Brachistocrona (Schwarzschild interna), Delta={delta_deg:.1f} deg, u={u:.4f}")
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
                        help="compattezza u = r_s/R = 2GM/(c^2 R), deve stare in (0, 8/9)")
    parser.add_argument("--R", type=float, default=R_TERRA,
                        help="raggio fisico R in metri, usato per convertire il tempo coordinato in secondi")
    parser.add_argument("--plot", action="store_true",
                        help="salva il grafico della traiettoria ottima")
    parser.add_argument("--out", type=str, default="brachistocrona_schwarzschild_interna.png",
                        help="nome file del grafico PNG")
    parser.add_argument("--earth-u", action="store_true",
                        help="ignora --u e usa la compattezza terrestre")
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
    ris = costruisci_traiettoria(delta_target=delta, u=u)

    T_coord = ris["tbar"] * args.R / C

    print()
    print("=== Brachistocrona in Schwarzschild interna ===")
    print(f"Delta target      = {args.delta_deg:.6f} deg")
    print(f"u = r_s/R         = {u:.10f}")
    print(f"q ottimo (r_min/R)= {ris['q']:.10f}")
    print(f"Delta numerico    = {math.degrees(ris['delta_num']):.6f} deg")
    print(f"Errore Delta      = {math.degrees(ris['err_delta']):.3e} deg")
    print(f"q candidati       = {len(ris['candidati_q'])}")
    print(f"c*T/R (ottimo)    = {ris['tbar']:.10f}")
    print(f"T coordinato      = {T_coord:.6f} s = {T_coord/60.0:.6f} min")

    if args.plot:
        path = visualizza_traiettoria(
            ris=ris,
            delta_deg=args.delta_deg,
            u=u,
            out_file=args.out,
        )
        print()
        print(f"Grafico salvato in: {path}")


if __name__ == "__main__":
    main()
