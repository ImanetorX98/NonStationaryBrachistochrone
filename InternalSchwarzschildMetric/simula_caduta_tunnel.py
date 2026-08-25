#!/usr/bin/env python3
"""
simula_caduta_tunnel.py

Scopo:
    Simulare la caduta senza attrito di un grave lungo diverse traiettorie
    interne a una sfera newtoniana a densità uniforme, e confrontare i tempi
    con la traiettoria brachistocrona.

Fisica usata:
    Dentro una sfera a densità uniforme il campo gravitazionale è

        g(r) = omega^2 r

    dove

        omega = sqrt(GM/R^3)

    Una particella che parte da ferma dalla superficie ha energia conservata:

        v(r) = omega R sqrt(1 - (r/R)^2)

    In coordinate adimensionali rho = r/R:

        omega * T = integral ds_adim / sqrt(1 - rho^2)

    Quindi il tempo lungo una traiettoria si ottiene sommando numericamente

        dt = ds / v(r)

Traiettorie confrontate:
    1. tunnel rettilineo;
    2. tunnel parabolico poco profondo;
    3. tunnel parabolico più profondo;
    4. tunnel sinusoidale/coseno;
    5. tunnel brachistocrono teorico.

Esecuzione:
    python simula_caduta_tunnel.py
    python simula_caduta_tunnel.py --delta-deg 90
    python simula_caduta_tunnel.py --delta-deg 120 --plot

Dipendenze:
    numpy
    matplotlib solo se usi --plot
"""

import argparse
import math
import numpy as np


# ---------------------------------------------------------------------
# Costanti terrestri, usate solo per convertire il tempo in secondi/minuti
# ---------------------------------------------------------------------

G = 6.67430e-11
M_TERRA = 5.9722e24
R_TERRA = 6.371e6


def omega_uniforme(M=M_TERRA, R=R_TERRA):
    """
    Frequenza naturale del moto armonico dentro una sfera uniforme.
    """
    return math.sqrt(G * M / R**3)


# ---------------------------------------------------------------------
# Calcolo numerico del tempo lungo una curva campionata
# ---------------------------------------------------------------------

def tempo_lungo_curva(x, y, omega):
    """
    x, y sono array adimensionali: il raggio della sfera vale 1.
    Restituisce:
        tau = omega*T, tempo adimensionale
        T   = tempo fisico in secondi

    Metodo:
        - divido la traiettoria in piccoli segmenti;
        - calcolo ds per ogni segmento;
        - calcolo rho al punto medio del segmento;
        - uso v/(omega R) = sqrt(1-rho^2);
        - sommo dtau = ds / sqrt(1-rho^2).
    """
    x = np.asarray(x)
    y = np.asarray(y)

    dx = np.diff(x)
    dy = np.diff(y)
    ds = np.sqrt(dx*dx + dy*dy)

    xm = 0.5 * (x[:-1] + x[1:])
    ym = 0.5 * (y[:-1] + y[1:])
    rho = np.sqrt(xm*xm + ym*ym)

    # Dentro la sfera uniforme: v/(omega R)
    v_adim = np.sqrt(np.maximum(1.0 - rho*rho, 1e-14))

    dtau = ds / v_adim

    tau = np.sum(dtau)
    T = tau / omega
    return tau, T


def profilo_temporale_curva(x, y, omega):
    """
    Come `tempo_lungo_curva`, ma restituisce anche il profilo cumulativo
    del tempo lungo la curva.

    Restituisce:
        tau_tot, T_tot, tau_cum, T_cum
    """
    x = np.asarray(x)
    y = np.asarray(y)

    dx = np.diff(x)
    dy = np.diff(y)
    ds = np.sqrt(dx*dx + dy*dy)

    xm = 0.5 * (x[:-1] + x[1:])
    ym = 0.5 * (y[:-1] + y[1:])
    rho = np.sqrt(xm*xm + ym*ym)

    v_adim = np.sqrt(np.maximum(1.0 - rho*rho, 1e-14))
    dtau = ds / v_adim

    tau_cum = np.concatenate(([0.0], np.cumsum(dtau)))
    T_cum = tau_cum / omega
    tau_tot = tau_cum[-1]
    T_tot = T_cum[-1]
    return tau_tot, T_tot, tau_cum, T_cum


def posizione_al_tempo(x, y, T_cum, t):
    """
    Interpola la posizione della pallina lungo la traiettoria al tempo t.
    """
    if t <= 0.0:
        return float(x[0]), float(y[0])
    if t >= T_cum[-1]:
        return float(x[-1]), float(y[-1])

    i = np.searchsorted(T_cum, t, side="right") - 1
    i = int(np.clip(i, 0, len(T_cum) - 2))

    t0 = T_cum[i]
    t1 = T_cum[i + 1]
    alpha = (t - t0) / (t1 - t0) if t1 > t0 else 0.0

    xt = x[i] + alpha * (x[i + 1] - x[i])
    yt = y[i] + alpha * (y[i + 1] - y[i])
    return float(xt), float(yt)


def minuti(secondi):
    return secondi / 60.0


# ---------------------------------------------------------------------
# Costruzione delle traiettorie
# ---------------------------------------------------------------------

def estremi_superficie(delta):
    """
    Due punti sulla circonferenza unitaria separati da angolo centrale delta.
    Li mettiamo simmetrici rispetto all'asse x:

        A = (cos(delta/2), -sin(delta/2))
        B = (cos(delta/2), +sin(delta/2))
    """
    a = delta / 2.0
    A = np.array([math.cos(a), -math.sin(a)])
    B = np.array([math.cos(a), +math.sin(a)])
    return A, B


def curva_retta(delta, N=5000):
    """
    Tunnel rettilineo: segmento A-B.
    """
    A, B = estremi_superficie(delta)
    s = np.linspace(0.0, 1.0, N)
    x = (1.0 - s) * A[0] + s * B[0]
    y = (1.0 - s) * A[1] + s * B[1]
    return x, y


def curva_parabolica(delta, profondita=0.25, N=5000):
    """
    Famiglia semplice di tunnel curvi.

    Usiamo phi in [-delta/2, +delta/2] e imponiamo

        rho(phi) = 1 - profondita * (1 - (phi/a)^2)

    con a = delta/2.

    Quindi:
        rho(+-a) = 1      sulla superficie;
        rho(0)   = 1-h    punto più profondo.
    """
    a = delta / 2.0
    phi = np.linspace(-a, a, N)
    z = phi / a
    rho = 1.0 - profondita * (1.0 - z*z)

    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return x, y


def curva_coseno(delta, profondita=0.35, N=5000):
    """
    Altra famiglia semplice e liscia:

        rho(phi) = 1 - h/2 * (1 + cos(pi phi/a))

    con:
        rho(+-a)=1
        rho(0)=1-h
    """
    a = delta / 2.0
    phi = np.linspace(-a, a, N)
    rho = 1.0 - 0.5 * profondita * (1.0 + np.cos(math.pi * phi / a))

    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return x, y


def curva_brachistocrona(delta, N=5000):
    """
    Brachistocrona teorica nella sfera uniforme.

    Formula ricavata dalla legge di Snell radiale:

        q = r_*/R = 1 - delta/pi

    Per metà tunnel, dal punto più profondo alla superficie:

        rho(u) = sqrt(q^2 + (1-q^2) sin^2 u)

        phi(u) = atan2(sin u, q cos u) - q u

    con u in [0, pi/2].

    L'intera curva si ottiene per simmetria.
    """
    q = 1.0 - delta / math.pi

    # Caso limite antipodale: la brachistocrona è il diametro.
    if abs(q) < 1e-12:
        return curva_retta(delta, N=N)

    Nh = N // 2
    u = np.linspace(0.0, math.pi/2, Nh)

    rho = np.sqrt(q*q + (1.0 - q*q) * np.sin(u)**2)
    phi = np.arctan2(np.sin(u), q * np.cos(u)) - q * u

    # Ramo sinistro: dalla superficie al fondo.
    rho_left = rho[::-1]
    phi_left = -phi[::-1]

    # Ramo destro: dal fondo alla superficie.
    rho_right = rho[1:]
    phi_right = phi[1:]

    rho_full = np.concatenate([rho_left, rho_right])
    phi_full = np.concatenate([phi_left, phi_right])

    x = rho_full * np.cos(phi_full)
    y = rho_full * np.sin(phi_full)
    return x, y


def rho_min_brachistocrona(delta):
    """
    Minimo r/R raggiunto dalla brachistocrona teorica:
        q = 1 - delta/pi
    """
    return max(0.0, 1.0 - delta / math.pi)


def curva_parabolica_sotto_brach(delta, fattore_q=0.75, N=5000):
    """
    Curva parabolica con minimo imposto sotto quello della brachistocrona.

    Se q è il minimo brachistocrono, imponiamo:
        rho_min = fattore_q * q, con 0 < fattore_q < 1.
    """
    q = rho_min_brachistocrona(delta)
    rho_min = np.clip(fattore_q * q, 0.0, 0.999999)
    profondita = np.clip(1.0 - rho_min, 0.0, 0.999999)
    return curva_parabolica(delta, profondita=profondita, N=N)


def curva_coseno_sotto_brach(delta, fattore_q=0.45, N=5000):
    """
    Curva coseno con minimo imposto sotto quello della brachistocrona.
    """
    q = rho_min_brachistocrona(delta)
    rho_min = np.clip(fattore_q * q, 0.0, 0.999999)
    profondita = np.clip(1.0 - rho_min, 0.0, 0.999999)
    return curva_coseno(delta, profondita=profondita, N=N)


def tempo_brachistocrona_teorico(delta, omega):
    """
    Formula teorica del tempo minimo nella sfera uniforme:

        T_min = pi/omega * sqrt(1 - q^2)

    con q = 1 - delta/pi.
    """
    q = 1.0 - delta / math.pi
    tau = math.pi * math.sqrt(max(0.0, 1.0 - q*q))
    return tau, tau / omega


def salva_animazione(
    risultati, delta_deg, nome_file="animazione_tunnel.gif", fps=50, durata_video_s=20.0
):
    """
    Salva un'animazione con più palline che percorrono simultaneamente
    le traiettorie confrontate.
    """
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation

    T_max = max(item["T"] for item in risultati)
    T_anim = 1.15 * T_max
    n_frames = max(2, int(fps * durata_video_s))
    t_frames = np.linspace(0.0, T_anim, n_frames)

    colori = {
        "retta": "#4e79a7",
        "parabolica h=0.20": "#59a14f",
        "parabolica h=0.45": "#f28e2b",
        "coseno h=0.35": "#e15759",
        "parab sotto q x0.75": "#76b7b2",
        "coseno sotto q x0.45": "#af7aa1",
        "brachistocrona": "#c9a227",
    }

    fig, ax = plt.subplots(figsize=(8, 8))

    # Superficie della sfera
    tt = np.linspace(0.0, 2.0 * math.pi, 1000)
    ax.plot(np.cos(tt), np.sin(tt), "k-", linewidth=1.0)

    punti = {}
    ordine_arrivo = sorted(risultati, key=lambda item: item["T"])
    ordine_text = "Ordine arrivo: " + "  <  ".join(
        f"{idx+1}) {item['nome']}" for idx, item in enumerate(ordine_arrivo)
    )

    for item in risultati:
        nome = item["nome"]
        x = item["x"]
        y = item["y"]
        c = colori.get(nome, None)
        lw = 3.0 if nome == "brachistocrona" else 1.7
        alpha_linea = 1.0 if nome == "brachistocrona" else 0.95
        ax.plot(x, y, color=c, linewidth=lw, alpha=alpha_linea)

        ms = 9 if nome == "brachistocrona" else 7
        punti[nome], = ax.plot(
            [x[0]], [y[0]], "o", color=c, markersize=ms,
            markeredgecolor="black", markeredgewidth=0.6
        )

    tempo_text = ax.text(
        0.02, 0.98, "", transform=ax.transAxes, va="top", ha="left",
        fontsize=11, bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.7", alpha=0.9)
    )
    ordine_arrivo_text = ax.text(
        0.02, 0.02, ordine_text, transform=ax.transAxes, va="bottom", ha="left",
        fontsize=9, color="0.15"
    )
    _ = ordine_arrivo_text  # evita warning su variabile non usata in alcuni ambienti

    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(-1.08, 1.08)
    ax.set_ylim(-1.08, 1.08)
    ax.set_xlabel("x/R")
    ax.set_ylabel("y/R")
    ax.set_title(f"Caduta in tunnel: Delta = {delta_deg:.1f}°")
    ax.grid(True, alpha=0.3)

    def update(frame_idx):
        t = t_frames[frame_idx]
        tempo_text.set_text(f"t = {t:6.1f} s   ({minuti(t):.2f} min)")

        for item in risultati:
            nome = item["nome"]
            x_t, y_t = posizione_al_tempo(item["x"], item["y"], item["T_cum"], t)
            punti[nome].set_data([x_t], [y_t])

        artists = [tempo_text]
        artists.extend(punti.values())
        return artists

    ani = FuncAnimation(
        fig, update, frames=n_frames, interval=1000.0 / fps, blit=True, repeat=False
    )

    try:
        from matplotlib.animation import PillowWriter
        ani.save(nome_file, writer=PillowWriter(fps=fps), dpi=150)
    except Exception:
        # Fallback MP4 se Pillow non è disponibile.
        from matplotlib.animation import FFMpegWriter
        if nome_file.lower().endswith(".gif"):
            nome_file = nome_file[:-4] + ".mp4"
        ani.save(nome_file, writer=FFMpegWriter(fps=fps), dpi=150)

    plt.close(fig)
    return nome_file


# ---------------------------------------------------------------------
# Simulazione tabellare
# ---------------------------------------------------------------------

def confronta(
    delta_deg=90.0,
    plot=False,
    animate=False,
    anim_file="animazione_tunnel.gif",
    fps=50,
    anim_dur=20.0,
):
    delta = math.radians(delta_deg)
    omega = omega_uniforme()

    print()
    print("=== Simulazione caduta in tunnel interni ===")
    print(f"Angolo centrale Delta = {delta_deg:.2f} gradi")
    print(f"omega = sqrt(GM/R^3) = {omega:.6e} 1/s")
    q_brach = rho_min_brachistocrona(delta)
    print(f"Minimo brachistocrona q = r_min/R = {q_brach:.6f}")
    print()

    traiettorie = [
        ("retta", curva_retta(delta)),
        ("parabolica h=0.20", curva_parabolica(delta, profondita=0.20)),
        ("parabolica h=0.45", curva_parabolica(delta, profondita=0.45)),
        ("coseno h=0.35", curva_coseno(delta, profondita=0.35)),
        ("brachistocrona", curva_brachistocrona(delta)),
    ]
    if q_brach > 1e-6:
        traiettorie.extend([
            ("parab sotto q x0.75", curva_parabolica_sotto_brach(delta, fattore_q=0.75)),
            ("coseno sotto q x0.45", curva_coseno_sotto_brach(delta, fattore_q=0.45)),
        ])
    else:
        print("Nota: q=0 (caso antipodale), non esistono traiettorie interne sotto il minimo brachistocrono.")
        print()

    risultati = []

    for nome, (x, y) in traiettorie:
        tau, T, tau_cum, T_cum = profilo_temporale_curva(x, y, omega)
        risultati.append({
            "nome": nome,
            "tau": tau,
            "T": T,
            "x": x,
            "y": y,
            "tau_cum": tau_cum,
            "T_cum": T_cum,
        })

    # Aggiungiamo anche il tempo teorico della brachistocrona.
    tau_th, T_th = tempo_brachistocrona_teorico(delta, omega)

    risultati_ordinati = sorted(risultati, key=lambda item: item["tau"])

    print("Tempi numerici:")
    print(f"{'traiettoria':22s} {'omega*T':>12s} {'T [s]':>12s} {'T [min]':>12s}")
    print("-" * 62)

    for item in risultati_ordinati:
        print(f"{item['nome']:22s} {item['tau']:12.6f} {item['T']:12.3f} {minuti(item['T']):12.3f}")

    print()
    print("Brachistocrona teorica:")
    print(f"omega*T_min teorico = {tau_th:.6f}")
    print(f"T_min teorico       = {T_th:.3f} s = {minuti(T_th):.3f} min")

    print()
    print("Rapporti rispetto alla brachistocrona numerica:")
    tau_b = [r["tau"] for r in risultati if r["nome"] == "brachistocrona"][0]
    for item in risultati:
        print(f"{item['nome']:22s} T/T_brach = {item['tau']/tau_b:.6f}")

    if plot:
        import matplotlib.pyplot as plt

        plt.figure(figsize=(8, 8))

        # superficie
        t = np.linspace(0, 2*math.pi, 1000)
        plt.plot(np.cos(t), np.sin(t), "k-", linewidth=1, label="superficie")

        for item in risultati:
            if item["nome"] == "brachistocrona":
                plt.plot(item["x"], item["y"], linewidth=3, label=f"{item['nome']}: {minuti(item['T']):.2f} min")
            else:
                plt.plot(item["x"], item["y"], linewidth=1.5, label=f"{item['nome']}: {minuti(item['T']):.2f} min")

        plt.gca().set_aspect("equal", adjustable="box")
        plt.xlabel("x/R")
        plt.ylabel("y/R")
        plt.title(f"Confronto tunnel, Delta = {delta_deg:.1f}°")
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.tight_layout()
        plt.savefig("confronto_tunnel.png", dpi=180)

        print()
        print("Grafico salvato in: confronto_tunnel.png")

    if animate:
        nome_salvato = salva_animazione(
            risultati=risultati,
            delta_deg=delta_deg,
            nome_file=anim_file,
            fps=fps,
            durata_video_s=anim_dur,
        )
        print()
        print(f"Animazione salvata in: {nome_salvato} (durata video ~{anim_dur:.1f} s)")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--delta-deg", type=float, default=90.0,
                        help="angolo centrale tra i due punti sulla superficie")
    parser.add_argument("--plot", action="store_true",
                        help="salva un grafico delle traiettorie")
    parser.add_argument("--animate", action="store_true",
                        help="salva un'animazione con palline in moto sulle traiettorie")
    parser.add_argument("--anim-file", type=str, default="animazione_tunnel.gif",
                        help="nome file per l'animazione (.gif o .mp4)")
    parser.add_argument("--fps", type=int, default=50,
                        help="frame per secondo dell'animazione")
    parser.add_argument("--anim-dur", type=float, default=20.0,
                        help="durata del video animato in secondi")
    args = parser.parse_args()

    if not (0.0 < args.delta_deg <= 180.0):
        raise ValueError("delta-deg deve stare tra 0 e 180 gradi.")

    if args.fps <= 0:
        raise ValueError("fps deve essere > 0.")
    if args.anim_dur <= 0.0:
        raise ValueError("anim-dur deve essere > 0.")

    confronta(
        delta_deg=args.delta_deg,
        plot=args.plot,
        animate=args.animate,
        anim_file=args.anim_file,
        fps=args.fps,
        anim_dur=args.anim_dur,
    )


if __name__ == "__main__":
    main()
