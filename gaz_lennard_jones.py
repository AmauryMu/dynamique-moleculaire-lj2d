"""
Dynamique moléculaire d'un gaz de Lennard-Jones en 2D.

N particules en interaction de paire (potentiel de Lennard-Jones) dans une
boîte à murs réfléchissants. Intégration de Verlet-vitesse, thermostat de
Berendsen, calcul de la température, de la pression (théorème du viriel)
et comparaison des vitesses à la distribution de Maxwell-Boltzmann 2D.
Un curseur permet de changer la température cible pendant la simulation.

Unités réduites de Lennard-Jones : m = 1, kB = 1, epsilon = 1, sigma = 1.

Projet M1 Physique, CY Cergy Paris Université (février 2026), réalisé en binôme.

Utilisation :
    python gaz_lennard_jones.py            # simulation interactive
    python gaz_lennard_jones.py --sauver   # enregistre un GIF dans figures/
                                           # (condensation, fusion, gaz, recondensation)
"""

import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.widgets import Slider

SAUVER = "--sauver" in sys.argv

# ---------------------------------------------------------------------------
# Paramètres (unités réduites LJ)
# ---------------------------------------------------------------------------

N = 30            # Nombre de particules
boxsize = 10.0    # Côté de la boîte
dt = 0.01         # Pas de temps
epsilon = 1.0     # Profondeur du puits de potentiel
sigma = 1.0       # Distance où le potentiel s'annule
temp_ini = 0.5    # Température cible initiale

np.random.seed(42)  # Graine fixe : même état initial à chaque lancement

# Positions aléatoires, vitesses gaussiennes (variance = température)
positions = np.random.rand(N, 2) * boxsize
vitesses = np.random.normal(0, np.sqrt(temp_ini), (N, 2))
accelerations = np.zeros((N, 2))

# Historiques pour les courbes d'énergie et la pression moyenne
hist = 200
ep_hist, ec_hist, et_hist = [0] * hist, [0] * hist, [0] * hist
P_hist = [0.0] * hist
x_list = list(range(hist))

v_bins = np.linspace(0, 4, 15)
v_centers = (v_bins[:-1] + v_bins[1:]) / 2

# Trajectoires de quelques particules pour la visualisation
nbr_traj = 5
long_traj = 400
traj_history = np.zeros((nbr_traj, long_traj, 2))
for i in range(nbr_traj):
    traj_history[i, :, :] = positions[i]

# ---------------------------------------------------------------------------
# Interface graphique
# ---------------------------------------------------------------------------

plt.rcParams["toolbar"] = "None"
fig = plt.figure(figsize=(11, 8))
gs = fig.add_gridspec(4, 2, width_ratios=[1.2, 1], height_ratios=[0.5, 1, 1, 0.1])

ax_sim = fig.add_subplot(gs[0:3, 0])
ax_sim.set_xlim(0, boxsize)
ax_sim.set_ylim(0, boxsize)
ax_sim.set_aspect("equal")
ax_sim.set_xticks([])
ax_sim.set_yticks([])
ax_sim.set_title("Gaz de Lennard-Jones 2D")

# Particules colorées selon la norme de leur vitesse
scat = ax_sim.scatter(positions[:, 0], positions[:, 1], s=150, c=np.linalg.norm(vitesses, axis=1),
                      cmap="coolwarm", edgecolors="black", vmin=0, vmax=3)
lines_traj = [ax_sim.plot([], [], "-", linewidth=1.5, alpha=0.7)[0] for _ in range(nbr_traj)]

ax_text = fig.add_subplot(gs[0, 1])
ax_text.axis("off")
time_text = ax_text.text(0.0, 0.5, "", fontsize=9, family="monospace", verticalalignment="center")

ax_graph = fig.add_subplot(gs[1, 1])
ax_graph.set_xlim(0, hist)
ax_graph.set_title("Énergies (sur les 200 derniers pas de temps)", fontsize=10)
ax_graph.grid(True, linestyle="--", alpha=0.5)
line_ec, = ax_graph.plot([], [], color="red", label="Ec", linewidth=1.5)
line_ep, = ax_graph.plot([], [], color="green", label="Ep", linewidth=1.5)
line_et, = ax_graph.plot([], [], color="black", linestyle=":", label="Etot", linewidth=1)
ax_graph.legend(loc="upper right", fontsize=8)

ax_hist = fig.add_subplot(gs[2, 1])
ax_hist.set_xlim(0, 4)
ax_hist.set_ylim(0, 1.5)
ax_hist.set_title("Vitesses (Maxwell-Boltzmann)", fontsize=10)
ax_hist.grid(True, linestyle="--", alpha=0.5)
line_hist, = ax_hist.plot([], [], drawstyle="steps-mid", color="blue", label="Simulation")
line_mb, = ax_hist.plot([], [], "r--", label="Théorie")
ax_hist.legend(loc="upper right", fontsize=8)

# Curseur pour changer la température cible en direct
ax_slider = fig.add_subplot(gs[3, :])
slider_temp = Slider(ax_slider, "Temp", 0.0, 5.0, valinit=temp_ini, color="orange")


def update_slider(val):
    global temp_ini
    temp_ini = slider_temp.val


slider_temp.on_changed(update_slider)

# ---------------------------------------------------------------------------
# Boucle d'intégration
# ---------------------------------------------------------------------------

def update(frame):
    global positions, vitesses, accelerations, ep_hist, ec_hist, et_hist, P_hist

    # Verlet-vitesse (partie 1) : demi-pas sur v, pas entier sur x
    vitesses += 0.5 * accelerations * dt
    positions += vitesses * dt

    # Conditions aux limites : murs réfléchissants
    for i in range(N):
        if positions[i, 0] < 0:
            positions[i, 0] = -positions[i, 0]
            vitesses[i, 0] *= -1
        elif positions[i, 0] > boxsize:
            positions[i, 0] = 2 * boxsize - positions[i, 0]
            vitesses[i, 0] *= -1

        if positions[i, 1] < 0:
            positions[i, 1] = -positions[i, 1]
            vitesses[i, 1] *= -1
        elif positions[i, 1] > boxsize:
            positions[i, 1] = 2 * boxsize - positions[i, 1]
            vitesses[i, 1] *= -1

    accelerations.fill(0)
    Ep = 0
    viriel = 0

    # Forces de Lennard-Jones, sur chaque paire (i < j)
    for i in range(N):
        for j in range(i + 1, N):
            dx = positions[i, 0] - positions[j, 0]
            dy = positions[i, 1] - positions[j, 1]
            dist_sq = dx ** 2 + dy ** 2
            r = np.sqrt(dist_sq)

            # Sécurité : on plafonne la répulsion si deux particules se
            # chevauchent (évite l'explosion numérique au démarrage)
            if r < 0.8:
                r = 0.8
                dist_sq = r * r

            # On part de (sigma/r)^2 pour éviter les puissances coûteuses
            sr2 = (sigma ** 2) / dist_sq
            sr6 = sr2 ** 3
            sr12 = sr6 ** 2

            # facteur_force = |F| / r = -(1/r) dV/dr, donc Fx = facteur_force * dx
            facteur_force = (24 * epsilon / dist_sq) * (2 * sr12 - sr6)
            fx = facteur_force * dx
            fy = facteur_force * dy

            # Action-réaction (3e loi de Newton) : une seule évaluation par paire
            accelerations[i, 0] += fx
            accelerations[i, 1] += fy
            accelerations[j, 0] -= fx
            accelerations[j, 1] -= fy

            Ep += 4 * epsilon * (sr12 - sr6)
            # Terme du viriel r_ij . F_ij = facteur_force * r^2
            viriel += facteur_force * dist_sq

    # Verlet-vitesse (partie 2)
    vitesses += 0.5 * accelerations * dt

    # Grandeurs thermodynamiques
    Ec = 0.5 * np.sum(vitesses ** 2)
    # Équipartition en 2D : Ec = N kB T, donc T = Ec / N en unités réduites
    temp_actu = Ec / N

    # Pression par le théorème du viriel en 2D : P = (Ec + W/2) / Aire
    P_inst = (Ec + 0.5 * viriel) / (boxsize ** 2)
    P_hist.pop(0)
    P_hist.append(P_inst)
    P_moy = np.mean(P_hist)

    # Thermostat de Berendsen (le facteur 0.1 amortit la correction)
    if temp_actu > 0:
        scale = np.sqrt(temp_ini / temp_actu)
        vitesses *= 1.0 + 0.1 * (scale - 1.0)

    # Mise à jour des graphiques
    Et = Ec + Ep
    ep_hist.pop(0)
    ep_hist.append(Ep)
    ec_hist.pop(0)
    ec_hist.append(Ec)
    et_hist.pop(0)
    et_hist.append(Et)

    norme_vit = np.linalg.norm(vitesses, axis=1)
    scat.set_offsets(positions)
    scat.set_array(norme_vit)

    line_ep.set_data(x_list, ep_hist)
    line_ec.set_data(x_list, ec_hist)
    line_et.set_data(x_list, et_hist)
    ax_graph.set_ylim(min(min(ep_hist), min(ec_hist)) - 5, max(max(ep_hist), max(ec_hist)) + 5)

    # density=True normalise l'histogramme pour le comparer à la densité théorique
    counts, _ = np.histogram(norme_vit, bins=v_bins, density=True)
    line_hist.set_data(v_centers, counts)

    # Maxwell-Boltzmann 2D : f(v) = (v / T) exp(-v^2 / 2T)
    if temp_actu > 0.01:
        f_v_theorie = (v_bins / temp_actu) * np.exp(-(v_bins ** 2) / (2 * temp_actu))
        line_mb.set_data(v_bins, f_v_theorie)
    else:
        line_mb.set_data([], [])

    for i in range(nbr_traj):
        traj_history[i, 1:] = traj_history[i, :-1]
        traj_history[i, 0] = positions[i]
        lines_traj[i].set_data(traj_history[i, :, 0], traj_history[i, :, 1])

    info = (f"Temp cible   : {temp_ini:.2f}\n"
            f"Temp réelle  : {temp_actu:.2f}\n"
            f"P inst       : {P_inst:.2f}\n"
            f"P moy        : {P_moy:.2f}\n\n"
            f"Ec   : {Ec:.1f}\n"
            f"Ep   : {Ep:.1f}\n"
            f"Etot : {Et:.1f}")
    time_text.set_text(info)

    return scat, time_text, line_ep, line_ec, line_et, line_hist, line_mb, *lines_traj


fig.tight_layout()

if SAUVER:
    # Démonstration sans souris, en quatre étapes de température :
    # condensation, liquide, gaz, puis recondensation.
    # On n'enregistre pas chaque pas de temps pour garder un GIF léger.
    os.makedirs("figures", exist_ok=True)
    scenario = [  # (température cible, nombre d'images, pas de temps par image, légende)
        (0.2, 100, 20, "T = 0.2 : condensation en amas"),
        (0.5, 75, 20, "T = 0.5 : l'amas fond et s'évapore en partie"),
        (5.0, 150, 4, "T = 5 : gaz"),
        (0.2, 100, 20, "Retour à T = 0.2 : recondensation"),
    ]
    etapes = []  # une entrée par image du GIF
    for T, n, pas, legende in scenario:
        etapes += [(T, pas, legende)] * n

    def demo(k):
        T, pas, legende = etapes[k]
        if k == 0 or etapes[k - 1] != etapes[k]:
            slider_temp.set_val(T)
            ax_sim.set_title(f"Gaz de Lennard-Jones 2D\n{legende}")
        for _ in range(pas - 1):
            update(k)
        return update(k)

    ani = animation.FuncAnimation(fig, demo, frames=len(etapes), blit=False)
    ani.save("figures/gaz_lennard_jones.gif", writer=animation.PillowWriter(fps=20), dpi=55)
    print("Animation enregistrée : figures/gaz_lennard_jones.gif")
else:
    ani = animation.FuncAnimation(fig, update, frames=200, interval=20, blit=False)
    plt.show()
