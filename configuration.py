import math
import numpy as np
from scipy.optimize import fsolve
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt
from cross_section import * 
from planform import *
import json


db = 0.5
b = 32.1632
b_cur = 0

b_pos = []
cross_sections = []

PLOT = False


wing = WingSizing(S_w=149.9, b=b, c_root=7.2855, c_tip=2.0372, taper_ratio=0.2796,
                 leading_sweep=37.537, quart_sweep=34.4871, dihedral=5)

CENTROID_X = []
CENTROID_Y = []
I_XX = []
I_YY = []
J_P = []
Q = []

def interp(object, location):
    y = np.arange(0, b / 2, db)
    return np.interp(location, y, object)

while b_cur < wing.b/2:
    stiffeners = [
       (0.25, 'up', 90, 10),
       (0.25, 'down', 90, 10),
       (0.30, 'up', 90, 10),
       (0.30, 'down', 90, 10),
       (0.35, 'up', 90, 10),
       (0.35, 'down', 90, 10),
       (0.40, 'up', 90, 10),
       (0.40, 'down', 90, 10),
       (0.45, 'up', 90, 10),
       (0.45, 'down', 90, 10),
       (0.50, 'up', 90, 10),
       (0.50, 'down', 90, 10),
       (0.55, 'up', 90, 10),
       (0.55, 'down', 90, 10),
       (0.60, 'up', 90, 10),
       (0.60, 'down', 90, 10),
    ]

    if PLOT:
        wing.plot(chord_position=b_cur)

    t_spar1 = 6
    t_spar2 = 6
    t_skin_up = 6
    t_skin_down = 6
    cs = CrossSection(xc_spar1=0.2, xc_spar2=0.65, chord=wing.chord(b_cur)*1000, b_cur=b_cur,
                                t_spar1=t_spar1, t_spar2=t_spar2,
                                t_skin_up=t_skin_up, t_skin_down=t_skin_down, stiffeners=stiffeners,
                                filepath="airfoils/NASA SC(2)-0414.dat", save_plot=PLOT)

    centroid_X, centroid_Y, I_xx, I_yy, J_p = cs.assembly_centroid_finder()
    b_pos.append(b_cur)
    CENTROID_X.append(centroid_X)
    CENTROID_Y.append(centroid_Y)
    I_XX.append(I_xx)
    I_YY.append(I_yy)
    J_P.append(J_p)
    Q.append(cs.find_Q())

    cross_sections.append(cs)

    b_cur = b_cur + db


plt.figure()

plt.plot(b_pos, CENTROID_X, label="Centroid X")
plt.plot(b_pos, CENTROID_Y, label="Centroid Y")

plt.xlabel("Position along blade (b_pos)")
plt.ylabel("Centroid coordinate")
plt.title("Centroid Position Along")
plt.legend()

plt.grid(True)
plt.show()

plt.figure()

plt.plot(b_pos, I_XX, label="Ixx")
plt.plot(b_pos, I_YY, label="Iyy")
plt.plot(b_pos, J_P, label="J_P")

plt.xlabel("Position along blade (b_pos)")
plt.ylabel("Moment of Inertia")
plt.title("Sectional Moments of Inertia Along Blade")
plt.legend()

plt.grid(True)
plt.show()

print(I_XX)
print(CENTROID_X)
print(J_P)

with open("i_xx.txt", "w") as f:
    json.dump(I_XX, f)

with open("centroid_x.txt", "w") as f:
    json.dump(CENTROID_X, f)

with open("j_p.txt", "w") as f:
    json.dump(J_P, f)

with open("Q.txt", "w") as f:
    json.dump(Q, f)

with open("spar_thickness.txt", "w") as f:
    spar_thickness = [t_spar1, t_spar2]
    json.dump(spar_thickness, f)

with open("stiffeners_areas.txt", "w") as f:
    stiffeners_areas = [stiffener[2] * 1e-6 for stiffener in stiffeners]
    json.dump(stiffeners_areas, f)