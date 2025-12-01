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

wing = WingSizing(S_w=149.9, b=b, c_root=7.2855, c_tip=2.0372, taper_ratio=0.2796,
                 leading_sweep=37.537, quart_sweep=34.4871, dihedral=5)

CENTROID_X = []
CENTROID_Y = []
I_XX = []
I_YY = []
J_P = []

def interp(object, location):
    y = np.arange(0, b / 2, db)
    return np.interp(location, y, object)

while b_cur < wing.b/2:
    stiffeners = [
        #(0.3, 'up', 40000, 10),
        #(0.3, 'down', 40000, 10),
        #(0.5, 'up', 40000, 10),
        #(0.5, 'down', 40000, 10),
    ]

    wing.plot(chord_position=b_cur)

    cs = CrossSection(xc_spar1=0.2, xc_spar2=0.6, chord=wing.chord(b_cur)*1000, b_cur=b_cur,
                                t_spar1=10, t_spar2=10,
                                t_skin_up=10, t_skin_down=10, stiffeners=stiffeners,
                                filepath="airfoils/NASA SC(2)-0414.dat", display_data=True)

    centroid_X, centroid_Y, I_xx, I_yy, J_p = cs.assembly_centroid_finder()
    b_pos.append(b_cur)
    CENTROID_X.append(centroid_X)
    CENTROID_Y.append(centroid_Y)
    I_XX.append(I_xx)
    I_YY.append(I_yy)
    J_P.append(J_p)

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