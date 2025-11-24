import math
import numpy as np
from scipy.optimize import fsolve
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt
from cross_section import * 
from planform import *


db = 1
b = 32.1632
b_cur = 0
cross_sections = []

wing = WingSizing(S_w=149.9, b=b, c_root=7.2855, c_tip=2.0372, taper_ratio=0.2796,
                 leading_sweep=37.537, quart_sweep=34.4871, dihedral=5)

CENTROID_X = []
CENTROID_Y = []
I_XX = []
I_YY = []

def interp(object, location):
    y = np.arange(0, b / 2, db)
    return np.interp(location, y, object)

while b_cur < wing.b/2:
    stiffeners = [
        (0.3, 'up', 400, 10),
        (0.4, 'up', 400, 5),
        (0.5, 'up', 400, 10),
        (0.3, 'down', 400, 8),
        (0.4, 'down', 400, 3),
        (0.5, 'down', 400, 8),
    ]

    wing.plot(chord_position=b_cur)

    cs = CrossSection(xc_spar1=0.2, xc_spar2=0.6, chord=wing.chord(b_cur)*1000, b_cur=b_cur,
                                t_spar1=50, t_spar2=50,
                                t_skin_up=50, t_skin_down=50, stiffeners=stiffeners,
                                filepath="airfoils/NASA SC(2)-0414.dat", display_data=False)

    centroid_X, centroid_Y, I_xx, I_yy = cs.assembly_centroid_finder()

    CENTROID_X.append(centroid_X)
    CENTROID_Y.append(centroid_Y)
    I_XX.append(I_xx)
    I_YY.append(I_yy)

    cross_sections.append(cs)

    b_cur = b_cur + db
