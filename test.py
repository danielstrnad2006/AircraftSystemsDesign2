import math
import numpy as np
from scipy.optimize import fsolve
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt
from cross_section import * 
from wing import *


wing = WingSizing(S_w=149.9, b=32.1632, c_root=7.2855, c_tip=2.0372, taper_ratio=0.2796,
                 leading_sweep=37.537, quart_sweep=34.4871, dihedral=5)


db = 1
b_cur = 0

while b_cur < wing.b/2:

    stiffeners = [
        (0.3, 'up', 400),
        (0.4, 'up', 400),
        (0.5, 'up', 400),
        (0.3, 'down', 400),
        (0.4, 'down', 400),
        (0.5, 'down', 400),
    ]

    cross_section1 = CrossSection(xc_spar1=0.2, xc_spar2=0.6, chord=wing.chord(b_cur)*1000,
                                t_spar1=50, t_spar2=50,
                                t_skin_up=50, t_skin_down=50, stiffeners=stiffeners,
                                filepath="airfoils/NASA SC(2)-0414.dat")

    cross_section1.assembly_centroid_finder()

    b_cur = b_cur + db
