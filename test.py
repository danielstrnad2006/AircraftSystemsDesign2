import math
import numpy as np
from scipy.optimize import fsolve
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt
from wewillfindaname import * 
from wing import *




stiffeners = [
    (0.3, 'up', 400),
    (0.3, 'down', 400),
    (0.4, 'down', 400),
    (0.5, 'down', 400),
    (0.6, 'down', 400),
]

wing = WingSizing

cross_section1 = CrossSection(xc_spar1=0.2, xc_spar2=0.6, chord=500,
                              t_spar1=5, t_spar2=5,
                              t_skin_up=5, t_skin_down=5, stiffeners=stiffeners,
                              filepath="airfoils/NASA SC(2)-0414.dat")

cross_section1.assembly_centroid_finder()