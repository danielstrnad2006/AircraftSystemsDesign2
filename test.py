import math
import numpy as np
from scipy.optimize import fsolve
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt
from wewillfindaname import * 

cross_section1 = CrossSection(xc_spar1=0.2, xc_spar2=0.6, chord=200,
                              t_spar1=10, t_spar2=10,
                              t_skin_up=5, t_skin_down=5,
                              filepath="airfoils/NASA SC(2)-0414.dat")

cross_section1.assembly_centroid_finder()
cross_section1.plot()