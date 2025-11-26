# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import internal_loads
from Wing_twist_deflection import  twist_plot, wing_tip_twist
from Bending_Deflection import BeamDeflection
import math
import numpy as np
from scipy.optimize import fsolve
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt
from cross_section import * 
from planform import *

db = 0.5
b = 32.1632

"""
internal_properties=internal_loads.halfWing
internal_properties.set_conditions(2.5, 103544, 120, 1.225, 100)
internal_properties.torque_plot()
twist_plot(internal_properties)
internal_properties.get_internal_plot()
moment_distribution = internal_properties.internal_bending
"""

#PLOT BEAM DEFLECTION
beam = BeamDeflection(b, db)

I_XX = [np.float64(33480308035.556774), np.float64(31856908251.130173), np.float64(30279526734.80252), np.float64(28747742163.1427), np.float64(27261133202.22146), np.float64(25819278507.381264), np.float64(24421756722.999924), np.float64(23068146482.247917), np.float64(21758026406.839333), np.float64(20490975106.77606), np.float64(19266571180.085144), np.float64(18084393212.549038), np.float64(16944019777.428473), np.float64(15845029435.177847), np.float64(14787000733.152805), np.float64(13769512205.30965), np.float64(12792142371.896643), np.float64(11854469739.136549), np.float64(10956072798.900415), np.float64(10096530028.372175), np.float64(1053064031.2885389), np.float64(930135986.3936577), np.float64(817167263.0431123), np.float64(713737138.6943815), np.float64(619424890.8049406), np.float64(533809796.8322648), np.float64(456471134.23382986), np.float64(386988180.4671123), np.float64(324940212.9895887), np.float64(269906509.25873303), np.float64(221466346.73202258), np.float64(179199002.8669331), np.float64(142683755.12093985)]

beam.assignI_XX(I_XX)
beam.plot()

max_v, max_y = beam.max_deflection()
print(f"Maximum deflection: {max_v:.6f} at y = {max_y:.3f}")
