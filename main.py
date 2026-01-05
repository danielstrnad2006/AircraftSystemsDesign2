# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import internal_loads
from Wing_twist_deflection import  twist_plot, wing_tip_twist
from Bending_Deflection import BeamDeflection
from Daniel_test_Torsion_Deflection import TorsionDeflection, plotTwists

import math
import numpy as np
import scipy as sp
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from cross_section import * 

from planform import *
import json

from Stach_tries_Shear_buckling import Shear_buckling
from column_buckling import Column_buckling
from skin_buckling import Skin_buckling

db = 0.5
b = 32.1632

internal_properties=internal_loads.halfWing
crit_conds = [[2.5, 103544, 308.7, 1.225, 100],[-1, 103544,  87.3, 0.433, 100], [2.5,  43807, 308.7, 1.225,   0], [2.5, 103544, 138.0, 1.225, 100]]

ribs_locations = [0, 0.5, 1, 1.5, 3.5, 4, 5, 6, 7.03893, 8, 9, 12, 16.0816]

#I_XX = [np.float64(33480308035.556774), np.float64(31856908251.130173), np.float64(30279526734.80252), np.float64(28747742163.1427), np.float64(27261133202.22146), np.float64(25819278507.381264), np.float64(24421756722.999924), np.float64(23068146482.247917), np.float64(21758026406.839333), np.float64(20490975106.77606), np.float64(19266571180.085144), np.float64(18084393212.549038), np.float64(16944019777.428473), np.float64(15845029435.177847), np.float64(14787000733.152805), np.float64(13769512205.30965), np.float64(12792142371.896643), np.float64(11854469739.136549), np.float64(10956072798.900415), np.float64(10096530028.372175), np.float64(1053064031.2885389), np.float64(930135986.3936577), np.float64(817167263.0431123), np.float64(713737138.6943815), np.float64(619424890.8049406), np.float64(533809796.8322648), np.float64(456471134.23382986), np.float64(386988180.4671123), np.float64(324940212.9895887), np.float64(269906509.25873303), np.float64(221466346.73202258), np.float64(179199002.8669331), np.float64(142683755.12093985)]
#CENTROID_X = [np.float64(2910.365816629366), np.float64(2845.284184294819), np.float64(2780.2000091008904), np.float64(2715.1132578730003), np.float64(2650.02389685698), np.float64(2584.9318917063597), np.float64(2519.8372074693193), np.float64(2454.739808575291), np.float64(2389.6396588212087), np.float64(2324.536721357386), np.float64(2259.4309586730137), np.float64(2194.3223325812705), np.float64(2129.2108042040186), np.float64(2064.0963339560954), np.float64(1998.9788815291624), np.float64(1933.8584058751196), np.float64(1868.7348651890532), np.float64(1803.6082168917148), np.float64(1738.478417611508), np.float64(1673.345423165974), np.float64(1598.4458874721213), np.float64(1533.657895028458), np.float64(1468.8699025847939), np.float64(1404.0819101411303), np.float64(1339.293917697467), np.float64(1274.5059252538033), np.float64(1209.7179328101392), np.float64(1144.9299403664754), np.float64(1080.141947922812), np.float64(1015.3539554791481), np.float64(950.5659630354845), np.float64(885.7779705918211), np.float64(820.9899781481573)]
#J_P = [np.float64(69483201639.66885), np.float64(64920647095.20532), np.float64(60562326996.2465), np.float64(56403564783.54618), np.float64(52439683897.85787), np.float64(48666007779.93528), np.float64(45077859870.531944), np.float64(41670563610.40155), np.float64(38439442440.297714), np.float64(35379819800.97405), np.float64(32487019133.18418), np.float64(29756363877.681732), np.float64(27183177475.220356), np.float64(24762783366.553623), np.float64(22490504992.43521), np.float64(20361665793.618725), np.float64(18371589210.857788), np.float64(16515598684.906036), np.float64(14789017656.51707), np.float64(13187169566.444551), np.float64(11705377855.442066), np.float64(10338965964.263279), np.float64(9083257333.661774), np.float64(7933575404.391199), np.float64(6885243617.205185), np.float64(5933585412.857342), np.float64(5073924232.101298), np.float64(4301583515.690682), np.float64(3611886704.3791294), np.float64(3000157238.920244), np.float64(2461718560.0676656), np.float64(1991894108.5750177), np.float64(1586007325.1959202)]

with open("i_xx.txt") as f:
    I_XX = json.load(f)
with open("centroid_x.txt") as f:
    CENTROID_X = json.load(f)
with open("j_p.txt") as f:
    J_P = json.load(f)
with open("Q.txt") as f:
    Q = json.load(f)
with open("spar_thickness.txt") as f:
    spar_thickness = json.load(f)
with open("stiffeners_areas.txt") as f:
    stiffeners_areas = json.load(f)

internal_properties.set_torsion_params(db, CENTROID_X, J_P)
internal_properties.set_buckling_params(db, Q, I_XX, spar_thickness=spar_thickness, ribs_locations = ribs_locations)
compressive_yield_strength = 4e8 #Pa

if input("Start with detailed analysis of critical conditions? (y)")=="y":
    for cond in crit_conds:
        print(cond) 
        print("going through load case with load factor: ", cond[0], ", mass", cond [1], "kg, Equivalent Air Speed: ", cond[2], "m/s, and density", cond[3], "kg/m^3, fuel percentage of:", cond[4], "%")
        internal_properties.set_conditions(load_factor=cond[0], weight=cond[1], v_EAS=cond[2], rho=cond[3], fuel_percentage=cond[4])
        #internal_properties.get_coefficient_plots()
        #internal_properties.get_forces_plot()
        internal_properties.get_internal_plot()
        #internal_properties.get_debugging_torsion_plot()
        internal_properties.get_internal_torsion_plot()
        #internal_properties.torque_plot()

        


        #twist_plot(internal_properties, CENTROID_X, J_P, b, db)
        moment_distribution = internal_properties.internal_bending
        torsion_noT_distribution = internal_properties.internal_torsion_noT
        torsion_fullT_distribution = internal_properties.internal_torsion_fullT
       
        #PLOT BEAM DEFLECTION
        beam = BeamDeflection(b, db, moment_distrib=moment_distribution)
        beam.assignI_XX(I_XX)
        beam.plot()

        #PLOT BEAM TWIST
        beam_twist_noT = TorsionDeflection(b, db, torsion_distrib=torsion_noT_distribution, J_distrib=internal_properties.J)
        beam_twist_fullT = TorsionDeflection(b, db, torsion_distrib=torsion_fullT_distribution, J_distrib=internal_properties.J)
        plotTwists(theta_fullT=beam_twist_fullT.theta, theta_noT=beam_twist_noT.theta)


        sigma_lst = internal_properties.get_normal_stress_at_sections()
        sigma_interp1d = sp.interpolate.interp1d(ribs_locations[:-1], sigma_lst, kind='previous', bounds_error=False, fill_value=(sigma_lst[0], sigma_lst[-1]))
        tau_lst = internal_properties.get_shear_at_sections()

        print("List of maximum shear stresses along the span in MPa:", [int(tau/1e6) for tau in tau_lst])
        print("List of maximum normal stresses along the span in MPa:", [int(sigma/1e6) for sigma in sigma_lst])


        max_v, max_y = beam.max_deflection()
        print(f"Maximum deflection: {max_v:.6f} meters at y = {max_y:.3f} meters")

        ### Buckling Analysis ###
        shear_buckling_safety = Shear_buckling(ribs_locations, [tau/1e6 for tau in tau_lst], spar_thickness[0])
        skin_buckling_safety = Skin_buckling()[1]

        shear_buckling_safety_interp1d = shear_buckling_safety 
        skin_buckling_safety_interp1d = sp.interpolate.interp1d(ribs_locations[:-1], skin_buckling_safety, kind='previous', bounds_error=False, fill_value=(skin_buckling_safety[0], skin_buckling_safety[-1]))
        column_buckling_safety_interp1d = Column_buckling(ribs_locations, [sigma/1e6 for sigma in sigma_lst], stringer_areas_input=stiffeners_areas)
    
        input("Press Enter to continue to graphing of Margins of Safety...") 
        x_plot = np.linspace(0, b/2, 300)
        shear_buckling_safety_func = shear_buckling_safety_interp1d(x_plot)
        compressive_stress_safety_func = compressive_yield_strength/sigma_interp1d(x_plot)
        
        plt.plot(x_plot, shear_buckling_safety_func, label="Shear Buckling Margin Of Safety")
        plt.plot(x_plot, column_buckling_safety_interp1d(x_plot), label="Column Buckling Margin Of Safety")
        plt.plot(x_plot, compressive_stress_safety_func, label="Compressive Yielding Margin Of Safety")
        plt.plot(x_plot, skin_buckling_safety_interp1d(x_plot), label="Skin Buckling Margin Of Safety")
        plt.ylim(0, 8)
        plt.xlabel("Spanwise location y [m]")
        plt.ylabel("Margin Of Safety[-]")
        plt.grid(True)
        plt.hlines(y=1, color='r', linestyle='--', label="Failure Threshold", xmin=0, xmax=b/2)
        plt.title("Margin Of Safety Along Span for different Failure Modes")
        plt.legend()
        plt.show()





#if input("Is the final cross section chosen and do you want to proceed to verify deflection at all loading conditions? (y)")=="y":
#    crit_conds = [
#    [2.5, 103544, 138.0, 0.433, 100],   # LC-1 MTOW FL310
#    [2.5,  62767, 107.5, 0.433,   0],   # LC-2 ZFW  FL310
#    [2.5,  43807,  89.8, 0.433,   0],   # LC-3 OEW  FL310
#
#    [2.5, 103544, 138.0, 1.225, 100],   # LC-4 MTOW FL0
#    [2.5,  62767, 107.5, 1.225,   0],   # LC-5 ZFW  FL0
#    [2.5,  43807,  89.8, 1.225,   0],   # LC-6 OEW  FL0
#
#
#    # --- Negative manoeuvring @ VS1 ---
#    [-1, 103544,  87.3, 0.433, 100],    # LC-7 MTOW FL310
#    [-1,  62767,  67.7, 0.433,   0],    # LC-8 ZFW  FL310
#    [-1,  43807,  56.8, 0.433,   0],    # LC-9 OEW  FL310
#
#    [-1, 103544,  87.3, 1.225, 100],    # LC-10 MTOW FL0
#    [-1,  62767,  67.7, 1.225,   0],    # LC-11 ZFW  FL0
#    [-1,  43807,  56.8, 1.225,   0],    # LC-12 OEW  FL0
#
#
#    # --- Positive manoeuvring @ VD ---
#    [2.5, 103544, 163.0, 0.433, 100],   # LC-13 MTOW FL310
#    [2.5,  62767, 163.0, 0.433,   0],   # LC-14 ZFW  FL310
#    [2.5,  43807, 163.0, 0.433,   0],   # LC-15 OEW  FL310
#
#    [2.5, 103544, 308.7, 1.225, 100],   # LC-16 MTOW FL0
#    [2.5,  62767, 308.7, 1.225,   0],   # LC-17 ZFW  FL0
#    [2.5,  43807, 308.7, 1.225,   0],   # LC-18 OEW  FL0
#
#
#    # --- Positive manoeuvring @ VC ---
#    [2.5, 103544, 154.0, 0.433, 100],   # LC-19 MTOW FL310
#    [2.5,  62767, 154.0, 0.433,   0],   # LC-20 ZFW  FL310
#    [2.5,  43807, 154.0, 0.433,   0],   # LC-21 OEW  FL310
#
#    [2.5, 103544, 154.0, 1.225, 100],   # LC-22 MTOW FL0
#    [2.5,  62767, 154.0, 1.225,   0],   # LC-23 ZFW  FL0
#    [2.5,  43807, 154.0, 1.225,   0]    # LC-24 OEW  FL0
#    ]
#    twist_tip_lst = []
#    deflection_tip_lst = []
#
#    for cond in crit_conds:
#        print("going through load case with load factor: ", cond[0], ", mass", cond [1], "kg, Equivalent Air Speed: ", cond[2], "m/s, and density", cond[3], "kg/m^3")
#        internal_properties.set_conditions(load_factor=cond[0], weight=cond[1], v_EAS=cond[2], rho=cond[3], fuel_percentage=cond[4])
#        #internal_properties.get_internal_plot()
#        #internal_properties.get_internal_torsion_plot()
#
#        #twist_plot(internal_properties, CENTROID_X, J_P, b, db)
#        moment_distribution = internal_properties.internal_bending
#        torsion_noT_distribution = internal_properties.internal_torsion_noT
#        torsion_fullT_distribution = internal_properties.internal_torsion_fullT
#
#        #PLOT BEAM DEFLECTION
#        beam = BeamDeflection(b, db, moment_distrib=moment_distribution)
#        beam.assignI_XX(I_XX)
#        #beam.plot()
#
#        #PLOT BEAM TWIST
#        beam_twist_noT = TorsionDeflection(b, db, torsion_distrib=torsion_noT_distribution, J_distrib=internal_properties.J)
#        beam_twist_fullT = TorsionDeflection(b, db, torsion_distrib=torsion_fullT_distribution, J_distrib=internal_properties.J)
#        #plotTwists(theta_fullT=beam_twist_fullT.theta, theta_noT=beam_twist_noT.theta)
#        print("At load case with load factor: ", cond[0], ", mass", cond [1], "kg, Equivalent Air Speed: ", cond[2], "m/s, and density", cond[3], "kg/m^3")
#        print("The Reaction Shear Force is:", internal_properties.internal_shear(0))
#        print("The Reaction Bending Moment is:", internal_properties.internal_bending(0))
#        print("The Reaction Torsion Moment at full throttle is:", internal_properties.internal_torsion_fullT(0))
#        print("The Reaction Torsion Moment at zero throttle is:", internal_properties.internal_torsion_noT(0))
#        print()
#        print("The tip deflection is:", beam.v(16.0816)* 1e15, "mm")
#        deflection_tip_lst.append(beam.v(16.0816)* 1e12)
#        print("The tip twist (zero T) is:", beam_twist_noT.theta(16.0816)*180/np.pi, "deg")
#        print("The tip twist (full T) is:", beam_twist_fullT.theta(16.0816)*180/np.pi, "deg")
#        if beam_twist_fullT.theta(16.0816)>0:
#            twist_tip_lst.append(max(beam_twist_fullT.theta(16.0816)*180/np.pi, beam_twist_noT.theta(16.0816)*180/np.pi))
#        else: 
#            twist_tip_lst.append(min(beam_twist_noT.theta(16.0816)*180/np.pi, beam_twist_fullT.theta(16.0816)*180/np.pi))
#        print("-----------------------------------------------------")
#        print()
#        print()
#    print("The tip deflections for all loading conditions in m are:", deflection_tip_lst)
#    print("The tip twists for all loading conditions in deg are:", twist_tip_lst)
#