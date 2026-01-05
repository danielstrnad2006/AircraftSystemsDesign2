import math
from planform import * 
import scipy as sp

def Shear_buckling(ribs_input_lst, shear_stresses_input_lst, thicknesses_input):

    ribs= ribs_input_lst[:-1]  #rib locations along span of beam
    shear_stresses= shear_stresses_input_lst  #maximum shear stress along that section between the ribs MPA
    thicknesses= [thicknesses_input*1e-3] * (len(ribs)) #metres

    wing = WingSizing(S_w=149.9, b=32.1632, c_root=7.2855, c_tip=2.0372, taper_ratio=0.2796,
                    leading_sweep=37.537, quart_sweep=34.4871, dihedral=5)

    k_s_values = [
            (1, 1.05, 15),
            (1.05, 1.1, 13.8),
            (1.1, 1.2, 13.2),
            (1.2, 1.3, 12.6),
            (1.3, 1.4, 12),
            (1.4, 1.5, 11.8),
            (1.5, 1.6, 11.5),
            (1.6, 1.7, 11.2),
            (1.7, 1.8, 10.9),
            (1.8, 1.9, 10.7),
            (1.9, 2.0, 10.5),
            (2.0, 2.1, 10.4),
            (2.1, 2.2, 10.2),
            (2.2, 2.3, 10.0),
            (2.3, 2.4, 9.8),
            (2.4, 2.5, 9.75),
            (2.5, 2.6, 9.7),
            (2.6, 3.1, 9.6),
            (3.1, 3.3, 9.55),
            (3.3, 3.6, 9.5),
            (3.6, 5.0, 9.4) 
        ]

    def aspect_to_K(aspect_ratio):
        for i in k_s_values:
            if i[0] <= aspect_ratio < i[1]:
                return i[2]
            elif aspect_ratio>=k_s_values[-1][1]:
                return 9.4
            elif aspect_ratio<k_s_values[0][0]:
                return 15
            

    def shear_buckling_along(ribs, shear_stresses, thicknesses):
        output=[]
        x_data=[0]
        poisson = 0.33
        E_mod = 72.4*10**9
        for i in range(len(ribs)):
            if i == len(ribs) - 1:
                length=wing.b/2-ribs[-1]
            else:
                length=ribs[i+1]-ribs[i]
            x_data.append(length+x_data[-1])
            if i!=0:
                b_spar =  0.14 * wing.chord(ribs[i-1])  #conservative assumption, takes the chord length at that location
            else: b_spar =0.14*wing.chord(0)
            aspect_ratio = (length/b_spar)
            K_value=aspect_to_K(aspect_ratio)
            tau_crit = (math.pi**2*K_value*E_mod/(12*(1-poisson**2))*(thicknesses[i]/b_spar)**2)*10**(-6)
            safety=tau_crit/shear_stresses[i]
            if safety>20:safety=20
            output.append(safety)
        x_data.pop(-1)
        
        graph = sp.interpolate.interp1d(x_data,output,kind="previous",bounds_error=False,fill_value=(output[0], output[-1]))
        return output, graph

    ##input all of this to find the graphs for safety buckling 


    #safety_Factor=1.2
    #count=50
    #tolerance=0.05
    #for _ in range(count):
    #    values, useless = shear_buckling_along(ribs, shear_stresses, thicknesses)
#
    #    for k, MS in enumerate(values):
#
    #        if MS > safety_Factor + tolerance:
    #            thicknesses[k] -= 0.001
    #        elif MS < safety_Factor:
    #            thicknesses[k] += 0.001

    
    return shear_buckling_along(ribs, shear_stresses, thicknesses)[0]

