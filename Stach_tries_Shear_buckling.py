import math
from planform import * 
import scipy as sp

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
ribs= [0, 2, 4, 6, 7.03893, 9, 12]  #rib locations along span of beam
shear_stresses= [np.float64(1334.2858661244402), np.float64(1142.9765962870895), np.float64(949.8923116195851), np.float64(765.7974376040061), np.float64(675.4898938320113), np.float64(457.8899650787215), np.float64(177.3924126996869)]  #maximum shear stress along that section between the ribs MPA
thicknesses=[0.05,0.05,0.05,0.03,0.03,0.02,0.02] #metres


safety_Factor=1.2
count=50
tolerance=0.05
for _ in range(count):
    values, useless = shear_buckling_along(ribs, shear_stresses, thicknesses)

    for k, MS in enumerate(values):

        if MS > safety_Factor + tolerance:
            thicknesses[k] -= 0.001
        elif MS < safety_Factor:
            thicknesses[k] += 0.001

values,graph= shear_buckling_along(ribs, shear_stresses, thicknesses)
x_plot = np.linspace(0, wing.b/2, 300)
sf_plot = graph(x_plot)
plt.plot(x_plot, sf_plot)
plt.xlabel("Spanwise location y [m]")
plt.ylabel("Shear Buckling Margin Of Safety[-]")
plt.grid(True)
plt.show()

print(thicknesses)

