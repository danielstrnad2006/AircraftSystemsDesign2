import math
from planform import * 

db = 0.0005      # m
b = 32.1632      # m
b_cur = 0

wing = WingSizing(S_w=149.9, b=b, c_root=7.2855, c_tip=2.0372, taper_ratio=0.2796,
                 leading_sweep=37.537, quart_sweep=34.4871, dihedral=5)

ribs = [0,3,6,9,12]
b_rib= 0      #m
drib = 0.0005 #m


while b_cur < wing.b/2:  
    poisson = 0.33
    k_s = 1 ##  !!! dummy value ###
    E = 72.4e9    # GPA
    t = 2.5e-3    # m, design value
    b_spar =  0.14 * wing.chord(b)    # m, short side of the plate, height of the spar

    i = 0
    for i in range(len(ribs)):
        a =  float(ribs(i+1) - ribs(i))  # m, distance between ribs


    aspect_ratio = a/b_spar

    k_s_values = [
        (1,   1.05,  15),
        (1.05, 1.1,  13.8),
        (1.1, 1.2,  13.2),
        (1.2, 1.3,  12.6),
        (1.3, 1.4,  12),
        (1.4, 1.5,  11.8),
        (1.5, 1.6,  11.5),
        (1.6, 1.7,  11.2),
        (1.7, 1.8,  10.9),
        (1.8, 1.9,  10.7),
        (1.9, 2.0,  10.5),
        (2.0, 2.1,  10.4),
        (2.1, 2.2,  10.2),
        (2.2, 2.3,  10.0),
        (2.3, 2.4,  9.8),
        (2.4, 2.5,  9.75),
        (2.5, 2.6,  9.7),
        (2.6, 3.1,  9.6),
        (3.1, 3.3,  9.55),
        (3.3, 3.6,  9.5),
        (3.6, 5.0,  9.4) 
    ]

    def get_value(aspect_ratio):
        for low, high, val in k_s_values:
            if low <= aspect_ratio < high:
                return val




    k_s = float(get_value(aspect_ratio))

    tau_crit = math.pi**2*k_s*E/(12*(1-poisson**2))*(t/b_spar)**2
    b_cur = b_cur + db

print("tau critical:", tau_crit)




k_s_values = [
        (1,   1.05,  15),
        (1.05, 1.1,  13.8),
        (1.1, 1.2,  13.2),
        (1.2, 1.3,  12.6),
        (1.3, 1.4,  12),
        (1.4, 1.5,  11.8),
        (1.5, 1.6,  11.5),
        (1.6, 1.7,  11.2),
        (1.7, 1.8,  10.9),
        (1.8, 1.9,  10.7),
        (1.9, 2.0,  10.5),
        (2.0, 2.1,  10.4),
        (2.1, 2.2,  10.2),
        (2.2, 2.3,  10.0),
        (2.3, 2.4,  9.8),
        (2.4, 2.5,  9.75),
        (2.5, 2.6,  9.7),
        (2.6, 3.1,  9.6),
        (3.1, 3.3,  9.55),
        (3.3, 3.6,  9.5),
        (3.6, 5.0,  9.4) 
    ]

def aspect_to_K(aspect_ratio):
    for i in k_s_values:
        if i[0] <= aspect_ratio < i[1]:
                return i[2]
        

def shear_buckling_along(rib_loc, shear_stresses, thicknesses):
    output=[]
    poisson2 = 0.33
    E_mod = 72.4*10**9
    for i in range(len(rib)-1):
        length=rib_loc[i+1]-rib_loc[i]
        b_spar =  0.14 * wing.chord(rib_loc[i])  #conservative assumption, takes the chord length at that location
        aspect_ratio = length/b_spar
        K_value=aspect_to_K(aspect_ratio)

        tau_crit = math.pi**2*K_value*E_mod/(12*(1-poisson2**2))*(thicknesses[i]/b_spar)**2
        output.append(tau_crit/shear_stresses[i])
    return output #returns list of safety factors
        
ribs= [0, 3, 6, 9, 12] #rib locations along span of beam
shear_stresses=[10**9,2*10**9,3*10**9,4*10**9,5*10**9]  #maximum shear stress along that section between the ribs
thicknesses=[0.001,0.001,0.001,0.001,0.001]

print(shear_buckling_along(ribs,shear_stresses, thicknesses))







shear_force = 1000 # N, shear force
h_f =   0.2    # m, height of the front spar    !!!dummy value
t_f =  2.5e-3  # m, thickness of the front spar !!! design dependant
h_r =   0.2    # m, height of the rear spar     !!!dummy value
t_r =  2.5e-3  # m, thickness of the rear spar  !!! design dependant
k_v =  1.2     # !!!!dummy value!!! 

tau_ave_shear = shear_force / (h_f*t_f + h_r*t_r)
tau_max_shear = k_v * tau_ave_shear

print("tau average shear:", tau_ave_shear)
print("tau max. shear:", tau_max_shear)