import math
import matplotlib.pyplot as plt
from planform import *





def Skin_buckling(ribs, normal_stresses, thickness_input):
    span = 32.1632
    poisson = 0.33
    E = 72.4e9
    ds = 0

    sigma_cr_tab = []
    safety_tab = []
    ds_tab = []
    t = thickness_input

    wing = WingSizing(S_w=149.9, b=span, c_root=7.2855, c_tip=2.0372, taper_ratio=0.2796,
                    leading_sweep=37.537, quart_sweep=34.4871, dihedral=5)


    def skin_buckling(b, kc):
        sigma_cr = ((((math.pi**2)*kc*E)/(12*(1-(poisson**2))))*((t/b)**2))/(10**6)
        return sigma_cr
    
    def skin_buckling_wide(a, kc):
        sigma_cr = ((((math.pi**2)*E)/(12*(1-(poisson**2))))*((t/a)**2))/(10**6)
        return sigma_cr
    
    def bending_stress(M, y, I_xx):
        sigma_ben = M * y / I_xx
        return sigma_ben

    def safety_margin(sigma_ben, sigma_cr):
        margin_of_safety = sigma_cr/sigma_ben
        return margin_of_safety

    c_root = wing.c_root
    c_tip = wing.c_tip
    b = wing.c_root
    a = wing.b/2

    print(f"c_root: {c_root}")    
    print(f"c_tip: {c_tip}")    

    i = 0
    for i in range(len(ribs)-1):

        if i == len(ribs) - 1:
            a=wing.b/2-ribs[-1]
        else:
            a=ribs[i+1]-ribs[i]

        b = (c_root - (2*(c_root - c_tip)/span) * ribs[i+1])*0.45

        kc_ratio = a/b

        kc_values = [
            (0.0, 0.1, 26.3),
            (0.1, 0.2, 24.7),
            (0.2, 0.3, 23.1),
            (0.3, 0.4, 21.5),
            (0.4, 0.5, 19.9),
            (0.5, 0.6, 18.3),          
            (0.6, 0.7, 16.7),
            (0.7, 0.8, 15.1),
            (0.8, 0.9, 13.5),
            (0.9, 1.0, 11.9),
            (1.0, 1.1, 10.3),
            (1.1, 1.15, 10.4),
            (1.15, 1.2, 10.5),
            (1.2, 1.3, 9.85),
            (1.3, 1.4, 9.66),
            (1.4, 1.5, 9.46),
            (1.5, 1.6, 9.31),
            (1.6, 1.7, 9.15),
            (1.7, 1.8, 8.58),
            (1.8, 1.9, 8),
            (1.9, 2, 7.96),
            (2, 2.1, 7.92),       
            (2.1, 2.6, 7.90),
            (2.6, 2.7, 7.77),                           
        ]

        def get_value(kc_ratio):
            for low, high, val in kc_values:
                if low <= kc_ratio < high:
                    return val
                
            if kc_ratio >= 2.7:
                return 7.5
            else:
                return 20
                        
        kc = float(get_value(kc_ratio))
        
 
        print(f"From {ribs[i]} to {ribs[i+1]}")
        print(f"delta_a: {a}")
        print(f"chord/b: {b}")
        print(f"kc_ratio: {kc_ratio}")
        print(f"kc: {kc}")  

        if kc_ratio < 0.5: 
            sigma_cr = skin_buckling_wide(a, kc)
        else:
            sigma_cr = skin_buckling(b, kc) 
            
        safety_factor = safety_margin(normal_stresses[i], sigma_cr)

        ds = float(ribs[i])

        sigma_cr_tab.append(sigma_cr)

        safety_tab.append(safety_factor)
        ds_tab.append(ds)
    return sigma_cr_tab, safety_tab

#rib = [0, 0.5, 1, 1.5, 3.5, 4, 5, 6, 7.03893, 8, 9, 12, 16.0816]
#normal_stresses = [257, 248, 239, 230, 193, 184, 166, 149, 132, 111, 90, 33]

#ds_tab, sigma_cr_tab, safety_tab = Skin_buckling(rib, normal_stresses, 22e-3)
#plt.scatter(ds_tab, sigma_cr_tab)
#plt.show()
#print(safety_tab)