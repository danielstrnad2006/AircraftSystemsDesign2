import math
import matplotlib.pyplot as plt
from planform import *





def Skin_buckling(ribs, normal_stresses):
    span = 32.1632
    poisson = 0.33
    E = 72.4e9
    ds = 0

    sigma_cr_tab = []
    safety_tab = []
    ds_tab = []
    t = 2.5e-3 #[m]

    wing = WingSizing(S_w=149.9, b=span, c_root=7.2855, c_tip=2.0372, taper_ratio=0.2796,
                    leading_sweep=37.537, quart_sweep=34.4871, dihedral=5)


    def skin_buckling(b, kc):
        sigma_cr = ((((math.pi**2)*kc*E)/(12*(1-(poisson**2))))*((t/b)**2))/10**6
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
            (0.7, 0.8, 15.1),
            (0.8, 0.9, 13.5),
            (0.9, 1.0, 11.9),
            (1.0, 1.1, 10.3),
            (1.1, 1.15, 10.4),
            (1.15, 1.2, 10.5),
            (1.2, 1.3, 9.85),
            (1.3, 1.4, ),
            (1.4, 1.5, ),
            (1.5, 1.6, ),
            (1.6, 1.7, ),
            (1.7, 1.8, ),
            (1.8, 1.9, ),
            (1.9, 2, ),
            (2, 2.1, ),       
            (2.1, 2.2, ),
            (2.2, 2.3, ),
            (2.3, 2.4, ),
            (2.4, 2.5, ),  
            (2.5, 2.6, ),                                
        ]

        def get_value(kc_ratio):
            for low, high, val in kc_values:
                if low <= kc_ratio < high:
                    return val
                
            if kc_ratio >= 1.7:
                return 7.5
            else:
                return 20
                        
        kc = float(get_value(kc_ratio))

        print(f"From {ribs[i]} to {ribs[i+1]}")
        print(f"delta_a: {a}")
        print(f"chord/b: {b}")
        print(f"kc_ratio: {kc_ratio}")
        print(f"kc: {kc}")    

        sigma_cr = skin_buckling(b, kc)
        safety_factor = safety_margin(normal_stresses, sigma_cr)

        ds = float(ribs[i])

        sigma_cr_tab.append(sigma_cr)

        safety_tab.append(safety_factor)
        ds_tab.append(ds)
    return sigma_cr_tab, safety_tab