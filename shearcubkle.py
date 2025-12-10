
import math
poisson = 0.33
k_s = 1 ##  !!! dummy value ###
E = 72.4e9 # GPA
t = 2.5e-3    # m, design value
b =  2.0    # m, short side of the plate, height of the spar
a =   1.0     # m, distance between ribs
aspect_ratio = a/b

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
        if low <= a < high:
            return val

k_s = get_value(aspect_ratio)
print(aspect_ratio)
print(k_s)

tau_crit = math.pi**2*k_s*E/(12*(1-poisson**2))*(t/b)**2

print(tau_crit)


shear_force = 1000 # N, shear force
h_f =   0.2    # m, height of the front spar    !!!dummy value
t_f =  2.5e-3  # m, thickness of the front spar !!! design dependant
h_r =   0.2    # m, height of the rear spar     !!!dummy value
t_r =  2.5e-3  # m, thickness of the rear spar  !!! design dependant
k_v =  1.2     # !!!!dummy value!!! 

tau_ave_shear = shear_force / (h_f*t_f + h_r*t_r)
tau_max_shear = k_v * tau_ave_shear

Torque = 1000 #Nm, torque
Area   = 0.3  #m^2, enclosed area of the wingbox
q_shear_flow = Torque / (2*Area)
