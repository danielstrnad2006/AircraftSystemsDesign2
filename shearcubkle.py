
import math
poisson = 0.33
k_s = 1 ##  !!! dummy value ###
E = 72.4e9 # GPA
t = 2.5e-3    # m, design value
b =  2.0    # m, short side of the plate, height of the spar
a =   1.0     # m, distance between ribs
fraction = a/b
print(fraction)

k_s = float(input("input k_s = "))
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