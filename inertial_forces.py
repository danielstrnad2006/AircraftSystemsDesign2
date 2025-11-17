import numpy as np



wing_weight=6215.1408
fuel_weight=0
engine_mass=0
loading=0
root_chord=7.2855
taper_ratio=0.2796
span=32.1632



def span_to_mass(span_location):
    const=wing_weight/(span*(1+taper_ratio)*root_chord/4)
    c=7.2855-0.31852*np.abs(span_location)
    inertial=const*c
    return inertial
sum=0
for i in np.arange(0,span/2,0.01):
    sum+=span_to_mass(i)*0.01
print(sum)