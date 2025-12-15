import numpy as np



def stress_at_location(Moment, y_distance_of_stringer, Ixx):
    return Moment*y_distance_of_stringer/Ixx



def column_buckling(stringer_length, Boundry_condition, I, A, E):
    stress_crit=Boundry_condition*np.pi**2*E*I/stringer_length**2/A
    return stress_crit

#stringer center of area
#ratio is L/B where L is the stringer part that is parallel to the x-axis and B is the stringer length 
def stringer_Ixx(area, ratio, thickness):
    L=area/(thickness*(1+1/ratio))
    B=L/ratio
    x_bar=B**2*thickness/(2*area)
    I_xx=1/12*L*thickness**3+L*thickness*x_bar**2+1/12*thickness*B**3+thickness*B*(B-x_bar)**2
    return I_xx

#input here

 #pascal - N/m^2
def buckling_yes_or_not(bending_stress, stringer_length, Boundry_condition, I, A): # bending_stress in N/m^2, all units in SI
    critical_stress= Boundry_condition*np.pi**2*E*I/stringer_length**2/A
    if critical_stress>bending_stress:
        return True # no buckling occurs
    else: return False

    