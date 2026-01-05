import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

def Column_buckling(ribs_input_lst, shear_stresses_input_lst, stringer_areas_input):
    #input section
    ribs=  ribs_input_lst[:-1] #rib locations along span of beam [m] # assuming the end of the wing has a rib

    normal_stresses=shear_stresses_input_lst  #maximum normal stress along that section between the ribs MPA [MPa]

    #this is for each stringer, so if you want to change the stringer from 0 to 1st rib, change the first value 
    stringer_area=stringer_areas_input #area of each stringer in m^2
    stringer_ratios=[1] * len(stringer_areas_input) #ratio is L/B where L is the stringer part that is parallel to the x-axis and B is the stringer length 
    stringer_thicknesses=[0.003] * len(stringer_areas_input) #m 

    #if there is no stringer from a certain point, watch out, don't input zeros, just remove the last value from x_data and output




    b=32.1632/2
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


    
    def column_buckling_along(ribs,normal_stresses,stringer_area,stringer_ratios,stringer_thicknesses):
        output=[]
        x_data=[0]
        for i in range(len(ribs)):
            if i == len(ribs) - 1:
                length=b-ribs[-1]
            else:
                length=ribs[i+1]-ribs[i]
            boundary_condition = 4
            x_data.append(length+x_data[-1])
            I_xx=stringer_Ixx(stringer_area[i], stringer_ratios[i], stringer_thicknesses[i])
            value=column_buckling(length, boundary_condition, I_xx, stringer_area[i], E=72.4*10**9)*10**(-6)
            print(value)
            safety=value/np.abs(normal_stresses[i])
            if safety>20:safety=20
            output.append(safety) #returns buckling safety factor
        x_data.pop(-1)
        column_buckling_interp1d = sp.interpolate.interp1d(x_data,output,kind="previous",bounds_error=False,fill_value=(output[0], output[-1]))
        
        return column_buckling_interp1d
    
    return column_buckling_along(ribs, normal_stresses, stringer_area, stringer_ratios, stringer_thicknesses)




