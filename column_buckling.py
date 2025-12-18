import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

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


#input section
ribs=  [0, 2, 4, 6, 7.03893, 9, 12] #rib locations along span of beam [m] # assuming the end of the wing has a rib

normal_stresses=[np.float64(-221.55436933426083), np.float64(-184.66710646372886), np.float64(-149.41670185138847), np.float64(-132.10259507082256), np.float64(-90.00072765002402), np.float64(-33.40312927820247), np.float64(3.5356739692271065e-13)] #max normal stress for each segment
#this is for each stringer, so if you want to change the stringer from 0 to 1st rib, change the first value 
stringer_area=[0.000181,0.000181,0.000181,0.000181,0.000181,0.000181,0.000181,0.000181] #area of each stringer in m^2
stringer_ratios=[1,1,1,1,1,1,1,1]  #ratio is L/B where L is the stringer part that is parallel to the x-axis and B is the stringer length 
stringer_thicknesses=[0.003,0.003,0.003,0.003,0.003,0.003,0.003,0.003] #m 

#if there is no stringer from a certain point, watch out, don't input zeros, just remove the last value from x_data and output

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
    graph = sp.interpolate.interp1d(x_data,output,kind="previous",bounds_error=False,fill_value=(output[0], output[-1]))
    x_plot = np.linspace(0, b, 300)
    sf_plot = graph(x_plot)
    plt.plot(x_plot, sf_plot)
    plt.xlabel("Spanwise location y [m]")
    plt.ylabel("Column Buckling Margin Of Safety[-]")
    plt.grid(True)
    plt.show()
    return output


letter= column_buckling_along(ribs,normal_stresses,stringer_area,stringer_ratios,stringer_thicknesses) #
print((letter))

