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

#Code to check safety factors across length of beam
ribs= [0, 3, 6, 9, 12] #rib locations along span of beam
normal_stresses=[5*10**9,10**9,2*10**9,4*10**9,5*10**9] #max normal stress for each segment

#here you can change the parameters of the stringers along the section

stringer_area=[0.005,0.005,0.005,0.005,0.005] #area of each stringer in m^2
stringer_ratios=[1.2,1.2,1.2,1.3,1.3]  #ratio is L/B where L is the stringer part that is parallel to the x-axis and B is the stringer length 
stringer_thicknesses=[0.001,0.001,0.001,0.001,0.001] #m 

def column_buckling_along(ribs,normal_stresses,stringer_area,stringer_ratios,stringer_thicknesses):
    output=[]
    x_data=[0]
    for i in range(len(ribs)):
        if i == len(ribs) - 1:
            boundary_condition = 0.25
            length=b-ribs[-1]
        else:
            boundary_condition = 4
            length=ribs[i+1]-ribs[i]
        x_data.append(length+x_data[-1])
        I_xx=stringer_Ixx(stringer_area[i], stringer_ratios[i], stringer_thicknesses[i])
        value=column_buckling(length, boundary_condition, I_xx, stringer_area[i], E=72.4*10**9)
        output.append(value/normal_stresses[i]) #returns buckling safety factor
        print(length)
    x_data.pop(-1)
    graph = sp.interpolate.interp1d(x_data,output,kind="previous",bounds_error=False,fill_value=(output[0], output[-1]))
    x_plot = np.linspace(0, b, 300)
    sf_plot = graph(x_plot)
    plt.plot(x_plot, sf_plot)
    plt.xlabel("Spanwise location y [m]")
    plt.ylabel("Column Buckling Margin Of Safety[-]")
    plt.grid(True)
    plt.show()


column_buckling_along(ribs,normal_stresses,stringer_area,stringer_ratios,stringer_thicknesses) #


