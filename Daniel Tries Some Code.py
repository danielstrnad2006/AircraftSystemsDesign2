import scipy as sp
import numpy as np  
from scipy import interpolate
from scipy import signal 
import matplotlib.pyplot as plt


# This code loads wing data exactly as it is described in the assignment appendix B
footer_lines = 1029 #number of lines to disregard at the end of the file
span = 32.1632 #m



y_pos_lst = np.genfromtxt('data/MainWing_a=0.00_v=10.00ms.txt', skip_header=40, skip_footer=footer_lines, usecols=0).tolist()

chord_lst = np.genfromtxt('data/MainWing_a=0.00_v=10.00ms.txt', skip_header=40, skip_footer=footer_lines, usecols=1).tolist()
print(chord_lst)
chord_intrpl =sp.interpolate.interp1d(y_pos_lst, chord_lst, kind='linear', fill_value="extrapolate")


Ai_0_lst = np.genfromtxt('data/MainWing_a=0.00_v=10.00ms.txt', skip_header=40, skip_footer=footer_lines, usecols=2).tolist()
Cl_0_lst = np.genfromtxt('data/MainWing_a=0.00_v=10.00ms.txt', skip_header=40, skip_footer=footer_lines, usecols=3).tolist()
lCd_0_lst = np.genfromtxt('data/MainWing_a=0.00_v=10.00ms.txt', skip_header=40, skip_footer=footer_lines, usecols=5).tolist()
Cm_0_lst = np.genfromtxt('data/MainWing_a=0.00_v=10.00ms.txt', skip_header=40, skip_footer=footer_lines, usecols=7).tolist()

Ai_10_lst = np.genfromtxt('data/MainWing_a=10.00_v=10.00ms.txt', skip_header=40, skip_footer=footer_lines, usecols=2).tolist()
Cl_10_lst = np.genfromtxt('data/MainWing_a=10.00_v=10.00ms.txt', skip_header=40, skip_footer=footer_lines, usecols=3).tolist()
lCd_10_lst = np.genfromtxt('data/MainWing_a=10.00_v=10.00ms.txt', skip_header=40, skip_footer=footer_lines, usecols=5).tolist()
Cm_10_lst = np.genfromtxt('data/MainWing_a=10.00_v=10.00ms.txt', skip_header=40, skip_footer=footer_lines, usecols=7).tolist()

Ai_0_intrpl =sp.interpolate.interp1d(y_pos_lst, Ai_0_lst, kind='cubic', fill_value="extrapolate")
Ai_grad_lst = [(Ai_10 - Ai_0)/10 for Ai_10, Ai_0 in zip(Ai_10_lst, Ai_0_lst)]
Ai_grad_intrpl = sp.interpolate.interp1d(y_pos_lst, Ai_grad_lst, kind='cubic', fill_value="extrapolate")

Cl_0_intrpl =sp.interpolate.interp1d(y_pos_lst, Cl_0_lst, kind='cubic', fill_value="extrapolate")
Cl_grad_lst = [(Cl_10 - Cl_0)/10 for Cl_10, Cl_0 in zip(Cl_10_lst, Cl_0_lst)]
Cl_grad_intrpl = sp.interpolate.interp1d(y_pos_lst, Cl_grad_lst, kind='cubic', fill_value="extrapolate")

lCd_0_intrpl =sp.interpolate.interp1d(y_pos_lst, lCd_0_lst, kind='cubic', fill_value="extrapolate")
lCd_grad_lst = [(lCd_10 - lCd_0)/10 for lCd_10, lCd_0 in zip(lCd_10_lst, lCd_0_lst)]
lCd_grad_intrpl = sp.interpolate.interp1d(y_pos_lst, lCd_grad_lst, kind='cubic', fill_value="extrapolate")

Cm_0_intrpl =sp.interpolate.interp1d(y_pos_lst, Cm_0_lst, kind='cubic', fill_value="extrapolate")
Cm_grad_lst = [(Cm_10 - Cm_0)/10 for Cm_10, Cm_0 in zip(Cm_10_lst, Cm_0_lst)]
Cm_grad_intrpl = sp.interpolate.interp1d(y_pos_lst, Cm_grad_lst, kind='cubic', fill_value="extrapolate")

params_intrpl = [chord_intrpl, Ai_0_intrpl, Ai_grad_intrpl, Cl_0_intrpl, Cl_grad_intrpl, lCd_0_intrpl, lCd_grad_intrpl, Cm_0_intrpl, Cm_grad_intrpl]


class HalfWing:
    def __init__(self, params, v_ref=None, aoa_ref=None, rho_ref=None, g_loading=1, fuel_percentage_ref=100):
        # store interpolators (intercept and gradient) passed in params
        self.chord    = params[0]
        self.Ai_0     = params[1]; self.Ai_grad = params[2]
        self.Cl_0     = params[3]; self.Cl_grad = params[4]
        self.lCd_0    = params[5]; self.lCd_grad = params[6]
        self.Cm_0     = params[7]; self.Cm_grad = params[8]

        self.S = 149.9 #m
        self.b = 32.1632 #m^2
        self.y_engine = 6 #m # to be determined
        self.m_engine_and_nacelle = 3989.45376 #kg
        self.velocity = v_ref #m s^-1
        self.aoa = aoa_ref #deg
        self.rho = rho_ref #kg m^-3

        self.johannes_fuel_constant = 0.05884365  #fraction of fuel cross-sectional area over chord^2
        self.mass = 6215.1408 #kg
        self.fuel_percentage = fuel_percentage_ref #percent
        self.kerosene_density = 800 #kg m^-3
        self.g_loading = g_loading
        self.g = -9.80665 #kg m s^-2

        self.Cl_0_total = float(0.192557) # to be continued, read from the xflr results directly
        self.CL_grad_total= float((0.962596-0.192557)/10)

        # Inertial forces:
        L = self.b / 2.0
        c_root = float(self.chord(0))
        c_tip = float(self.chord(L))
        # integral of chord(y)^2 over 0..L for a linear chord: I = L*(c_root^2 + c_root*c_tip + c_tip^2)/3
        self.massConsts = self.mass/2 / (L * (c_root**2 + c_root*c_tip + c_tip**2) / 3.0)
        self.wing_mu = lambda y: self.massConsts * (float(self.chord(y))**2)
        self.fuel_volume_distribution = lambda y: ((self.chord(y)**2 * self.johannes_fuel_constant) if y<(self.b/2)*0.65  else 0.0)
        
        

    def compute_internal_forces(self):
        self.fuel_mass_distribution = lambda y:(self.fuel_percentage/100)*self.fuel_volume_distribution(y)*self.kerosene_density
        self.wing_mass_distribution = lambda y: self.fuel_mass_distribution(y)+self.wing_mu(y)


        self.x_cp_ratio = lambda y: (1/4) - self.get_Cm(y)/self.get_Cl(y)
        self.x_cp_distance = lambda y: self.x_cp_ratio(y)*self.chord(y)

        self.x_centroid_distance = lambda y: 1 #meter
                                        

        self.Lift = lambda y: 0.5 * self.rho * self.velocity**2 * self.chord(y) * self.get_Cl(y) #up positive lift
        # quick check: integrate to get total lift on the half wing
        self.total_lift_half, _ = sp.integrate.quad(self.Lift, 0, self.b / 2)
        print("Total lift is: ",2*self.total_lift_half,"[N]")


        cont_normal_force = lambda y: self.Lift(y) + self.g*self.g_loading*(self.wing_mass_distribution(y))

        reaction_shear = -(sp.integrate.quad(cont_normal_force, 0, self.b/2)[0])
        self.internal_shear = lambda y: -sp.integrate.quad(cont_normal_force, 0, y)[0] + (self.g * self.g_loading * self.m_engine_and_nacelle if y < self.y_engine else 0) - reaction_shear
        #self.total_bendingMoment, _ = sp.integrate.quad(lambda y: self.Lift(y)*y, 0, self.b/2)
        #self.bendingMoment = lambda y: y*self.ShearForce(y)-self.total_bendingMoment

        #self.torsionMoment = lambda y: (self.x_cp_distance(y)-self.x_centroid_distance(y)) * self.Lift(y) #Nm/m
        #self.torsionMoment_total, _ = sp.integrate.quad(lambda y: self.torsionMoment(y), 0, self.b/2)
        #self.internalTorsion = lambda y: -self.torsionMoment+sp.integrate.quad(lambda y: self.torsionMoment(y), 0, y)





    def set_conditions(self, velocity, CL_des, rho):
        self.velocity = velocity #m s^-1

        self.aoa = (CL_des-self.Cl_0_total)/self.CL_grad_total #deg
        print("target angle of attack is: ", self.aoa)
        self.rho = rho #kg m^-3
        self.compute_internal_forces()




    def _eval(self, intercept, grad, y, aoa=0.0):
        # aoa: same units used when computing gradients (here degrees)
        return float(intercept(y)) + aoa * float(grad(y))

    # convenience getters
    def get_Ai(self, y):
        return self._eval(self.Ai_0, self.Ai_grad, y, self.aoa)

    def get_Cl(self, y):
        return self._eval(self.Cl_0, self.Cl_grad, y, self.aoa)

    def get_lCd(self, y):
        return self._eval(self.lCd_0, self.lCd_grad, y, self.aoa)

    def get_Cm(self, y):
        return self._eval(self.Cm_0, self.Cm_grad, y, self.aoa)
    


    def get_coefficient_plots(self):
        ax = np.linspace(0, self.b/2, 200)
        Cl_plot = [self.Cl_0(y_pos) for y_pos in ax]
        Cm_plot = [self.Cm_0(y_pos) for y_pos in ax]
        x_cp_ratio_plot = [self.x_cp_ratio(y_pos) for y_pos in ax]
        x_cp_plot = [self.x_cp_distance(y_pos) for y_pos in ax ]

        #plt.plot(ax, x_cp_plot, label = "centre of pressure position [m] from LE")
        plt.plot(ax, Cl_plot, label = "lift coefficient distribution at zero aoa")
        plt.plot(ax, Cm_plot, label="moment coefficient distribution at zero aoa")
        plt.plot(ax, x_cp_ratio_plot, label="position of c.p. as ratio of chord")
        plt.legend()
        plt.show()
    


    def get_internal_shear_plot(self):
        y = np.linspace(0, self.b/2, 400)
        dry_wing_plot = [self.g*self.g_loading*self.wing_mu(y_pos) for y_pos in y]
        fuel_plot = [self.g*self.g_loading*self.fuel_mass_distribution(y_pos) for y_pos in y]
        lift_plot = [self.Lift(y_pos) for y_pos in y]
        Vz_plot = [self.internal_shear(y_pos) for y_pos in y]

        fig, ax1 = plt.subplots()
        l1, = ax1.plot(y, dry_wing_plot, label="wing structure weight distribution [N]")
        l2, = ax1.plot(y, fuel_plot, label="fuel weight distribution [N]")
        l3, = ax1.plot(y, lift_plot, label="Lift distribution [N]")
        ax1.set_xlabel("y [m]")
        ax1.set_ylabel("force per unit span [N]")

        ax2 = ax1.twinx()
        l4, = ax2.plot(y, Vz_plot, color='k', label="internal shear distribution [N]")
        ax2.set_ylabel("internal shear [N]")
        handles = [l1, l2, l3, l4]
        labels = [h.get_label() for h in handles]
        ax1.legend(handles, labels, loc='lower right')

        plt.show()

    def update_fuel(self, percentage):
        self.fuel_percentage = percentage
        self.compute_internal_forces()

    def update_g_loading(self, g_loading=1.0):
        self.g_loading = g_loading
        self.compute_internal_forces()   


halfWing = HalfWing(params_intrpl)
halfWing.set_conditions(velocity=120, CL_des=0.54, rho=1.225)
halfWing.get_internal_shear_plot()
halfWing.update_fuel(percentage=10)
halfWing.get_internal_shear_plot()
halfWing.update_fuel(percentage=100)
halfWing.update_g_loading(g_loading=-1.5)
halfWing.get_internal_shear_plot()

