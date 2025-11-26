import scipy as sp
import numpy as np  
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




def centroid_position(y):
    # this is a function that will output the position of the centriod
    return 2,0.2

class HalfWing:
    def __init__(self, params, v_ref=None, aoa_ref=None, rho_ref=None, g_loading=1, fuel_percentage_ref=100):
        # store interpolators (intercept and gradient) passed in params
        self.CL_des   = None
        self.chord    = params[0]
        self.Ai_0     = params[1]; self.Ai_grad = params[2]
        self.Cl_0     = params[3]; self.Cl_grad = params[4]
        self.lCd_0    = params[5]; self.lCd_grad = params[6]
        self.Cm_0     = params[7]; self.Cm_grad = params[8]

        self.S = 149.9 #m
        self.b = 32.1632 #m^2
        self.G=28*10**9 #change this to the G value of aluminium we are using
        
        self.y_engine = 7.03893 #m # to be determined
        self.m_engine_and_nacelle = 3989.45376 #kg
        self.engine_thrust=363800/2 #N check this value
        self.y_engine = 7.04 #m # to be determined
        self.engine_z_pos=-1.9939/2-0.14*self.chord(self.y_engine) #from centerline 
        self.engine_x_pos=0 #from centerline leading edge, so we assume the center of mass of the engine coincides with the leading edge

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
        
    def polar_moment_inertia(self, y): #this code still has to be written properly
        
        cs = 1   # get geometry at this span station
        integral =1# ( cs.topspar_length    / cs.topspar_t +
          #  cs.bottomspar_length / cs.bottomspar_t +
          #  cs.left_wall_length  / cs.left_wall_t +
           # cs.right_wall_length / cs.right_wall_t)
        Area= 1#cs.Area
        J = 4 * (Area ** 2) / integral
        return 0.01  #this still has to be changed so it outputs the correct polar moment of inertia at the correct position

    def compute_internal_forces(self):
        self.fuel_mass_distribution = lambda y:(self.fuel_percentage/100)*self.fuel_volume_distribution(y)*self.kerosene_density
        self.wing_mass_distribution = lambda y: self.fuel_mass_distribution(y)+self.wing_mu(y)


        self.x_cp_ratio = lambda y: (1/4) - self.get_Cm(y)/self.get_Cl(y)
        self.x_cp_distance = lambda y: self.x_cp_ratio(y)*self.chord(y)

        self.x_centroid_distance = lambda y: 1 #meter
                                        
        Ai_case = lambda y: self.get_Ai(y)
        self.Lift = lambda y: 0.5 * self.rho * self.velocity**2 * self.chord(y) * self.get_Cl(y, Ai=Ai_case(y))   #up positive lift
        self.Drag = lambda y: self.Lift(y)*np.sin(np.deg2rad(Ai_case(y)))
        self.aerodynamic_normal= lambda y: np.cos(np.rad2deg(self.aoa))*self.Lift(y) + np.sin(np.rad2deg(self.aoa))*self.Drag(y)
        # quick check: integrate to get total lift on the half wing
        self.total_lift_half, _ = sp.integrate.quad(self.Lift, 0, self.b / 2)
        print("Total lift is: ",2*self.total_lift_half,"[N]")


        cont_normal_force = lambda y: self.Lift(y) + self.g*self.g_loading*(self.wing_mass_distribution(y))

        self.reaction_shear = -(sp.integrate.quad(cont_normal_force, 0, self.b/2)[0])
        self.integral_of_normal_force = lambda y: self.integrate_halfspan(cont_normal_force)(y)

        self.internal_shear = lambda y: -self.integral_of_normal_force(y) + (self.g * self.g_loading * self.m_engine_and_nacelle if y < self.y_engine else 0) - self.reaction_shear 
        self.internal_shear = self.function_to_intrp1d(self.internal_shear)
        #self.internal_shear = lambda y: -sp.integrate.quad(cont_normal_force, 0, y)[0] + (self.g * self.g_loading * self.m_engine_and_nacelle if y < self.y_engine else 0) - reaction_shear
        self.reaction_bending = self.integrate_halfspan(self.internal_shear)(self.b/2)
        #self.reaction_bending = self.function_to_intrp1d(self.reaction_bending)
        self.internal_bending = lambda y: self.integrate_halfspan(self.internal_shear)(y) - self.reaction_bending
        self.internal_bending = self.function_to_intrp1d(self.internal_bending)
        #print(self.reaction_bending)

        #this stopped working so I commented it out



        #self.internal_bending = lambda y: -sp.integrate.quad(self.internal_shear, y, self.b/2)[0] + self.reaction_bending
        #self.torsionMoment = lambda y: (self.x_cp_distance(y)-self.x_centroid_distance(y)) * self.Lift(y) #Nm/m
        #self.torsionMoment_total, _ = sp.integrate.quad(lambda y: self.torsionMoment(y), 0, self.b/2)
        #self.internalTorsion = lambda y: -self.torsionMoment+sp.integrate.quad(lambda y: self.torsionMoment(y), 0, y)


    def torque_per_span(self,y):
        # different torques: torque due to lift, no torque due to wing weight, point torque due to engine thrust and weight, 
        #  # this is to be changed\
        #positive torque is pointed towards the negative y axis

        x_pos, z_pos = centroid_position(y)
        return (self.x_cp_distance(y) - x_pos) * self.Lift(y)
    def internal_torque(self, bound):
        torque_lift, err = sp.integrate.quad(lambda y: self.torque_per_span(y), 0, bound)
        if bound >= self.y_engine:
            x_pos_e, z_pos_e = centroid_position(self.y_engine)
            engine_pos = np.array([
                self.engine_x_pos - x_pos_e,
                0,
                self.engine_z_pos - z_pos_e
            ])
            engine_force = np.array([
                -self.engine_thrust,
                0,
                self.m_engine_and_nacelle * self.g
            ])
            engine_torque = np.cross(engine_pos, engine_force)[1]
        else:
            engine_torque = 0
        return torque_lift + engine_torque
    
    def torque_plot(self):
        y_vals = np.linspace(0, self.b/2, 1000)
        torques = []
        sum=self.internal_torque(self.b/2)
        for yi in y_vals:
            cum_torque= self.internal_torque(yi)
            torques.append(cum_torque)
        plt.plot(y_vals, torques-sum)
        plt.xlabel("Spanwise position y [m]")
        plt.ylabel("Cumulative torque [Nm]")
        plt.grid(True)
        plt.show()   
        


 




    def _eval(self, intercept, grad, y, aoa_eff=0.0):
        # aoa: same units used when computing gradients (here degrees)
        return float(intercept(y)) + aoa_eff * float(grad(y))

    # convenience getters
    def get_Ai(self, y):
        return self._eval(self.Ai_0, self.Ai_grad, y, aoa_eff=self.aoa)

    def get_Cl(self, y, Ai=0):
        return self._eval(self.Cl_0, self.Cl_grad, y, aoa_eff=self.aoa+Ai)

    def get_lCd(self, y, Ai=0):
        return self._eval(self.lCd_0, self.lCd_grad, y, aoa_eff=self.aoa+Ai)

    def get_Cm(self, y, Ai=0):
        return self._eval(self.Cm_0, self.Cm_grad, y, aoa_eff=self.aoa+Ai)
    


    def get_coefficient_plots(self):
        ax = np.linspace(0, self.b/2, 200)
        Ai_plot = [self.Ai_0(y_pos) for y_pos in ax]
        #Cl_plot = [self.Cl_0(y_pos) for y_pos in ax]
        #Cm_plot = [self.Cm_0(y_pos) for y_pos in ax]
        #x_cp_ratio_plot = [self.x_cp_ratio(y_pos) for y_pos in ax]
        x_cp_plot = [self.x_cp_distance(y_pos) for y_pos in ax ]

        #plt.plot(ax, x_cp_plot, label = "centre of pressure position [m] from LE")
        plt.plot(ax, Ai_plot, label = "induced angle of attack at zero aoa")
        #plt.plot(ax, Cl_plot, label = "lift coefficient distribution at zero aoa")
        #plt.plot(ax, Cm_plot, label="moment coefficient distribution at zero aoa")
        #plt.plot(ax, x_cp_ratio_plot, label="position of c.p. as ratio of chord")
        plt.legend()
        plt.show()
    


    def get_forces_plot(self):
        y = np.linspace(0, self.b/2, 120)
        dry_wing_plot = [self.g*self.g_loading*self.wing_mu(y_pos) for y_pos in y]
        fuel_plot = [self.g*self.g_loading*self.fuel_mass_distribution(y_pos) for y_pos in y]
        lift_plot = [self.Lift(y_pos) for y_pos in y]
        fig, ax1 = plt.subplots()
        l1, = ax1.plot(y, dry_wing_plot, label="dry wing structure weight distribution [N/m]")
        l2, = ax1.plot(y, fuel_plot, label="fuel weight distribution [N/m]")
        l3, = ax1.plot(y, lift_plot, label="Lift distribution [N/m]")
        ax1.set_xlabel("y [m]")
        ax1.set_ylabel("force per unit span [N/m]")
        handles = [l1, l2, l3]
        labels = [h.get_label() for h in handles]
        ax1.legend(handles, labels, loc='lower right')

        plt.show()

    def get_internal_plot(self):
        y = np.linspace(0, self.b/2, 120)
        Vz_plot = [self.internal_shear(y_pos) for y_pos in y]
        Mx_plot = [self.internal_bending(y_pos) for y_pos in y]
        fig, ax1 = plt.subplots()
        l1, = ax1.plot(y, Vz_plot, label="Internal Shear Force [N]")
        ax1.set_xlabel("y [m]")
        ax1.set_ylabel("Internal Force [N]")

        ax2 = ax1.twinx()
        l2, = ax2.plot(y, Mx_plot, color='k', label="internal bending moment distribution [Nm]")
        ax2.set_ylabel("Internal Moment [N]")
        handles = [l1, l2]
        labels = [h.get_label() for h in handles]
        ax1.legend(handles, labels, loc='lower right')

        plt.show()

    def update_fuel(self, percentage):
        self.fuel_percentage = percentage
        self.compute_internal_forces()

    def update_g_loading(self, g_loading=1.0):
        self.g_loading = g_loading
        self.compute_internal_forces()


    # def set_conditions(self, velocity, CL_des, rho):
    #     self.velocity = velocity #m s^-1
    #
    #     self.aoa = (CL_des-self.Cl_0_total)/self.CL_grad_total #deg
    #     print("target angle of attack is: ", self.aoa)
    #     self.rho = rho #kg m^-3
    #     self.compute_internal_forces()



    def set_conditions(self, load_factor, weight, v_EAS, rho, fuel_percentage):
        self.fuel_percentage = fuel_percentage
        self.velocity = v_EAS*np.sqrt(1.225/rho) #m s^-1, convert EAS to TAS

        self.CL_des = (-weight * self.g * load_factor * 2) / (rho * self.velocity**2 * self.S)
        print("target design lift coefficient is: ", self.CL_des)
        self.aoa = (self.CL_des-self.Cl_0_total)/self.CL_grad_total #deg
        print("target angle of attack is: ", self.aoa)
        self.rho = rho #kg m^-3
        self.compute_internal_forces()



    def integrate_halfspan(self, integrand):
        """
        Numerically integrates the given integrand function over the half-span of the wing (from y=0 to y=b/2)
        using the trapezoidal rule, and returns the interpolated value of the integral.
        """
        dy = 0.1
        length = self.b/2
        y_pos_lst = np.arange(0, length, dy)
        #print(y_pos_lst)
        integrand_lst = [integrand(y_pos) for y_pos in y_pos_lst]
        #print(integrand_lst)
        integral_lst = [0]
        for i in range(1, len(integrand_lst)):
            new_integrand = integrand_lst[:i]
            integral_lst.append(sum(new_integrand)*dy)
        integral_cont = sp.interpolate.interp1d(y_pos_lst, integral_lst, kind='cubic', fill_value="extrapolate")
        return integral_cont

    def function_to_intrp1d(self, complicated_function):
        dy = 0.1
        length = self.b/2
        y_pos_lst = np.arange(0, length, dy)
        complicated_function_lst = [complicated_function(y_pos) for y_pos in y_pos_lst]
        interpolated = sp.interpolate.interp1d(y_pos_lst, complicated_function_lst, kind='cubic', fill_value="extrapolate")
        return interpolated

halfWing = HalfWing(params_intrpl)



def findCritical():
    
    #set conditions for loadcase 1: MTOW=103 544 kilograms, 
    # ZFW = 62 767.459 kilograms
    # OEW = 43 807.4567 kilograms
    conds = [[2.5, 103544, 138.016, 1.225, 100],
             [2.5, 103544, 163, 1.225, 100],  
             [-1, 103544, 87.29, 1.225, 100], 
             [-1, 103544, 154, 1.225, 100],
             
             [2.5, 43807, 89.772, 1.225, 0],
             [2.5, 43807, 163, 1.225, 0],  
             [-1, 43807, 56.777, 1.225, 0], 
             [-1, 43807, 154, 1.225, 0],

             [2.5, 62767, 107.46, 1.225, 0],
             [2.5, 62767, 163, 1.225, 0],  
             [-1, 62767, 67.96, 1.225, 0], 
             [-1, 62767, 154, 1.225, 0]]
    moment_lst = []
    for cond in conds: 
        halfWing.set_conditions(load_factor=cond[0], weight=cond[1], v_EAS=cond[2], rho=cond[3], fuel_percentage=cond[4])
        print()
        print("At load case with load factor: ", cond[0], ", mass", cond [1], "kg, Equivalent Air Speed: ", cond[2], "m/s, and density", cond[3], "kg/m^3")
        print("The Reaction Shear Force is:", halfWing.reaction_shear)
        print("The Reaction Bending Moment is:", halfWing.reaction_bending)
        moment_lst.append(halfWing.reaction_bending)
        print()
    print(max(moment_lst))
#findCritical()
