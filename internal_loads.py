import scipy as sp
import numpy as np  
import matplotlib.pyplot as plt


# This code loads wing data exactly as it is described in the assignment appendix B
footer_lines = 1029 #number of lines to disregard at the end of the file
span = 32.1632 #m



y_pos_lst = np.genfromtxt('data/MainWing_a=0.00_v=10.00ms.txt', skip_header=40, skip_footer=footer_lines, usecols=0).tolist()

chord_lst = np.genfromtxt('data/MainWing_a=0.00_v=10.00ms.txt', skip_header=40, skip_footer=footer_lines, usecols=1).tolist()
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
    def __init__(self, params, v_ref=None, aoa_ref=None, rho_ref=None, g_loading=1, fuel_percentage_ref=100, ribs_locations=None):
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

        self.ribs_locations = ribs_locations if ribs_locations is not None else np.linspace(0, self.b/2, 5)
        
        self.y_engine = 7.03893 #m # to be determined
        self.m_engine_and_nacelle = 3989.45376 #kg
        self.engine_TO_thrust= 185500 #N (https://stands.aero/blog/aircraft-engines/pw2000-overview-and-specifications/)
        self.engine_z_pos=-1.9939/2-0.14*self.chord(self.y_engine) #from centerline 
        self.engine_x_pos=0 #from centerline leading edge, so we assume the center of mass of the engine coincides with the leading edge
        self.throttle_percentage =100

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
        self.x_centroid_distance = 0

        

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


                                    
        Ai_case = lambda y: self.get_Ai(y)

        self.x_cp_ratio = lambda y: (1/4) - self.get_Cm(y, Ai=Ai_case(y))/self.get_Cl(y, Ai=Ai_case(y))
        self.x_cp_distance = lambda y: self.x_cp_ratio(y)*self.chord(y)

        self.Lift = lambda y: 0.5 * self.rho * self.velocity**2 * self.chord(y) * self.get_Cl(y)   #up positive lift
        self.Drag = lambda y: np.abs(self.Lift(y)*np.sin(np.deg2rad(Ai_case(y))))
        self.aerodynamic_normal = lambda y: np.cos(np.deg2rad(self.aoa))*self.Lift(y) + np.sin(np.deg2rad(self.aoa))*self.Drag(y)
        self.aerodynamic_normal = self.aerodynamic_normal
        # quick check: integrate to get total lift on the half wing
        #self.total_lift_half, _ = sp.integrate.quad(self.Lift, 0, self.b / 2)
        #print("Total lift is: ",2*self.total_lift_half,"[N]")


        cont_normal_force = lambda y: self.aerodynamic_normal(y) + self.g*self.g_loading*(self.wing_mass_distribution(y))

        self.reaction_shear = -(self.integrate_halfspan(cont_normal_force)(self.b/2))

        

        self.integral_of_normal_force = lambda y: self.integrate_halfspan(cont_normal_force)(y)

        self.thrust = self.engine_TO_thrust*(self.rho/1.225) #*(1-self.velocity/400) #assuming effective exhaust velocity of 400 m/s

        self.internal_shear = lambda y: -self.integral_of_normal_force(y) + (self.g * self.g_loading * self.m_engine_and_nacelle if y < self.y_engine else 0) - self.reaction_shear 
        self.internal_shear = self.function_ribs_discretization(self.internal_shear) 
        #self.internal_shear = lambda y: -sp.integrate.quad(cont_normal_force, 0, y)[0] + (self.g * self.g_loading * self.m_engine_and_nacelle if y < self.y_engine else 0) - reaction_shear
        self.reaction_bending = self.integrate_halfspan(self.internal_shear)(self.b/2)

        self.internal_bending = lambda y: self.integrate_halfspan(self.internal_shear)(y) + (((self.thrust*np.abs(self.engine_z_pos)*np.sin(np.deg2rad(34))) if self.g_loading <0 else 0) if y < self.y_engine else 0) - self.reaction_bending
        self.internal_bending = self.function_to_intrp1d(self.internal_bending)


        self.torsion_distribtuion =  lambda y:  (self.x_centroid_distance(y)-self.x_cp_distance(y))*self.aerodynamic_normal(y)
        

        #print(self.thrust)
        #print(self.thrust/self.engine_TO_thrust*100)
        
        self.reaction_torsion = self.integrate_halfspan(self.torsion_distribtuion)(self.b/2)
        self.internal_torsion_noT = lambda y: self.integrate_halfspan(self.torsion_distribtuion)(y) - (self.x_centroid_distance(self.y_engine)*self.g * self.g_loading * self.m_engine_and_nacelle if y < self.y_engine else 0) - self.reaction_torsion
        self.internal_torsion_fullT = lambda y: self.integrate_halfspan(self.torsion_distribtuion)(y) - (self.x_centroid_distance(self.y_engine)*self.g * self.g_loading * self.m_engine_and_nacelle if y < self.y_engine else 0) + ((self.throttle_percentage/100)*self.thrust*self.engine_z_pos if y < self.y_engine else 0)- self.reaction_torsion
        
        self.internal_torsion_noT = self.function_ribs_discretization(self.internal_torsion_noT)
        self.internal_torsion_fullT = self.function_ribs_discretization(self.internal_torsion_fullT)
        

        self.sigma = lambda y: self.internal_bending(y) * self.y_max(y) / self.Q_buckling(y)

        self.q_torsion_noT = lambda y: self.internal_torsion_noT(y) / (2 * self.G * self.J(y))
        self.q_torsion_fullT = lambda y: self.internal_torsion_fullT(y) / (2 * self.G * self.J(y))
        self.q_bending = lambda y: self.internal_bending(y) * self.Q_buckling(y) / (self.I_xx(y))
        


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
                -self.engine_TO_thrust,
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

        


 




    def _eval(self, intercept, grad, y, aoa_eff):
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
        ax = np.linspace(0, self.b/2, 100)
        Ai_plot = [self.get_Ai(y_pos) for y_pos in ax]
        Cl_plot = [self.get_Cl(y_pos) for y_pos in ax]
        Cm_plot = [self.get_Cm(y_pos) for y_pos in ax]
        x_cp_ratio_plot = [self.x_cp_ratio(y_pos) for y_pos in ax]
        x_cp_plot = [self.x_cp_distance(y_pos) for y_pos in ax ]

        #plt.plot(ax, x_cp_plot, label = "centre of pressure position [m] from LE")
        plt.plot(ax, Ai_plot, label = "induced angle of attack at zero aoa")
        plt.plot(ax, Cl_plot, label = "lift coefficient distribution at zero aoa")
        plt.plot(ax, Cm_plot, label="moment coefficient distribution at zero aoa")
        plt.plot(ax, x_cp_ratio_plot, label="position of c.p. as ratio of chord")
        plt.legend()
        plt.show()
    


    def get_forces_plot(self):
        y = np.linspace(0, self.b/2, 120)
        dry_wing_plot = [self.g*self.g_loading*self.wing_mu(y_pos) for y_pos in y]
        fuel_plot = [self.g*self.g_loading*self.fuel_mass_distribution(y_pos) for y_pos in y]
        aero_plot = [self.aerodynamic_normal(y_pos) for y_pos in y]
        lift_plot = [self.Lift(y_pos) for y_pos in y]
        drag_plot = [self.Drag(y_pos) for y_pos in y]
        fig, ax1 = plt.subplots()
        l1, = ax1.plot(y, dry_wing_plot, label="dry wing structure weight distribution [N/m]")
        l2, = ax1.plot(y, fuel_plot, label="fuel weight distribution [N/m]")
        l3, = ax1.plot(y, aero_plot, label="aerodynamic normal force distribution [N/m]")
        l4, = ax1.plot(y, lift_plot, label="lift force distribution [N/m]")
        l5, = ax1.plot(y, drag_plot, label="drag force distribution [N/m]")
        ax1.set_xlabel("y [m]")
        ax1.set_ylabel("force per unit span [N/m]")
        handles = [l1, l2, l3, l4, l5]
        labels = [h.get_label() for h in handles]
        ax1.legend(handles, labels, loc='lower right')

        plt.show()

    def get_internal_plot(self):
        y = np.linspace(0, self.b/2, 200)
        Vz_plot = [self.internal_shear(y_pos)/1000 for y_pos in y]
        Mx_plot = [self.internal_bending(y_pos)/1000000 for y_pos in y]
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.5, 3.0))
        
        # Top subplot: Shear Force
        ax1.plot(y, Vz_plot, linewidth=2.5, label="Internal Shear Force [kN]", color='#1f77b4')
        ax1.set_ylabel("Internal Shear Force [kN]", fontsize=10, fontweight='bold', color='black')
        ax1.tick_params(axis='y', labelsize=9)
        ax1.tick_params(axis='x', labelsize=9)
        ax1.grid(True, alpha=0.3, linestyle='--')
        ax1.axhline(y=0, color='k', linewidth=0.8)
        ax1.legend(loc='best', fontsize=9, framealpha=0.9)
        ax1.set_xlim(0, self.b/2)
        
        # Bottom subplot: Bending Moment
        ax2.plot(y, Mx_plot, linewidth=2.5, color='#ff7f0e', label="Internal Bending Moment [MNm]")
        ax2.set_xlabel("Spanwise Position y [m]", fontsize=10, fontweight='bold')
        ax2.set_ylabel("Internal Bending Moment [MNm]", fontsize=10, fontweight='bold', color='black')
        ax2.tick_params(axis='y', labelsize=9)
        ax2.tick_params(axis='x', labelsize=9)
        ax2.grid(True, alpha=0.3, linestyle='--')
        ax2.axhline(y=0, color='k', linewidth=0.8)
        ax2.legend(loc='best', fontsize=9, framealpha=0.9)
        ax2.set_xlim(0, self.b/2)
        
        plt.tight_layout()

        plt.show()

    def get_internal_torsion_plot(self):
        y = np.linspace(0, self.b/2, 120)
        My_plot_noT = [self.internal_torsion_noT(y_pos) for y_pos in y]
        My_plot_fullT = [self.internal_torsion_fullT(y_pos) for y_pos in y]
        fig, ax1 = plt.subplots(figsize=(10, 6))
        ax1.plot(y, My_plot_noT, linewidth=2.5, label="Internal Torsional Moment at zero throttle [Nm]", color='#1f77b4')
        ax1.plot(y, My_plot_fullT, linewidth=2.5, label="Internal Torsional Moment at full throttle [Nm]", color='#ff7f0e')
        ax1.set_xlabel("Spanwise Position y [m]", fontsize=12, fontweight='bold')
        ax1.set_ylabel("Internal Torsion Moment M_y [Nm]", fontsize=12, fontweight='bold')
        ax1.tick_params(labelsize=11)
        ax1.axhline(y=0, color='k', linewidth=0.8)
        ax1.legend(fontsize=11, loc='best', framealpha=0.9)
        ax1.grid(True, alpha=0.3, linestyle='--')
        plt.tight_layout()
        plt.show()

    def get_debugging_torsion_plot(self):
        y = np.linspace(0, self.b/2, 120)
        x_cp_plot = [self.x_cp_distance(y_pos) for y_pos in y]
        x_centroid_plot = [self.x_centroid_distance (y_pos) for y_pos in y]
        fig, ax1 = plt.subplots()
        l1, = ax1.plot(y, x_cp_plot, label="Distance of centre of pressure to the LE [m]")
        l2, = ax1.plot(y, x_centroid_plot, label="Distance of centroid to the LE [m]")
        ax1.set_xlabel("y [m]")
        ax1.set_ylabel("Distance from the LE [m]")

        ax2 = ax1.twinx()
        aero_plot = [self.aerodynamic_normal(y_pos) for y_pos in y]
        l3, = ax2.plot(y, aero_plot, color='k', label="aerodynamic force distribution [N/m]")
        ax2.set_ylabel("Force per unit span [N/m]")
        handles = [l1, l2, l3]
        labels = [h.get_label() for h in handles]
        ax1.legend(handles, labels, loc='lower right')

        plt.show()

        J_plot = [self.J(y_pos) for y_pos in y]
        fig, ax1 = plt.subplots()
        l1, = ax1.plot(y, J_plot, label="torsion constant [m^4]")
        ax1.set_xlabel("y [m]")
        ax1.set_ylabel("torsion constant[m^4]")
        ax1.legend()

        plt.show()


    def update_fuel(self, percentage):
        self.fuel_percentage = percentage
        self.compute_internal_forces()

    def update_g_loading(self, g_loading=1.0):
        self.g_loading = g_loading
        self.compute_internal_forces()





    def set_conditions(self, load_factor, weight, v_EAS, rho, fuel_percentage):
        self.fuel_percentage = fuel_percentage
        self.velocity = v_EAS*np.sqrt(1.225/rho) #m s^-1, convert EAS to TAS
        self.g_loading = load_factor

        self.CL_des = (-weight * self.g * self.g_loading * 2) / (rho * self.velocity**2 * self.S)
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
    
    def set_torsion_params (self, db, x_centroid_arr_mm, J_arr_mm4):
        y = np.arange(0, self.b / 2, db)
        x_centroid_arr = [el*1e-3 for el in x_centroid_arr_mm]
        self.x_centroid_distance = sp.interpolate.interp1d(y, x_centroid_arr, kind='cubic', fill_value="extrapolate")
        J_arr = [el*1e-12 for el in J_arr_mm4]
        self.J = sp.interpolate.interp1d(y, J_arr, kind='cubic', fill_value="extrapolate")

    def set_buckling_params(self, db, Q_arr_mm, I_xx_arr_mm, spar_thickness, ribs_locations=None):
        y = np.arange(0, self.b / 2, db)
        Q_arr = [el*1e-9 for el in Q_arr_mm]
        I_xx_arr = [el*1e-12 for el in I_xx_arr_mm]
        self.Q_buckling = sp.interpolate.interp1d(y, Q_arr, kind='linear', fill_value="extrapolate")
        self.I_xx = sp.interpolate.interp1d(y, I_xx_arr, kind='linear', fill_value="extrapolate")
        self.y_max = lambda y: 0.0637 * self.chord(y)## Change when known!
        self.A_m = lambda y: self.johannes_fuel_constant * (self.chord(y)**2)
        self.thickness_front_spar = spar_thickness[0]*1e-3 #m
        self.thickness_rear_spar = spar_thickness[1]*1e-3 #m
        self.ribs_locations = ribs_locations

    def function_ribs_discretization(self, function):
        avg_values = []
        for i in range(len(self.ribs_locations)-1):
            domain = np.linspace(self.ribs_locations[i], self.ribs_locations[i+1], 10)
            function_arr = [function(y_pos) for y_pos in domain]
            integral = sp.integrate.trapezoid(function_arr, x=domain)
            segment_length = self.ribs_locations[i+1]-self.ribs_locations[i]
            avg_value = integral/segment_length
            #store avg_value for this segment
            avg_values.append(avg_value)
        #create a piecewise constant function based on avg_values
        piecewise_function = sp.interpolate.interp1d(self.ribs_locations[:-1], avg_values, kind='next', fill_value="extrapolate")
        return piecewise_function
    
    def get_shear_at_sections(self):
        tau_front_spar_noT = []
        tau_rear_spar_noT = []
        tau_front_spar_fullT = []
        tau_rear_spar_fullT = []
        shear_values_sections = []
        for y_pos in self.ribs_locations:
            tau_front_spar_noT.append(self.q_bending(y_pos) + self.q_torsion_noT(y_pos))/self.thickness_front_spar
            tau_rear_spar_noT.append(self.q_bending(y_pos) + self.q_torsion_noT(y_pos))/self.thickness_rear_spar
            tau_front_spar_fullT.append(self.q_bending(y_pos) + self.q_torsion_fullT(y_pos))/self.thickness_front_spar
            tau_rear_spar_fullT.append(self.q_bending(y_pos) + self.q_torsion_fullT(y_pos))/self.thickness_rear_spar

        for i in range(len(self.ribs_locations)-1):
            shear_values_sections.append(max(
                np.abs(tau_front_spar_noT[i]),
                np.abs(tau_rear_spar_noT[i]),
                np.abs(tau_front_spar_fullT[i]),
                np.abs(tau_rear_spar_fullT[i]),
                np.abs(tau_front_spar_noT[i+1]),
                np.abs(tau_rear_spar_noT[i+1]),
                np.abs(tau_front_spar_fullT[i+1]),
                np.abs(tau_rear_spar_fullT[i+1])
            ))
        return shear_values_sections
    
    def get_normal_stress_at_sections(self):
        normal_stress_values_ribs = []
        normal_stress_values_sections = []
        for y_pos in self.ribs_locations:
            normal_stress_values_ribs.append(self.sigma(y_pos))
        for i in range(len(normal_stress_values_ribs)-1):
            normal_stress_values_sections.append(max(normal_stress_values_ribs[i], normal_stress_values_ribs[i+1]))
        return normal_stress_values_sections
            
        

halfWing = HalfWing(params_intrpl)



def goThroughAll():
    
    #set conditions for loadcase 1: 
    # MTOW=103 544 kilograms, 
    # ZFW = 62 767.459 kilograms
    # OEW = 43 807.4567 kilograms
            #[n [-], mass[kg], v_EAS [m/s], rho, fuel percentage (100 for MTOW and 0 for OEW and ZFW)]
    conds = [ ]
    
    

    for cond in conds: 
        halfWing.set_conditions(load_factor=cond[0], weight=cond[1], v_EAS=cond[2], rho=cond[3], fuel_percentage=cond[4])
        print()
        print("At load case with load factor: ", cond[0], ", mass", cond [1], "kg, Equivalent Air Speed: ", cond[2], "m/s, and density", cond[3], "kg/m^3")
        print("The Reaction Shear Force is:", halfWing.internal_shear(0))
        print("The Reaction Bending Moment is:", halfWing.internal_bending(0))
        print("The Reaction Torsion Moment at full throttle is:", halfWing.internal_torsion_fullT(0))
        print("The Reaction Torsion Moment at zero throttle is:", halfWing.internal_torsion_noT(0))
        #halfWing.get_coefficient_plots()
        #halfWing.get_forces_plot()
        print()
        print()
    print()
    print()
    print()