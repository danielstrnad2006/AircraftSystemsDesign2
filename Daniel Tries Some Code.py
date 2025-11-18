import scipy as sp
import numpy as np  
from scipy import interpolate
import matplotlib.pyplot as plt

# This code loads wing data exactly as it is described in the assignment appendix B
footer_lines = 1042 #number of lines to disregard at the end of the file
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




class HalfWing:
    def __init__(self, params, v_ref=None, aoa_ref=None, rho_ref=None):
        # store interpolators (intercept and gradient) passed in params
        self.chord    = params[0]
        self.Ai_0     = params[1]; self.Ai_grad = params[2]
        self.Cl_0     = params[3]; self.Cl_grad = params[4]
        self.lCd_0    = params[5]; self.lCd_grad = params[6]
        self.Cm_0     = params[7]; self.Cm_grad = params[8]

        self.S = 149.9 #m
        self.b = 32.1632 #m^2
        self.velocity = v_ref
        self.aoa = aoa_ref
        self.rho = rho_ref

        # integrate using wrappers because interp1d may return numpy types
        area_half, _ = sp.integrate.quad(self.chord, 0, self.b / 2)

        self.totalCl_0, _ = sp.integrate.quad(self.Cl_0, 0, self.b / 2)
        self.totalCl_grad, _ = sp.integrate.quad(self.Cl_grad, 0, self.b / 2)


    def set_conditions(self, velocity, aoa, rho):
        self.velocity = velocity
        self.aoa = aoa
        self.rho = rho
        self.velocity = velocity
        self.aoa = aoa
        self.rho = rho

        # create a callable that evaluates local lift (per unit span) at spanwise location y
        # interpolators must be called with y (they are callables)
        # get_Cl(y) already accounts for aoa via the stored gradient/intercept
        self.Lift = lambda y: 0.5 * self.rho * self.velocity**2 * self.chord(y) * self.get_Cl(y)
        # example: integrate to get total lift on the half wing
        self.total_lift_half, _ = sp.integrate.quad(self.Lift, 0, self.b / 2)
        self.ShearForce = lambda y: self.Lift(y)*y-self.total_lift_half
        self.total_bendingMoment, _ = sp.integrate.quad(lambda y: self.Lift(y)*y, 0, self.b/2)
        self.bendingMoment = lambda y: y*self.ShearForce(y)-self.total_bendingMoment


        ax = np.linspace(0, self.b/2, 200)
        Mz = [self.bendingMoment(y_pos) for y_pos in ax ]
        plt.plot(ax, Mz)
        plt.show()

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

halfWing = HalfWing(params_intrpl)
halfWing.set_conditions(velocity=150, aoa=8, rho=1.225)
