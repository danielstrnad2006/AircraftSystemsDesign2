import numpy as np
import scipy.integrate as sp
import matplotlib.pyplot as plt



class TorsionDeflection:

    def __init__(self, b, db, torsion_distrib, J_distrib):
        self.b = b
        self.db = db
        self.I_XX = 0
        self.torsion_distribution = torsion_distrib
        self.J_distribution = J_distrib
        self.G = 28000000

    def dtheta_dy (self, y):
        integrand = self.torsion_distribution(y)/(self.G*self.J_distribution(y))
        return integrand
    
    def theta(self, y):
        result, _ = sp.quad(self.dtheta_dy, 0, y)
        return result

    # MAX DEFLECTION
    def max_rotation(self, N=200):
        y_span = np.linspace(0, 16.0816, N)
        theta_span = np.array([self.v(y) for y in y_span])

        max_theta = np.max(theta_span)
        max_y = y_span[np.argmax(theta_span)]

        return max_theta, max_y

    # PLOTS
    def plot(self, N=200):
        y_span = np.linspace(0, 16.0816, N)
        theta_span = np.array([self.theta(y)*180/np.pi for y in y_span])

        plt.figure()
        plt.plot(y_span, theta_span)
        plt.xlabel("y")
        plt.ylabel("theta(y) [deg]")
        plt.title("Deflection Distribution due to Torsion")
        plt.grid(True)

        plt.show()
