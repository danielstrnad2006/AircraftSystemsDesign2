import numpy as np
import scipy.integrate as sp
import matplotlib.pyplot as plt
import internal_loads


internal_properties=internal_loads.halfWing
internal_properties.set_conditions(2.5, 103544, 120, 1.225, 100)
#internal_properties.get_internal_plot()
moment_distribution = internal_properties.internal_bending
#print(moment_distribution)



class BeamDeflection:

    def __init__(self, b, db):
        self.b = b
        self.db = db
        self.I_XX = 0
    
    def assignI_XX(self, I_xx):
        self.I_XX = I_xx

    def findI_XX(self, location):
        y = np.arange(0, self.b / 2, self.db)
        return np.interp(location, y, self.I_XX)
    
    def Mx(self, y):
        return moment_distribution(y)


    def EIxx(self, y):
        return 71e9 * self.findI_XX(y)

    # INTEGRALS
    def dv_dy(self, y):
        integrand = lambda s: self.Mx(s) / self.EIxx(s)
        result, _ = sp.quad(integrand, 0, y)
        return -result

    def v(self, y):
        result, _ = sp.quad(self.dv_dy, 0, y)
        return result

    # MAX DEFLECTION
    def max_deflection(self, N=200):
        y_span = np.linspace(0, 16.0816, N)
        v_span = np.array([self.v(y) for y in y_span])*10**12

        max_v = np.max(v_span)
        max_y = y_span[np.argmax(v_span)]

        return max_v, max_y

    # PLOTS
    def plot(self, N=200):
        y_span = np.linspace(0, 16.0816, N)
        v_span = np.array([self.v(y) for y in y_span])
        v_span *= 1e12

        plt.figure()
        plt.plot(y_span, v_span)
        plt.xlabel("y")
        plt.ylabel("v(y)")
        plt.title("Deflection Distribution")
        plt.grid(True)

        plt.show()


#beam = BeamDeflection()
#beam.plot()

#max_v, max_y = beam.max_deflection()
#print(f"Maximum deflection: {max_v:.6f} at y = {max_y:.3f}")