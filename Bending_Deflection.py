import numpy as np
import scipy.integrate as sp
import matplotlib.pyplot as plt

class BeamDeflection:

    def __init__(self):
        pass

    def Mx(self, y):
        return 10000 * (1 - y)

    def EIxx(self, y):
        return 71e9 * (0.7 + 0.3 * y)

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
        y_span = np.linspace(0, 1, N)
        v_span = np.array([self.v(y) for y in y_span])

        max_v = np.max(v_span)
        max_y = y_span[np.argmax(v_span)]

        return max_v, max_y

    # PLOTS
    def plot(self, N=200):
        y_span = np.linspace(0, 1, N)
        v_span = [self.v(y) for y in y_span]

        plt.figure()
        plt.plot(y_span, v_span)
        plt.xlabel("y")
        plt.ylabel("v(y)")
        plt.title("Deflection Distribution")
        plt.grid(True)

        plt.show()


beam = BeamDeflection()
beam.plot()

max_v, max_y = beam.max_deflection()
print(f"Maximum deflection: {max_v:.6f} at y = {max_y:.3f}")