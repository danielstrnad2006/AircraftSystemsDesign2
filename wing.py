import math
import numpy as np
from scipy.optimize import fsolve
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt



# ---------------------------
# Wing Sizing Class
# ---------------------------
class WingSizing:
    def __init__(self, S_w, b, X_LE, X_TE, c_root, c_tip, taper_ratio, MAC, leading_sweep, quart_sweep, trailing_sweep, dihedral, b1, b2):
        # Wing planform parameter
        self.S_w = S_w
        self.b = b
        self.X_LE = X_LE
        self.X_TE = X_TE
        self.c_root = c_root
        self.c_tip = c_tip
        self.taper_ratio = taper_ratio
        self.MAC = MAC
        self.leading_sweep = leading_sweep
        self.quart_sweep = quart_sweep
        self.trailing_sweep = trailing_sweep
        self.dihedral = dihedral
        self.b1 = b1
        self.b2 = b2

    # ---------------------------
    # Helper function
    # ---------------------------
    def leading_edge(self, y):
        return math.tan(math.radians(self.leading_sweep)) * y

    def trailing_edge(self, y):
        return math.tan(math.radians(self.trailing_sweep)) * y + self.c_root

    def chord(self, y):
        return self.trailing_edge(y) - self.leading_edge(y)

    # ---------------------------
    # Plotting method
    # ---------------------------
    def plot(self):
        y = np.linspace(0, self.b/2, 200)
        LE = np.array([self.leading_edge(yi) for yi in y])
        TE = np.array([self.trailing_edge(yi) for yi in y])

        # Fill area first
        plt.fill_between(y, -LE, -TE, color='gray', alpha=0.3)

        # Plot lines on top
        plt.plot(y, -LE, label='Leading Edge', color='gray', linewidth=2)
        plt.plot(y, -TE, label='Trailing Edge', color='gray', linewidth=2)

        plt.xlabel('y [m]')
        plt.ylabel('x [m]')
        plt.title('Wing Planform')
        plt.legend()

        # Fixed axes
        plt.axis('equal')
        plt.xlim(0, max(y)*1.1)
        plt.ylim(-max(TE)*1.1,0)

        plt.show()