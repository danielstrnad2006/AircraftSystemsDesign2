import math
import numpy as np
from scipy.optimize import fsolve
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt

# ---------------------------
# Wing Sizing Class
# ---------------------------
class WingSizing:
    def __init__(self, S_w, b, c_root, c_tip, taper_ratio,
                 leading_sweep, quart_sweep, dihedral):
        # Wing planform parameter
        self.S_w = S_w
        self.b = b
        self.c_root = c_root
        self.c_tip = c_tip
        self.taper_ratio = taper_ratio
        self.leading_sweep = leading_sweep
        self.quart_sweep = quart_sweep
        self.trailing_sweep = None
        self.dihedral = dihedral

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
    def plot(self, chord_position=0):

        self.trailing_sweep = math.degrees(math.atan(math.tan(math.radians(self.quart_sweep)) + 3 *
                                                     (self.c_root / (2 * self.b)) * (self.taper_ratio - 1)))  # [deg]

        y = np.linspace(0, self.b/2, 200)
        LE = np.array([self.leading_edge(yi) for yi in y])
        TE = np.array([self.trailing_edge(yi) for yi in y])

        print(LE, TE, y)



        # Fill area first
        plt.fill_between(y, -LE, -TE, color='gray', alpha=0.3)

        # Plot lines on top
        plt.plot(y, -LE, label='Leading Edge', color='gray', linewidth=2)
        plt.plot(y, -TE, label='Trailing Edge', color='gray', linewidth=2)

        # Fixed axes
        plt.xlim(0, max(y)*1.1)
        plt.ylim(-max(TE)*1.1,0)

        # Generate vertical line
        plt.plot([chord_position]*2, [-self.leading_edge(chord_position), -self.trailing_edge(chord_position)], linewidth=3, color='red')
        plt.savefig('temp/planform.png', bbox_inches='tight')


# wing = WingSizing(S_w=149.9, b=32.1632, c_root=7.2855, c_tip=2.0372, taper_ratio=0.2796,
#                  leading_sweep=37.537, quart_sweep=34.4871, dihedral=5)
#
# wing.plot(chord_position=0)