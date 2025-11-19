import math
import numpy as np
from pygments.styles.dracula import purple

class VnDiagram:
    def __init__(self, weight, label, color, g=9.81, rho=1.225, S=149.9,
                 CL_max=1.452, CL_max_takeoff=2.547, V_C=255.85, V_D=308):
        self.W = weight              # weight in lbs
        self.label = label
        self.color = color
        self.g = g
        self.rho = rho
        self.S = S
        self.CL_max = CL_max
        self.CL_max_takeoff = CL_max_takeoff
        self.V_C = V_C
        self.V_D = V_D
        self.calculate_speeds_and_limits()

    def calculate_speeds_and_limits(self):
        # Stall speed
        self.V_s = math.sqrt(2 * self.W * 0.453592 * self.g / (self.CL_max * self.S * self.rho))
        self.n_max = max(2.1 + 24000/(self.W + 10000),2.5)
        self.n_min = -1
        self.V_A = self.V_s * math.sqrt(self.n_max)
        self.V_F = 1.6 * self.V_s  # flap speed for MTOW
        self.V_s0 = math.sqrt(2 * self.W * 0.453592 * self.g / (self.CL_max_takeoff * self.S * self.rho))
        self.upper_lim = self.V_s0 * math.sqrt(2)

    def plot(self, ax):
        # Positive stall curve
        V_up = np.linspace(0, self.V_A, 300)
        n_up = (V_up / self.V_s)**2

        # Flap curve (only for MTOW)
        V_flap = np.linspace(0, self.upper_lim, 100)
        n_flap = (V_flap / self.V_s0)**2

        V_extended = [self.upper_lim, min(self.V_F, math.sqrt(2)*self.V_s)]
        n_extended = [2, 2]

        # Positive limit
        V_pos = [self.V_A, self.V_D]
        n_pos = [self.n_max, self.n_max]

        # Vertical dive
        V_drop = [self.V_D, self.V_D]
        n_drop = [self.n_max, 0]

        # Linear from VD to VC
        V_lin = np.linspace(self.V_D, self.V_C, 200)
        n_lin = (1/(self.V_D-self.V_C))*(V_lin-self.V_D)  # linear from n=0 at VD to n=0 at VC

        # Negative stall
        V_neg = np.linspace(0, self.V_s, 300)
        n_neg = self.n_min * (V_neg / self.V_s)**2

        # Negative limit
        V_neg_lim = [self.V_s, self.V_C]
        n_neg_lim = [self.n_min, self.n_min]

        # Plot curves

        ax.plot(V_up, n_up, color=self.color)
        ax.plot(V_pos, n_pos, color=self.color)
        ax.plot(V_drop, n_drop, color=self.color)
        ax.plot(V_lin, n_lin, color=self.color)
        ax.plot(V_neg, n_neg, color=self.color)
        ax.plot(V_neg_lim, n_neg_lim, color=self.color)

        # Only plot flap curve for MTOW
        if self.W == 228275.445:
            ax.plot(V_flap, n_flap, color=purple, label="Flaps Extended Configuration")
            ax.plot(V_extended, n_extended, color=purple)

        ax.set_title(self.label)
        ax.set_xlabel("V_EAS (m/s)")
        ax.set_ylabel("Load Factor n")
        ax.grid(True)
        ax.set_ylim(-1.2, 3)
        ax.legend()


