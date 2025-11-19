import math
import numpy as np
from scipy.optimize import fsolve
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# ---------------------------
# Airfoil Class
# ---------------------------
class Component:
    def __init__(self, t, L, x_pos, y_pos):
        self.t = t
        self.L = L
        self.x_pos = x_pos
        self.y_pos = y_pos

    @property
    def area(self):
        return self.t * self.L

class Spar(Component): 
    def __init__(self, t, h, x_pos, y_pos):
        super().__init__(t, h, x_pos, y_pos)


class Skin(Component): 
    def __init__(self, t, l, angle, x_pos, y_pos):
        super().__init__(t, l, x_pos, y_pos)
        self.angle =angle


class CrossSection:
    def __init__(self, xc_spar1, xc_spar2, chord, t_spar1, t_spar2, t_skin_up, t_skin_down, filepath):
        self.filepath = filepath
        self.xc_spar1 = xc_spar1
        self.xc_spar2 = xc_spar2
        self.chord = chord
        self.t_spar1 = t_spar1
        self.t_spar2 = t_spar2
        self.t_skin_up = t_skin_up
        self.t_skin_down = t_skin_down
        self.x, self.y, self.x_upper, self.y_upper, self.x_lower, self.y_lower = self._import_airfoil(filepath)


    def _import_airfoil(self, filepath):
        x, y = [], []
        with open(filepath, 'r') as file:
            for line in file:
                parts = line.strip().split()
                if len(parts) == 2:
                    try:
                        x_val, y_val = map(float, parts)
                        x_val = x_val*self.chord
                        y_val = y_val*self.chord
                        x.append(x_val)
                        y.append(y_val)
                    except ValueError:
                        continue

        x = np.array(x)
        y = np.array(y)

        # --- Find leading edge (minimum x) ---
        le_index = np.argmin(x)

        # Split into upper and lower surfaces
        x_upper = x[:le_index + 1]
        y_upper = y[:le_index + 1]
        x_lower = x[le_index:]
        y_lower = y[le_index:]

        # Ensure both are ordered left-to-right (increasing x)
        if x_upper[0] > x_upper[-1]:
            x_upper = np.flip(x_upper)
            y_upper = np.flip(y_upper)
        if x_lower[0] > x_lower[-1]:
            x_lower = np.flip(x_lower)
            y_lower = np.flip(y_lower)

        return x, y, x_upper, y_upper, x_lower, y_lower
    
    def get_thickness_at(self, x_c):
        y_u = np.interp(x_c, self.x_upper, self.y_upper)
        y_l = np.interp(x_c, self.x_lower, self.y_lower)
        return y_u - y_l
    
    def get_coordinate_at(self, x_c, top):
        if top: 
            y = np.interp(x_c, self.x_upper, self.y_upper)
        else:
            y = np.interp(x_c, self.x_lower, self.y_lower)
        return y
    
    def find_centroid(self, components):
        total_area = 0.0
        x_sum = 0.0
        y_sum = 0.0

        for comp in components:
            A = comp.area
            total_area += A
            x_sum += A * comp.x_pos
            y_sum += A * comp.y_pos
            
            print(f"x: {comp.x_pos:0.2f}\ty: {comp.y_pos:0.2f}\t A*x: {x_sum:0.2f} \tA*y: {y_sum:0.2f}\t A: {A:0.2f}")

        x_centroid = x_sum / total_area
        y_centroid = y_sum / total_area

        return x_centroid, y_centroid

    def assembly_centroid_finder(self):

        # upper and lower y coordinates of the spars
        y_spar1_up = self.get_coordinate_at(self.xc_spar1, True)
        y_spar1_down = self.get_coordinate_at(self.xc_spar1, False)
        y_spar2_up = self.get_coordinate_at(self.xc_spar2, True)
        y_spar2_down = self.get_coordinate_at(self.xc_spar2, False)

        # Height of the spars
        h_spar1 = y_spar1_up - y_spar1_down
        h_spar2 = y_spar2_up - y_spar2_down

        # Length and angle of the upper skin
        l_skin_up = math.sqrt((self.xc_spar1 - self.xc_spar2)**2 + (y_spar1_up - y_spar2_up)**2)
        theta_skin_up = np.atan((y_spar2_up - y_spar1_up)/(self.xc_spar2 - self.xc_spar1))
        
        # Length and angle of the lower skin
        l_skin_down = math.sqrt((self.xc_spar1 - self.xc_spar2)**2 + (y_spar1_down - y_spar2_down)**2)
        theta_skin_down = np.atan((y_spar2_down - y_spar1_down)/(self.xc_spar2 - self.xc_spar1))

        # Centroid of the first spar
        y_bar_spar1 = h_spar1/2
        x_bar_spar1 = self.xc_spar1

        # Centroid of the second spar
        y_bar_spar2 = h_spar2/2
        x_bar_spar2 = self.xc_spar2

        # Centroid of the upper skin
        y_bar_skin_up = l_skin_up * np.sin(theta_skin_up)
        x_bar_skin_up = l_skin_up * np.cos(theta_skin_up)
        
        # Centroid of the lower skin
        y_bar_skin_down = l_skin_down * np.sin(theta_skin_down)
        x_bar_skin_down = l_skin_down * np.cos(theta_skin_down)


        spar1 = Spar(self.t_spar1, h_spar1, x_bar_spar1, y_bar_spar1)
        spar2 = Spar(self.t_spar2, h_spar2, x_bar_spar2, y_bar_spar2 )

        skin_up = Skin(self.t_skin_up, l_skin_up, theta_skin_up ,x_bar_skin_up, y_bar_skin_up)
        skin_down = Skin(self.t_skin_down, l_skin_down, theta_skin_down, x_bar_skin_down, y_bar_skin_down)

        components = spar1, spar2, skin_up, skin_down
        print(self.find_centroid(components))


        print(y_spar1_up)
        print(y_spar1_down)
        print(y_spar2_up)
        print(y_spar2_down)
        
        print(h_spar1)
        print(h_spar2)

        print(l_skin_up)
        print(np.rad2deg(theta_skin_up))
        print(l_skin_down)
        print(np.rad2deg(theta_skin_down))
        return 

    def plot(self):
        plt.xlabel("x/c []")
        plt.ylabel("y/c []")
        plt.title("Airfoil cross section []")
        plt.plot(self.x, self.y, 'k')  # plot the airfoil outline

        plt.plot()
        #rectangle = patches.Rectangle((0.1, 0.1), 0.5, 0.3, linewidth=2, edgecolor='r', facecolor='none', angle=45)

        plt.axis("equal")
        plt.show()
        return 