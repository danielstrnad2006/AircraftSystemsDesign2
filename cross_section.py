import math
import numpy as np
from scipy.optimize import fsolve
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import os
from Big_Q import *

# ---------------------------
# Airfoil Class
# ---------------------------
class Component:
    def __init__(self, t, L, x_pos, y_pos, angle=0):
        self.t = t
        self.L = L
        self.x_pos = x_pos
        self.y_pos = y_pos
        self.angle = angle

    @property
    def area(self):
        return self.t * self.L
    
    @property
    def lower_left_corner(self):
        w = self.t
        h = self.L
        cx = self.x_pos
        cy = self.y_pos
        theta = self.angle

        # Convert centroid to lower-left corner
        x0 = cx - (w/2)*np.cos(theta) + (h/2)*np.sin(theta)
        y0 = cy - (w/2)*np.sin(theta) - (h/2)*np.cos(theta)

        return x0, y0
    
class Spar(Component): 
    def __init__(self, t, h, x_pos, y_pos):
        super().__init__(t, h, x_pos, y_pos)


class Skin(Component): 
    def __init__(self, t, l, angle, x_pos, y_pos):
        super().__init__(t, l, x_pos, y_pos)
        self.angle =angle

class Stiffener(Component):
    def __init__(self, area, x_pos, y_pos, end_pos):
        super().__init__(0, 0, x_pos, y_pos)
        self._area = area
        self.x_pos = x_pos
        self.y_pos = y_pos
        self._end_pos = end_pos

    @property
    def area(self):
        return self._area

    @property
    def end_pos(self):
        return self._end_pos



    
class CrossSection:
    def __init__(self, xc_spar1, xc_spar2, chord, b_cur, t_spar1, t_spar2, t_skin_up, t_skin_down, stiffeners, filepath, display_data=False, display_plot=False, save_plot=False):
        self.filepath = filepath
        self.xc_spar1 = xc_spar1
        self.xc_spar2 = xc_spar2
        self.chord = chord
        self.b_cur = b_cur
        self.x_spar1 = xc_spar1*chord
        self.x_spar2 = xc_spar2*chord
        self.x, self.y, self.x_upper, self.y_upper, self.x_lower, self.y_lower = self._import_airfoil(filepath)
        self.display_data = display_data
        self.display_plot = display_plot
        self.save_plot = save_plot
        
        self.assembly_centroid_x = 0
        self.assembly_centroid_y = 0
        self.t_spar1 = t_spar1
        self.t_spar2 = t_spar2
        self.stiffener = stiffeners

        # upper and lower y coordinates of the spars
        self.y_spar1_up = self.get_coordinate_at(self.x_spar1, True)
        self.y_spar1_down = self.get_coordinate_at(self.x_spar1, False)
        self.y_spar2_up = self.get_coordinate_at(self.x_spar2, True)
        self.y_spar2_down = self.get_coordinate_at(self.x_spar2, False)

        # Height of the spars
        self.h_spar1 = self.y_spar1_up - self.y_spar1_down
        self.h_spar2 = self.y_spar2_up - self.y_spar2_down

        self.t_skin_up = t_skin_up
        self.t_skin_down = t_skin_down

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

        print(f"Generated for cross section y={self.b_cur:0.2f}m")

        if self.display_data:
            print(f"Determine Centroid Location: ") 

        if self.save_plot or self.display_plot:
            plt.xlabel("x []")
            plt.ylabel("y []")
            plt.title(f"Airfoil cross section y={self.b_cur}m")
            plt.plot(self.x, self.y, 'k')  # plot the airfoil outline

        for comp in components:
            A = 0
            x_indv = 0
            y_indv = 0

            if not isinstance(comp, Stiffener) or (isinstance(comp, Stiffener) and self.b_cur < comp.end_pos):
                A = comp.area
                x_indv = A * comp.x_pos
                y_indv = A * comp.y_pos

            total_area += A
            x_sum += x_indv
            y_sum += y_indv 

            if isinstance(comp, Stiffener) and self.b_cur < comp.end_pos:
                if self.save_plot or self.display_plot:
                    plt.plot(comp.x_pos, comp.y_pos, 'ro', markersize=6)  # stiffener as dot

            if not isinstance(comp, Stiffener):
                rectangle = patches.Rectangle(comp.lower_left_corner, comp.t, comp.L, linewidth=0, edgecolor='gray',
                                              facecolor='gray', angle=np.rad2deg(comp.angle))
                # draw centroid as a dot
                if self.save_plot  or self.display_plot:
                    plt.plot(comp.x_pos, comp.y_pos, 'ko', markersize=4)  # black dot ("k"), size 4
                    plt.gca().add_patch(rectangle)  # <- adds to current plot
                    plt.axis("equal")

            if self.display_data:
                print(f"L: {comp.L:0.2f}mm\th: {comp.t:0.2f}mm\t x: {comp.x_pos:0.2f}\ty: {comp.y_pos:0.2f}\t A*x: {x_indv:0.2f} \tA*y: {y_indv:0.2f}\t A: {A:0.2f}")

        self.assembly_centroid_x = x_sum / total_area
        self.assembly_centroid_y = y_sum / total_area
        if self.save_plot or self.display_plot:
            plt.plot(self.assembly_centroid_x, self.assembly_centroid_y, 'o', markersize=10)

        #Add spanwise location to the plot
        if self.save_plot or self.display_plot:
            try:
                spanwise_img = plt.imread("temp/planform.png")
                scaled_img = OffsetImage(spanwise_img, zoom=0.2)
                ab = AnnotationBbox(scaled_img, (1, 0), xycoords='axes fraction', box_alignment=(1.1, -0.1))
                plt.gca().add_artist(ab)
            except FileNotFoundError:
                print("ERROR: Planform image not found")

            plt.axis("equal")

            if not os.path.exists('Plots') and self.save_plot:
                os.mkdir('Plots')
            plt.savefig(f'Plots/cross_section_{self.b_cur:0.2f}.png', dpi=300)

            if self.display_plot:
                plt.show()
            else:
                plt.close()
            return
    
    #Function to find I_xx and I_yy values for the whole assembly (cross section)
    def find_AMOI(self, components):

        #Assing value to variables so that they are defined
        I_xx_sum = 0.0
        I_yy_sum = 0.0

        if self.display_data:
            print(f"Determine Area Moment of Inertia Location: ") 

        #Go through the list of component including the spars, skins, and stiffeners
        for comp in components:            
            A = comp.area
        
            I_xx = 0
            I_yy = 0
            
            #No need to compute I_xx, I_yy of stiffeners around their own centoridal axis (Assumed as 0)
            #If the component is not a stiffener, compute the I_xx, I_yy of it around its own centrodial axis
            if not isinstance(comp, Stiffener):
                I_xx = (comp.t * comp.L**3 * (np.sin(comp.angle + np.pi/2)**2)) / 12
                I_yy = (comp.t * comp.L**3 * (np.cos(comp.angle + np.pi/2)**2)) / 12


            I_xx_parallel_axis = 0
            I_yy_parallel_axis = 0
            #Find the parallel axis term of the component
            if not isinstance(comp, Stiffener) or (isinstance(comp, Stiffener) and self.b_cur < comp.end_pos):
                I_xx_parallel_axis = A * (comp.y_pos - self.assembly_centroid_y)**2
                I_yy_parallel_axis = A * (comp.x_pos - self.assembly_centroid_x)**2   

            #Add the I_xx, I_yy about own centroidal axis of the component with the parallel axis term
            I_xx_sum = I_xx_sum + I_xx + I_xx_parallel_axis
            I_yy_sum = I_yy_sum + I_yy + I_yy_parallel_axis
         
            if self.display_data:
                print(f"x: {comp.x_pos:0.2f}\ty: {comp.y_pos:0.2f}\t Ixx: {I_xx:0.2f} \tIyy: {I_yy:0.2f}\t A*y^2: {I_xx_parallel_axis:0.2f}\t A*x^2: {I_yy_parallel_axis:0.2f}")

        print(f"Ixx_sum: {I_xx_sum:0.2f} \tIyy_sum: {I_yy_sum:0.2f}")
        return I_xx_sum, I_yy_sum
    
    def find_polar_inertia(self, components):
        A = (self.x_spar2 - self.x_spar1) * (self.h_spar2 + self.h_spar1) * (1/2)
        ds_over_t = 0

        for comp in components:
            if not isinstance(comp, Stiffener):
                ds_over_t += comp.L / comp.t
            else:
                continue

        J = (4 * A**2)/ds_over_t
        return J
    
    def assembly_centroid_finder(self):
        # Length and angle of the upper skin
        l_skin_up = math.sqrt((self.x_spar1 - self.x_spar2)**2 + (self.y_spar1_up - self.y_spar2_up)**2)
        theta_skin_up = np.atan((self.y_spar2_up - self.y_spar1_up)/(self.x_spar2 - self.x_spar1))

        # Length and angle of the lower skin
        l_skin_down = math.sqrt((self.x_spar1 - self.x_spar2)**2 + (self.y_spar1_down - self.y_spar2_down)**2)
        theta_skin_down = np.atan((self.y_spar2_down - self.y_spar1_down)/(self.x_spar2 - self.x_spar1))

        # Centroid of the first spar
        y_bar_spar1 = self.y_spar1_down + self.h_spar1/2

        # Centroid of the second spar
        y_bar_spar2 = self.y_spar2_down + self.h_spar2/2

        # Centroid of the upper skin
        y_bar_skin_up = l_skin_up/2 * np.sin(theta_skin_up) + self.y_spar1_up
        x_bar_skin_up = l_skin_up/2 * np.cos(theta_skin_up) + self.x_spar1 
        
        # Centroid of the lower skin
        y_bar_skin_down = l_skin_down/2 * np.sin(theta_skin_down) + self.y_spar1_down
        x_bar_skin_down = l_skin_down/2 * np.cos(theta_skin_down) + self.x_spar1


        spar1 = Spar(self.t_spar1, self.h_spar1, self.x_spar1, y_bar_spar1)
        spar2 = Spar(self.t_spar2, self.h_spar2, self.x_spar2, y_bar_spar2)

        skin_up = Skin(self.t_skin_up, l_skin_up, theta_skin_up + np.pi/2, x_bar_skin_up, y_bar_skin_up)
        skin_down = Skin(self.t_skin_down, l_skin_down, theta_skin_down + np.pi/2, x_bar_skin_down, y_bar_skin_down)

        components = [spar1, spar2, skin_up, skin_down]

        # now create Stiffener objects from your list of tuples
        stiffener_objects = []
        for xc_rel, side, area, end_pos in self.stiffener:  # self.stiffener is list of (xc, 'up/down', area, end position)
            x_stiffener = xc_rel * self.chord  # convert relative chord to absolute x

            # get y position on the airfoil
            if side == 'up':
                y_stiffener = l_skin_up*(x_stiffener/(self.x_spar2-self.x_spar1))*np.sin(theta_skin_up) + self.y_spar1_up - self.t_skin_up/2
            else:
                y_stiffener = l_skin_down*(x_stiffener/(self.x_spar2-self.x_spar1))*np.sin(theta_skin_down) + self.y_spar1_down + self.t_skin_down/2

            # create Stiffener object
            stiff = Stiffener(area=area, x_pos=x_stiffener, y_pos=y_stiffener, end_pos=end_pos)
            stiffener_objects.append(stiff)

        # add stiffeners to components
        components.extend(stiffener_objects)

        self.find_centroid(components)

        I_xx_sum, I_yy_sum = self.find_AMOI(components)
        J_p = self.find_polar_inertia(components)
        return self.assembly_centroid_x, self.assembly_centroid_y, I_xx_sum, I_yy_sum, J_p

