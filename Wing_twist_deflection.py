#torsion deflection
import matplotlib.pyplot as plt
import scipy as sp
from internal_loads import HalfWing, centroid_position
import numpy as np


def wing_tip_twist(object, bound, J_P_intrp, centroid_X_intrp):
    sum=object.internal_torque(object.b/2)

    def integrand(y):
        return (object.torque_per_span(y)-sum) / (object.G * J_P_intrp(y))
    
    twist_continuous, err = sp.integrate.quad(integrand, 0, bound)
    
    # point engine torque
    if bound >= object.y_engine:
        x_pos_e = centroid_X_intrp(object.y_engine)
        z_pos_e = 0 # FIX THIS if it needs fixing!! 
        engine_pos = np.array([
            object.engine_x_pos - x_pos_e,
            0,
            object.engine_z_pos - z_pos_e
        ])
        engine_force = np.array([
            -object.engine_thrust,
            0,
            object.m_engine_and_nacelle * object.g
        ])

        T_engine = np.cross(engine_pos, engine_force)[1]
        twist_point = T_engine / (object.G * object.polar_moment_inertia(object.y_engine))
    else:
        twist_point = 0

    return twist_continuous + twist_point


def twist_plot(half_wing, centroid_X_arr, J_P_arr, b, db, N=60):

    y = np.arange(0, b / 2, db)
    centroid_X_intrp = np.interp(y, centroid_X_arr, kind='cubic', fill_value="extrapolate")
    J_P_intrp = np.interp(y, J_P_arr, kind='cubic', fill_value="extrapolate")


    # span stations
    y_vals = np.linspace(0, (half_wing.b)/2, N)
    twist_vals = []

    for yi in y_vals:
        twist_vals.append(wing_tip_twist(half_wing, yi, J_P_intrp, centroid_X_intrp))

    plt.plot(y_vals, twist_vals)
    plt.xlabel("Spanwise position y [m]")
    plt.ylabel("Torsional twist θ(y) [rad]")
    plt.title("Torsional twist distribution along the span")
    plt.grid(True)
    plt.show()

