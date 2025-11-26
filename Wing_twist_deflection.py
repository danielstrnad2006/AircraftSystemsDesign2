#torsion deflection
import matplotlib.pyplot as plt
import scipy as sp
from internal_loads import HalfWing, centroid_position
import numpy as np



def wing_tip_twist(object, bound):
    sum=object.internal_torque(object.b/2)
    def integrand(y):
        return (object.torque_per_span(y)-sum) / (object.G * object.polar_moment_inertia(y))
    twist_continuous, err = sp.integrate.quad(integrand, 0, bound)
    # point engine torque
    if bound >= object.y_engine:
        x_pos_e, z_pos_e = centroid_position(object.y_engine)
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

def twist_plot(half_wing, N=400):
    import matplotlib.pyplot as plt
    import numpy as np

    # span stations
    y_vals = np.linspace(0, (half_wing.b)/2, N)
    twist_vals = []

    for yi in y_vals:
        twist_vals.append(wing_tip_twist(half_wing, yi))

    plt.plot(y_vals, twist_vals)
    plt.xlabel("Spanwise position y [m]")
    plt.ylabel("Torsional twist θ(y) [rad]")
    plt.title("Torsional twist distribution along the span")
    plt.grid(True)
    plt.show()

