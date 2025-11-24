#torsion deflection
import matplotlib.pyplot as plt
import scipy as sp
from internal_loads import HalfWing, centroid_position
import numpy as np


G=10**7

def wing_tip_twist(self, bound):
    sum=self.internal_torque(self.b/2)
    def integrand(y):
        return (self.torque_per_span(y)-sum) / (self.G * self.polar_moment_inertia(y))

    twist_continuous, err = sp.integrate.quad(integrand, 0, bound)

    # point engine torque
    if bound >= self.y_engine:
        x_pos_e, z_pos_e = centroid_position(self.y_engine)
        engine_pos = np.array([
            self.engine_x_pos - x_pos_e,
            0,
            self.engine_z_pos - z_pos_e
        ])
        engine_force = np.array([
            -self.engine_thrust,
            0,
            self.m_engine_and_nacelle * self.g
        ])

        T_engine = np.cross(engine_pos, engine_force)[1]
        twist_point = T_engine / (self.G * self.polar_moment_inertia(self.y_engine))
    else:
        twist_point = 0

    return twist_continuous + twist_point

def twist_plot(HalfWing, N=400):
    import matplotlib.pyplot as plt
    import numpy as np

    # span stations
    y_vals = np.linspace(0, HalfWing.b/2, N)
    twist_vals = []

    for yi in y_vals:
        twist_vals.append(HalfWing.wing_tip_twist(yi))

    plt.plot(y_vals, twist_vals)
    plt.xlabel("Spanwise position y [m]")
    plt.ylabel("Torsional twist θ(y) [rad]")
    plt.title("Torsional twist distribution along the span")
    plt.grid(True)
    plt.show()

twist_plot(HalfWing)