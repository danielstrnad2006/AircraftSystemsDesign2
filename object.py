import math
from typing import Any

import numpy as np
import scipy

from main import panel_size, _MISSING


class span():

    def __init__(self, y_pos, chord, Ai, Cl, lCd, Cm):
        self.y_pos = y_pos # span-wise location
        self.chord = scipy.interpolate.UnivariateSpline(y_pos, chord)
        self.Ai = scipy.interpolate.UnivariateSpline(y_pos, Ai) # AOA, induced
        self.Cl = scipy.interpolate.UnivariateSpline(y_pos, Cl) # Cl
        self.lCd = scipy.interpolate.UnivariateSpline(y_pos, lCd) # induced cd
        self.Cm = scipy.interpolate.UnivariateSpline(y_pos, Cm) # quarter chord moment

        self.S = 149.9 # m^2
        self.b = 32.1632 # m

        self.CL, _ = scipy.integrate.quad(self.Cl, 0, self.b/2)

    def set_velocity_and_rho(self, velocity, rho):
        self.velocity = velocity
        self.rho = rho
        wing_force = lambda y: 0.5 * velocity ** 2 * self.rho * self.Cl(y) * self.chord(y)
        self.lift = lambda y: wing_force(y) * np.cos(self.Ai(y))
        self.drag = lambda y: wing_force(y) * np.sin(self.Ai(y))

        self.total_lift, _ = scipy.integrate.quad(self.lift, 0, self.b/2)
        self.total_drag, _ = scipy.integrate.quad(self.drag, 0, self.b/2)

    def set_centroid(self, centroid):
        self.centroid = scipy.interpolate.UnivariateSpline(self.y_pos, centroid)

        self.moment = lambda y: 0.5 * self.velocity ** 2 * self.rho * self.Cl(y) * self.chord(y) ** 2

        self.total_moment, _ = scipy.integrate.quad(self.moment, 0, self.b/2)



