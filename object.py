import math

class span_pos(self):
    def __init__(self, span, chord, Ai, Cl, lCd, Cm):
        self.span = span # span-wise location
        self.chord = chord
        self.Ai = Ai # AOA, induced
        self.Cl = Cl # Cl
        self.lCd = lCd # induced cd
        self.Cm = Cm # quarter chord moment

    def set_velocity(self, velocity):
        self.velocity = velocity
        self.f_l = 0.5 * chord