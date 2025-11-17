import math
import scipy

from main import panel_size, _MISSING

class MissingDataError(Exception):
    pass

class span():
    def __init__(self, y_pos, chord, Ai, Cl, lCd, Cm):
        self.y_pos = scipy.interpolate(y_pos) # span-wise location
        self.chord = scipy.interpolate(chord)
        self.Ai = scipy.interpolate(Ai) # AOA, induced
        self.Cl = scipy.interpolate(Cl) # Cl
        self.lCd = scipy.interpolate(lCd) # induced cd
        self.Cm = scipy.interpolate(Cm) # quarter chord moment

    def cl(self):


