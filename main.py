# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import internal_loads
from Wing_twist_deflection import  twist_plot, wing_tip_twist
from Bending_Deflection import BeamDeflection

internal_properties=internal_loads.halfWing
internal_properties.set_conditions(2.5, 103544, 120, 1.225, 100)
internal_properties.torque_plot()
twist_plot(internal_properties)
internal_properties.get_internal_plot()
moment_distribution = internal_properties.internal_bending

#PLOT BEAM DEFLECTION
beam = BeamDeflection()
beam.plot()
