import numpy as np
import scipy.integrate as sp
import matplotlib.pyplot as plt

#PARAMETERS/FUNCTION DEF
def Mx(y):

    return 10000 * (1 - y)   # Random value

def EIxx(y):

    return 5e6 * (0.7 + 0.3*y)   # Random value

#FIRST INTEGRAND
def dv_dy(y):
    integrand = lambda s: Mx(s) / EIxx(s)
    result, _ = sp.quad(integrand, 0, y)
    return -result

#SECOND INTEGRAND
def v(y):
    result, _ = sp.quad(dv_dy, 0, y)
    return result

# Span (normalised 0–1)
y_span = np.linspace(0, 1, 200)

# Compute slope and deflection
dv_span = [dv_dy(y) for y in y_span]
v_span = [v(y) for y in y_span]

# Plot
plt.figure()
plt.plot(y_span, dv_span)
plt.xlabel("y")
plt.ylabel("dv/dy")
plt.title("Slope distribution")

plt.figure()
plt.plot(y_span, v_span)
plt.xlabel("y")
plt.ylabel("v(y)")
plt.title("Deflection distribution")

plt.show()
