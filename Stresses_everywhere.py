import numpy as np
import matplotlib.pyplot as plt

def mohrs_circle(sigma_x, sigma_y, tau_xy, plot=True):
  

    # Average normal stress (circle center)
    sigma_avg = 0.5 * (sigma_x + sigma_y)

    # Radius of Mohr circle
    R = np.sqrt(((sigma_x - sigma_y) / 2)**2 + tau_xy**2)

    # Principal stresses
    sigma_1 = sigma_avg + R
    sigma_2 = sigma_avg - R

    # Maximum shear stress
    tau_max = R

    # Principal angle (physical plane)
    theta_p = 0.5 * np.arctan2(2 * tau_xy, sigma_x - sigma_y)
    theta_p_deg = np.degrees(theta_p)

    results = {
        "sigma_1": sigma_1,
        "sigma_2": sigma_2,
        "tau_max": tau_max,
        "theta_p_deg": theta_p_deg
    }

    if plot:
        theta = np.linspace(0, 2*np.pi, 400)
        sigma_circle = sigma_avg + R * np.cos(theta)
        tau_circle = R * np.sin(theta)

        plt.figure()
        plt.plot(sigma_circle, tau_circle)
        plt.scatter([sigma_x, sigma_y], [tau_xy, -tau_xy])
        plt.axhline(0)
        plt.axvline(0)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.xlabel("Normal stress σ")
        plt.ylabel("Shear stress τ")
        plt.title("Mohr's Circle")
        plt.grid(True)
        plt.show()

    return results

# Example stresses (MPa)
sigma_x = 120e6    #Make an input
sigma_y = -30e6    #Make an input
tau_xy  = 25e6     #Make an input 

res = mohrs_circle(sigma_x, sigma_y, tau_xy)

print(f"Principal stress σ1 = {res['sigma_1']/1e6:.2f} MPa")
print(f"Principal stress σ2 = {res['sigma_2']/1e6:.2f} MPa")
print(f"Max shear stress τmax = {res['tau_max']/1e6:.2f} MPa")
print(f"Principal angle θp = {res['theta_p_deg']:.2f} deg")
