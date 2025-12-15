import math
from planform import *

# Wing definition
wing = WingSizing(S_w=149.9, b=32.1632, c_root=7.2855, c_tip=2.0372, 
                  taper_ratio=0.2796, leading_sweep=37.537, quart_sweep=34.4871, dihedral=5)

# Rib positions along the span (meters from root)
rib_positions = [1.0, 3.0, 6.0]

# Spar height function
def spar_height(span_pos):
    return 0.65 * wing.chord(span_pos)

# k_s lookup table
k_s_values = [
    (1, 1.05, 15), (1.05, 1.1, 13.8), (1.1, 1.2, 13.2), (1.2, 1.3, 12.6),
    (1.3, 1.4, 12), (1.4, 1.5, 11.8), (1.5, 1.6, 11.5), (1.6, 1.7, 11.2),
    (1.7, 1.8, 10.9), (1.8, 1.9, 10.7), (1.9, 2.0, 10.5), (2.0, 2.1, 10.4),
    (2.1, 2.2, 10.2), (2.2, 2.3, 10.0), (2.3, 2.4, 9.8), (2.4, 2.5, 9.75),
    (2.5, 2.6, 9.7), (2.6, 3.1, 9.6), (3.1, 3.3, 9.55), (3.3, 3.6, 9.5),
    (3.6, 5.0, 9.4)
]

def get_value(aspect_ratio):
    for low, high, val in k_s_values:
        if low <= aspect_ratio < high:
            return val
    # If no match is found, raise an error
    raise ValueError(f"No k_s value for aspect ratio = {aspect_ratio:.3f}")

# Material properties
poisson = 0.33
E = 72.4e9      # Pa
t = 2.5e-3      # m, plate thickness

for i, y in enumerate(rib_positions[:-1]):
    a = rib_positions[i+1] - rib_positions[i]      # spacing to next rib
    b_spar = spar_height(y)                        # spar height at this rib
    aspect_ratio = a / b_spar                      # aspect ratio

    k_s_val = get_value(aspect_ratio)             # will raise error if missing
    tau_crit_val = math.pi*2 * k_s_val * E / (12 * (1 - poisson2)) * (t / b_spar)*2

    print(f"Rib at {y} m → b_spar={b_spar:.3f} m, a={a:.3f} m, AR={aspect_ratio:.3f}, k_s={k_s_val}, tau_crit={tau_crit_val:.2f} Pa")