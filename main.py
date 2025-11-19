# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.


#global panel_size = 0.8994

import matplotlib.pyplot as plt
import matplotlib
from vn_diagram import VnDiagram
matplotlib.use('TkAgg')

# V-n Diagram USAGE
weights = {
    "MTOW": 228275.445,
    "OEW": 96578.91,
    "ZFW": 138378.56
}
colors = ['blue', 'red', 'green']

fig, ax = plt.subplots(1, 3, figsize=(15, 5))

for axis, (name, W), color in zip(ax, weights.items(), colors):
    diagram = VnDiagram(W, f"V–n Diagram Cruise – {name}", color)
    diagram.plot(axis)

plt.tight_layout()
plt.show()