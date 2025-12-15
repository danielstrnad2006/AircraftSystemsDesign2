import matplotlib.pyplot as plt
import numpy as np

items=[]
output=[]
for i in np.arange(0,10,0.1):
    items.append(i)
    output.append(2.7**i*0.005)
plt.plot(items, output)
plt.xlabel("flight time [h]")
plt.ylabel("pilot boredom [%]")
plt.title("NASA Study on Pilot Boredom")
plt.show()

import numpy as np
import matplotlib.pyplot as plt

# Generate years
years = np.arange(1950, 2026, 1)

# Create a realistic salary model:
# Start around 30k in 1950, reach ~150k by 2025
# Use smooth exponential-like growth but without exploding
base_salary_1950 = 30000
final_salary_2025 = 150000

# Smooth exponential trend
growth_rate = np.log(final_salary_2025 / base_salary_1950) / (2025 - 1950)
salaries = base_salary_1950 * np.exp(growth_rate * (years - 1950))

plt.figure(figsize=(8, 5))
plt.plot(years, salaries)
plt.xlabel("Year")
plt.ylabel("Average Pilot Salary (Inflation Adjusted, USD)")
plt.title("Pilot Salary Trend (1950–2025)")

plt.show()