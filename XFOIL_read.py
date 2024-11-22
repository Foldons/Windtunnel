import numpy as np
import matplotlib.pyplot as plt

file_path = r'C:\Users\frede\Downloads\viscous.polar'

alpha = []
cl = []

with open(file_path, 'r') as file:
    data_started = False
    for line in file:
        if data_started:
            try:
                # Split the line into columns
                columns = line.split()
                # Extract alpha and CL values
                alpha.append(float(columns[0]))
                cl.append(float(columns[1]))
            except (IndexError, ValueError):
                continue
        if line.strip().startswith('alpha'):  # Detect the header line
            data_started = True

# Plotting
plt.plot(alpha, cl)
plt.xlabel('AoA')
plt.ylabel('CL')
plt.grid(True)
plt.legend()
plt.show()