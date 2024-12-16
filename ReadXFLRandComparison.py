import matplotlib.pyplot as plt
import math as ma

with open("3D Simulation Raw Data.txt", "r") as file:
    simulation = file.readlines()

with open("raw_testG93D.txt", "r") as file:
    experiment = file.readlines()

alphaSim = []
CLSim = []
CDSim = []
CMSim = []
CLCDSim = []

alphaExp = []
CLExp = []
CDExp = []
CLCDExp = []

Fx = []
Fy = []


for i in range(49):
    alphaSim.append(float(simulation[8+i].split()[0].strip(',')))
    CLSim.append(float(simulation[8+i].split()[2].strip(',')))
    CDSim.append(float(simulation[8+i].split()[5].strip(',')))
    CLCDSim.append((float(simulation[8+i].split()[2].strip(',')))/(float(simulation[8+i].split()[5].strip(','))))
    CMSim.append(float(simulation[8+i].split()[8].strip(',')))

for i in range(31):
    alphaExp.append(float(experiment[2+i].split()[2].strip(',')))
    Fx.append(float(experiment[2+i].split()[6].strip(',')))
    Fy.append(float(experiment[2+i].split()[7].strip(',')))

for i in range(31):
    CLExp.append(Fy[i] / (0.5 * 1.177 * 22**2 * 0.06616))
    CDExp.append(Fx[i] / (0.5 * 1.177 * 22**2 * 0.06616))
    CLCDExp.append(CLExp[i]/CDExp[i])

plt.rcParams["figure.figsize"] = (8,4.5)

plt.subplot(2, 2, 1)
plt.plot(alphaSim,CLSim,'o--', label='Simulated $\mathregular{C_L}$', markersize=3, linewidth=1, color = 'orange', markerfacecolor='orange', markeredgecolor='black')
#plt.plot(alphaExp,CLExp,'^-', label='Experimental $\mathregular{C_L}$', markersize=4, linewidth=1, color = 'cornflowerblue', markerfacecolor='cornflowerblue', markeredgecolor='black')
plt.legend(loc = 'lower right')
plt.grid()
plt.xlabel(r'$\mathregular{\alpha}$ $[^\circ]$')
plt.ylabel("$\mathregular{C_L}$ $[-]$")

plt.subplot(2,2,2)
plt.plot(alphaSim,CMSim,'o--', label='Simulated $\mathregular{C_M}$', markersize=3, linewidth=1, color = 'orange', markerfacecolor='orange', markeredgecolor='black')
plt.legend(loc = 'upper right')
plt.grid()
plt.xlabel(r'$\mathregular{\alpha}$ $[^\circ]$')
plt.ylabel("$\mathregular{C_M}$ $[-]$")

#plt.subplot(2,2,2)
#plt.plot(alphaSim,CLCDSim,'o--', label='Simulated $\mathregular{C_L/C_D}$', markersize=3, linewidth=1, color = 'orange', markerfacecolor='orange', markeredgecolor='black')
#plt.plot(alphaExp,CLCDExp,'^-', label='Experimental $\mathregular{C_L/C_D}$', markersize=4, linewidth=1, color = 'cornflowerblue', markerfacecolor='cornflowerblue', markeredgecolor='black')
#plt.legend(loc = 'lower right')
#plt.grid()
#plt.xlabel(r'$\mathregular{\alpha}$ $[^\circ]$')
#plt.ylabel("$\mathregular{C_L/C_D}$ $[-]$")

plt.subplot(2,2,3)
plt.plot(CDSim,CLSim,'o--', label='Simulated $\mathregular{C_L}$', markersize=3, linewidth=1, color = 'orange', markerfacecolor='orange', markeredgecolor='black')
#plt.plot(CDExp,CLExp,'^-', label='Experimental $\mathregular{C_L}$', markersize=4, linewidth=1, color = 'cornflowerblue', markerfacecolor='cornflowerblue', markeredgecolor='black')
plt.legend(loc = 'lower right')
plt.grid()
plt.xlabel("$\mathregular{C_D}$ $[-]$")
plt.ylabel("$\mathregular{C_L}$ $[-]$")

plt.subplot(2,2,4)
plt.plot(alphaSim,CDSim,'o--', label='Simulated $\mathregular{C_D}$', markersize=3, linewidth=1, color = 'orange', markerfacecolor='orange', markeredgecolor='black')
#plt.plot(alphaExp,CDExp,'^-', label='Experimental $\mathregular{C_D}$', markersize=4, linewidth=1, color = 'cornflowerblue', markerfacecolor='cornflowerblue', markeredgecolor='black')
plt.legend(loc = 'upper left')
plt.grid()
plt.xlabel(r'$\mathregular{\alpha}$ $[^\circ]$')
plt.ylabel("$\mathregular{C_D}$ $[-]$")
plt.show()