import matplotlib.pyplot as plt
import math

AOA_2D=[-6.0, -4.0, -2.0, -0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 11.0, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 15.5, 16.0]
CL_2D=[-0.40133471911397883, -0.19117962328851557, -0.04723691990466346, 0.12858164159281982, 0.3542836949317292, 0.5162908307623671, 0.7139343778945355, 0.8701333853966505, 0.9922392166956446, 1.0313203976618917, 1.0320165595722868, 0.9981623803782504, 0.9725333514471046, 0.9395180847998562, 0.8998516025084831, 0.9769520505583902, 0.9293572431137997, 0.9074604341463027, 0.7880163007591591]
CD_2D=[0.028161658036154634, 0.01500470609000878, 0.014430675306175183, 0.013284018243395786, 0.015404227852663455, 0.015947648341586237, 0.014952668193025966, 0.012994380809345465, 0.015056945929885601, 0.027821147597072188, 0.025874793233874704, 0.027140881170715864, 0.036729767435627725, 0.05019473011209663, 0.06673729697305666, 0.19677801941128592, 0.2239255681057234, 0.24059666149189285, 0.23095701805023655]
AOA_3D=[]
CL_3D=[]
CD_3D=[]
with open(r"raw_testG93D[1].txt") as file:
    for line in file:
        line=line.split()
        if line[0] in ["Run_nr","/"]:
            continue
        else:
            AOA_3D.append(float(line[2]))
            CL_3D.append(float(line[7])*2/(1.177*22**2*0.066146))
            CD_3D.append(float(line[6])*2/(1.177*22**2*0.066146)) 

#plot lift curve
#################################################################################

plt.plot(AOA_3D, CL_3D, marker='o', linestyle='-', color="red", label="3D wing")
plt.plot(AOA_2D, CL_2D, marker='^', linestyle='-', color="blue", label="2D wing") 
plt.xlabel("Angle of Attack (deg)")
plt.ylabel("Lift Coefficient (-)")
plt.title("Lift Curve")
plt.grid(True)
plt.legend()
ax = plt.gca()  # Get current axes
ax.spines['left'].set_position('zero')  # Y-axis at x=0
ax.spines['bottom'].set_position('zero')  # X-axis at y=0

# Keep ticks and data visible
ax.spines['left'].set_color('black')
ax.spines['bottom'].set_color('black')
ax.spines['right'].set_color('none')  # Hide right spine
ax.spines['top'].set_color('none')  # Hide top spine
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

# Show the plot
plt.show()

#plot drag polar 
##################################################################################

plt.plot(CD_3D, CL_3D, marker='o', linestyle='-', color="red", label="3D wing")
plt.plot(CD_2D, CL_2D, marker='^', linestyle='-', color="blue", label="2D wing") 
plt.xlabel("Drag Coefficient (-)")
plt.ylabel("Lift Coefficient (-)")
plt.title("Drag Polar")
plt.grid(True)
plt.legend()
ax = plt.gca()  # Get current axes
ax.spines['left'].set_position('zero')  # Y-axis at x=0
ax.spines['bottom'].set_position('zero')  # X-axis at y=0

# Keep ticks and data visible
ax.spines['left'].set_color('black')
ax.spines['bottom'].set_color('black')
ax.spines['right'].set_color('none')  # Hide right spine
ax.spines['top'].set_color('none')  # Hide top spine
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

# Show the plot
plt.show()
