import matplotlib.pyplot as plt
import scipy as sp
from scipy.integrate import quad
import numpy as np

file_path= r'C:\Users\frede\Downloads\raw_testG9.txt'
Line= 8
chord = 16
single = False
List=np.genfromtxt(file_path,delimiter='')
rho = List[Line, 7]
v=22


def pressure(Line):
    pos = []
    p_top = []
    p_bottom = []
    p_back = []
    bottom = []
    pos2 = []
    x = 32
    b = 0
    for i in range(7,x):
        p_top.append(List[Line,i])
        pos.append(b)
        b+=0.667
    b=0
    for i in range(x,57):
        pos2.append(b)
        p_bottom.append(List[Line,i])
        b+=0.667
    for i in range(57,105):
        p_back.append(List[Line,i])
        bottom.append(i)
    return(p_top,pos,p_bottom,p_back,bottom)



def area_between_curves(y1, y2, x_positions):
    if len(y1) != len(y2) or len(y1) != len(x_positions):
        raise ValueError("Input lists must have the same length.")

    y1_interp = np.interp
    y2_interp = np.interp

    def integrand(x):
        y1_val = y1_interp(x, x_positions, y1)
        y2_val = y2_interp(x, x_positions, y2)
        return y1_val - y2_val

    total_area, _ = quad(integrand, x_positions[0], x_positions[-1])
    return total_area


if single == False:
    AOA = []
    CN_List = []
    for i in range(3,35):
        p_top, pos, p_bottom, p_back, bottom = pressure(i)
        area = area_between_curves(p_bottom, p_top, pos)
        CN = area / (0.5 * rho * v ** 2 * chord)
        AOA.append(List[i,2])
        CN_List.append(CN)
        #print(CN)
    plt.plot(AOA,CN_List)

    plt.show()
if single == True:
    p_top, pos, p_bottom, p_back, bottom = pressure(Line)
    area = area_between_curves(p_bottom, p_top, pos)
    CN = area/(0.5* rho * v**2 * chord)
    print( "AOA = "+ str(List[Line, 2]))
    print("normal coefficient = " +str(CN))
    plt.subplot(2, 1, 1)
    plt.plot(bottom,p_back)
    plt.subplot(2,1,2)
    plt.plot(pos,p_top)
    plt.plot(pos2,p_bottom)
    plt.gca().invert_yaxis()
    plt.grid()
    plt.show()
