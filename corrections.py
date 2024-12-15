import matplotlib.pyplot as plt
import math
import numpy as np
import os
from scipy.integrate import quad
from scipy.interpolate import interp1d

import read_windtunnel_edit_2 as windtunnel

from dataclasses import dataclass

import matplotlib.pyplot as plt
import scipy as sp
from scipy.integrate import quad
import numpy as np
import os
from scipy.interpolate import interp1d



file_path = os.path.join(os.getcwd(), 'raw_testG9.txt')
map_path = os.path.join(os.getcwd(), 'mapping.txt')  # r'C:\Users\frede\Downloads\raw_testG9.txt'
List = np.genfromtxt(file_path, delimiter='')
v = 22
chord = 0.16
with open(map_path, 'r') as file:
    locs = [float(line.strip()) for line in file]
airfoilUpperPerc = np.array(locs[:25]) / 100 * chord
airfoilLowerPerc = np.array(locs[25:49]) / 100 * chord
totalWakemm = np.array(locs[49:96]) / 1000
staticWakemm = np.array(locs[96:109]) / 1000


def ambient(Line, List):
    # rho=List[Line, 7]
    pt = List[Line, 104]
    p_bar = List[Line, 4]
    T = List[Line, 5] + 273.15
    pb = List[Line, 3]

    # rho = 28.96 * p_bar /(T* 8.314)
    rho = p_bar * 100 / (287 * T)

    S = 110.4
    T0 = 273.15
    viscosity = 1.716 * 10 ** (-5) * (T / T0) ** (3 / 2) * ((T0 + S) / (T + S))

    q_inf = 0.211804 + 1.928442 * pb + 1.879374 * 10 ** (-4) * pb ** 2
    ps_inf = p_bar * 100 - q_inf  # (pt-q_inf)/2p_bar*100 - q_inf

    U_inf = np.sqrt(2 * q_inf / rho)
    Rey = rho * U_inf * chord / viscosity

    return rho, viscosity, U_inf, Rey, ps_inf, p_bar


def pressure(Line, p_inf):
    pos = []
    p_top = []
    p_bottom = []
    p_backtotal = []
    p_backstatic = []
    pos2 = []
    x = 33
    b = 0
    for i in range(8, x):
        p_top.append(List[Line, i] + p_inf)
        pos.append(b)
        b += chord / (25 - 1)  # 0.667
    b = 0
    for i in range(x, 57):
        pos2.append(b)
        p_bottom.append(List[Line, i] + p_inf)
        b += chord / (24 - 1)

    backtotal = np.linspace(0, 1, len(range(57, 104)))
    backstatic = np.linspace(0, 1, len(range(105, 117)))
    for i in range(57, 104):
        p_backtotal.append(List[Line, i] + p_inf)
    for i in range(105, 117):
        p_backstatic.append(List[Line, i] + p_inf)

    return p_top, pos, p_bottom, pos2, p_backtotal, backtotal, p_backstatic, backstatic


def cn(upper, lower, upperPos, lowerPos):
    lower_interp = interp1d(lowerPos, lower, bounds_error=False, fill_value="extrapolate")
    upper_interp = interp1d(upperPos, upper, bounds_error=False, fill_value="extrapolate")
    cn, _ = quad(lambda x: (lower_interp(x) - upper_interp(x)) / chord, 0, chord)
    return cn


def cm(upper, lower, upperPos, lowerPos):
    lower_interp = interp1d(lowerPos, lower, bounds_error=False, fill_value="extrapolate")
    upper_interp = interp1d(upperPos, upper, bounds_error=False, fill_value="extrapolate")
    cm, _ = quad(lambda x: (upper_interp(x) - lower_interp(x)) * x / chord ** 2, 0, chord)
    return cm


def dragRake(Uinf, Pinf, rho, Uwake, Pwake, Upos, pPos):
    Uw_interpd = interp1d(Upos, Uwake, bounds_error=False, fill_value="extrapolate")
    Pw_interpd = interp1d(pPos, Pwake, bounds_error=False, fill_value="extrapolate")
    firstTerm, _ = quad(lambda y: (Uinf - Uw_interpd(y)) * Uw_interpd(y), Upos[0], Upos[-1])
    thirdTerm, _ = quad(lambda y: Pinf - Pw_interpd(y), Upos[0], Upos[-1])
    cd = (rho * firstTerm + thirdTerm) / (0.5 * rho * Uinf ** 2 * chord)
    return cd

def get_plots():
    Line = 18
    single = False
    if single == False:
        AOA = []
        CL_List = []
        CD_List = []
        CM_List = []
        for j in range(3, 22):

            rho, viscous, U_inf, Rey, p_inf, pATM = ambient(j, List)
            p_top, pos, p_bottom, pos2, p_backtotal, backtotal, p_backstatic, backstatic = pressure(j, p_inf)

            for i in range(len(pos)):
                p_top[i] = (p_top[i] - p_inf) / (0.5 * rho * U_inf ** 2)
            for i in range(len(pos2)):
                p_bottom[i] = (p_bottom[i] - p_inf) / (0.5 * rho * U_inf ** 2)
            for i in range(len(backtotal)):
                cpt = (p_backtotal[i] - p_inf) / (0.5 * rho * U_inf ** 2)
                if cpt > 1: cpt = 1
                p_backtotal[i] = cpt * 0.5 * rho * U_inf ** 2 + p_inf

            backstatic_interpd = interp1d(staticWakemm, p_backstatic, bounds_error=False, fill_value="extrapolate")
            velocity_wake = []
            velocity_pos = []
            for i, x in enumerate(totalWakemm):
                if x < staticWakemm[0] or x > staticWakemm[-1]: continue
                vw = (2 * (p_backtotal[i] - backstatic_interpd(x)) / rho) ** 0.5
                velocity_wake.append(vw)
                velocity_pos.append(x)

            CN = cn(p_top, p_bottom, airfoilUpperPerc, airfoilLowerPerc)
            CM = cm(p_top, p_bottom, airfoilUpperPerc, airfoilLowerPerc)
            CM25 = CM + 0.25 * CN
            xCP = -CM / CN
            cdRake = dragRake(U_inf, p_inf, rho, velocity_wake, p_backstatic, velocity_pos, staticWakemm)
            alpha = List[j, 2]
            ar = np.deg2rad(alpha)
            CL = CN * (np.cos(ar) + (np.sin(ar) ** 2) / np.cos(ar)) - cdRake * np.tan(ar)
            AOA.append(alpha)
            CL_List.append(CL)
            CD_List.append(cdRake)
            CM_List.append(CM25)


        return AOA, CL_List, CD_List, CM_List


@dataclass
class Corrections:

    h = 0.9 # windtunnel height, m WRONGWRONGWRONGWRONGWRONGWRONGWRONGWRONGWRONGWRONGWRONGWRONGWRONGWRONGWRONGWRONGWRONG
    c = 0.16 # chord length, m
    V_freestream = 20 # WRONGWRONGWRONGWRONGWRONGWRONGWRONGWRONGWRONGWRONGWRONGWRONGWRONGWRONGWRONGWRONGWRONG
    M_ref = math.sqrt(1.40 * 287 * 288.15) # mach 1 speed, m/s
    M = V_freestream / M_ref
    t_c = 0.104



    def get_sigma(self):
        sigma = ((math.pi**2) / 48) * ((self.c / self.h)**2)
        return sigma

    def get_tau(self):
        tau = (1/4) * (self.c/self.h)
        return tau

    def correct_Cl(self, C_l, C_d):
        Lambda = self.get_Lambda()
        tau = self.get_tau()
        sigma = self.get_sigma()
        part_1 = (sigma / (1 - (self.M**2)))
        part_2 = ((2 - (self.M**2)) / ((1 - (self.M**2))**(3/2)))  * Lambda * sigma
        part_3 = (((2 - (self.M**2)) * (1 + 0.4*(self.M**2))) / (1 - (self.M**2))) * tau * C_d
        C_l_corrected = C_l * (1 - part_1 - part_2 - part_3)
        return C_l_corrected

    def correct_Cd(self, C_d):
        Lambda = self.get_Lambda()
        tau = self.get_tau()
        sigma = self.get_sigma()
        part_1 = ((3 - 0.6 * (self.M**2)) / ((1 - (self.M**2))**(3/2))) * Lambda * sigma
        part_2 = (((2 - (self.M**2)) * (1 + 0.4*(self.M**2))) / (1 - (self.M**2))) * tau * C_d
        c_d_corrected = C_d * (1 - part_1 - part_2)
        return c_d_corrected

    def correct_Cm(self, C_l, C_d, C_m):
        tau = self.get_tau()
        sigma = self.get_sigma()
        Lambda = self.get_Lambda()
        part_1 = ((2 - (self.M**2)) / ((1 - (self.M**2))**(3/2)))  * Lambda * sigma
        part_2 = (((2 - (self.M**2)) * (1 + 0.4*(self.M**2))) / (1 - (self.M**2))) * tau * C_d
        part_3 = C_l * (sigma / (4 * (1 - (self.M**2))))
        C_m_corrected = C_m * (1 - part_1 - part_2) + part_3
        return C_m_corrected


    def correct_alpha(self, C_l, C_m, alpha):
        sigma = self.get_sigma()
        alpha_corrected = alpha + ((57.3 * sigma) / (2 * math.pi * math.sqrt(1 - (self.M**2)))) * (C_l + 4*C_m)
        return alpha_corrected

    def get_Lambda(self):
        path_upper = os.path.join(os.getcwd(), 'Airfoil Upper.txt')
        path_lower = os.path.join(os.getcwd(), 'Airfoil Lower.txt')
        points_upper = np.flip(np.genfromtxt(path_upper, delimiter=''), 0)
        points_lower = np.genfromtxt(path_lower, delimiter='')
        slopes_upper = np.empty(len(points_upper) - 1)
        slopes_lower = np.empty(len(points_lower) - 1)

        for i in range(0, len(slopes_upper)):
            dy = points_upper[i + 1][1] - points_upper[i][1]
            dx = points_upper[i + 1][0] - points_upper[i][0]
            slope = dy / dx
            slopes_upper[i] = slope

        for i in range(0, len(slopes_lower)):
            dy = points_lower[i + 1][1] - points_lower[i][1]
            dx = points_lower[i + 1][0] - points_lower[i][0]
            slope = dy / dx
            slopes_lower[i] = slope

        return 0.200

if __name__ == '__main__':
    AOA_List, CL_List, CD_List, CM_List = get_plots()

    print(AOA_List, len(CL_List), len(CD_List), len(CM_List))

    correct = Corrections()

    AOA_corrected = []
    CL_corrected = []
    CD_corrected = []
    CM_corrected = []

    for i in range(len(AOA_List)):
        AOA_corrected.append(correct.correct_alpha(CL_List[i], CM_List[i], AOA_List[i]))
        CL_corrected.append(correct.correct_Cl(CL_List[i], CD_List[i]))
        CD_corrected.append(correct.correct_Cd(CD_List[i]))
        CM_corrected.append(correct.correct_Cm(CL_List[i], CD_List[i], CM_List[i]))

    plt.plot(AOA_List, CL_List, color = "blue", label = "original data", linestyle="-", marker="o")
    plt.plot(AOA_corrected, CL_corrected, color = "red", label = "corrected data", linestyle="-", marker="^")
    plt.xlabel('Angle of Attack [degrees]', fontweight='medium')
    plt.ylabel('Lift Coefficient [-]', fontweight='medium')
    plt.legend(loc = 'lower right')
    plt.grid()
    plt.show()

    plt.plot(AOA_List, CD_List, color="blue", label = "original data", linestyle="-", marker="o")
    plt.plot(AOA_corrected, CD_corrected, color="red", label = "corrected data", linestyle="-", marker="^")
    plt.xlabel('Angle of Attack [degrees]', fontweight='medium')
    plt.ylabel('Drag Coefficient [-]', fontweight='medium')
    plt.legend(loc = 'upper left')
    plt.grid()
    plt.show()

    plt.plot(AOA_List, CM_List, color="blue", label = "original data", linestyle="-", marker="o")
    plt.plot(AOA_corrected, CM_corrected, color="red", label = "corrected data", linestyle="-", marker="^")
    plt.xlabel('Angle of Attack [degrees]', fontweight='medium')
    plt.ylabel('Moment Coefficient [-]', fontweight='medium')
    plt.legend(loc = 'lower left')
    plt.grid()
    plt.show()

    plt.plot(CD_List, CL_List, color="blue", label = "original data", linestyle="-", marker="o")
    plt.plot(CD_corrected, CL_corrected, color="red", label = "corrected data", linestyle="-", marker="^")
    plt.xlabel('Drag Coefficient [-]', fontweight='medium')
    plt.ylabel('Lift Coefficient [-]', fontweight='medium')
    plt.legend(loc = 'lower right')
    plt.grid()
    plt.show()