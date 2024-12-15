import matplotlib.pyplot as plt
import math
import numpy as np
import os

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
        for j in range(3, 35):

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
        plt.subplot(3, 1, 1)
        plt.plot(AOA, CL_List)
        plt.grid()
        plt.xlabel('alpha')
        plt.ylabel("cl")
        plt.subplot(3, 1, 2)
        plt.plot(AOA, CM_List)
        plt.grid()
        plt.xlabel('alpha')
        plt.ylabel("cm")
        plt.subplot(3, 1, 3)
        plt.plot(CD_List, CL_List)
        plt.xlabel('cd')
        plt.ylabel("cl")
        plt.grid()
        plt.show()

    if single == True:

        rho, viscous, U_inf, Rey, p_inf, pATM = ambient(Line, List)
        p_top, pos, p_bottom, pos2, p_backtotal, backtotal, p_backstatic, backstatic = pressure(Line, p_inf)

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
        alpha = List[Line, 2]
        ar = np.deg2rad(alpha)
        CL = CN * (np.cos(ar) + (np.sin(ar) ** 2) / np.cos(ar)) - cdRake * np.tan(ar)

        print(f"AOA = {alpha}")
        print(f"normal coefficient = {CN}")
        print(f"moment coefficient = {CM25}")
        print(f"lift coefficient = {CL}")
        print(f"center of pressure = {xCP}")
        print(f"wake drag coefficient = {cdRake}")
        print(f"Reynolds numver = {Rey:.1e}")

        velocity = True
        if velocity:
            # plt.scatter(velocity_pos, velocity_wake, color='red', s=10, zorder=5)
            plt.rcParams.update({
                'font.size': 14,  # Set global font size
                'font.weight': 'medium'  # Set global font weight to bold
            })
            plt.plot(np.array(velocity_pos) * 1000, velocity_wake, 'o-', label='Experimental data', markersize=4,
                     linewidth=1, color='black', markerfacecolor='orange', markeredgecolor='black')
            plt.xlabel('Total pressure probe position [mm]', fontweight='medium')
            plt.ylabel('Wake velocity [m/s]', fontweight='medium')
            plt.ylim(16, 23)
            plt.grid()
            plt.legend(loc='lower right')
        else:
            plt.plot(pos, p_top)
            plt.plot(pos2, p_bottom)
            plt.gca().invert_yaxis()
            plt.grid()
        plt.show()

        return AOA, CL_List, CD_List, CM_List


@dataclass
class Corrections:

    h = 1 # windtunnel height, m WRONGWRONGWRONGWRONGWRONGWRONGWRONGWRONGWRONGWRONGWRONGWRONGWRONGWRONGWRONGWRONGWRONG
    c = 0.16 # chord length, m
    V_freestream = 20 # WRONGWRONGWRONGWRONGWRONGWRONGWRONGWRONGWRONGWRONGWRONGWRONGWRONGWRONGWRONGWRONGWRONG
    M_ref = math.sqrt(1.40 * 287 * 288.15) # mach 1 speed, m/s
    M = V_freestream / M_ref



    def get_sigma(self):
        sigma = ((math.pi**2) / 48) * ((self.c / self.h)**2)
        return sigma

    def get_tau(self):
        tau = (1/4) * (self.c/self.h)
        return tau

    #def get_Lambda(self):
        #for x in x_pos:

    #def correct_Cl(self):

    def correct_Cd(self, c_d):
        Lambda = self.get_Lambda()
        tau = self.get_tau()
        sigma = self.get_sigma()
        part_1 = ((3 - 0.6(self.M**2)) / ((1 - (self.M**2))**(3/2))) * Lambda * sigma
        part_2 = (((2 - (self.M**2)) * (1 + 0.4*(self.M**2))) / (1 - (self.M**2))) * tau * c_d
        c_d_corrected = c_d * (1 - part_1 - part_2)
        return c_d_corrected

    #def correct_Cm(self):


    #def correct_alpha(self):





if __name__ == '__main__':
    AOA, CL_List, CD_List, CM_List = get_plots()

