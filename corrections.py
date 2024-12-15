import matplotlib.pyplot as plt
import math

import read_windtunnel_edit as windtunnel

from dataclasses import dataclass

@dataclass
class Corrections:

    h = ??? # windtunnel height, m
    c = 0.16 # chord length, m
    V_freestream = ???
    M_ref = math.sqrt(1.40 * 287 * 288.15) # mach 1 speed, m/s
    M = V_freestream / M_ref



    def get_sigma(self):
        sigma = ((math.pi**2) / 48) * ((self.c / self.h)**2)
        return sigma

    def get_tau(self):
        tau = (1/4) * (self.c/self.h)
        return tau

    def get_Lambda(self):
        for x in x_pos:

    def correct_Cl(self):

    def correct_Cd(self, c_d):
        Lambda = self.get_Lambda()
        tau = self.get_tau()
        sigma = self.get_sigma()
        part_1 = ((3 - 0.6(self.M^**2)) / ((1 - (self.M**2))**(3/2))) * Lambda * sigma
        part_2 = (((2 - (self.M**2)) * (1 + 0.4*(self.M**2))) / (1 - (self.M**2))) * tau * c_d
        c_d_corrected = c_d * (1 - part_1 - part_2)
        return c_d_corrected

    def correct_Cm(self):


    def correct_alpha(self):



