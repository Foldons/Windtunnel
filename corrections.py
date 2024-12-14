import matplotlib.pyplot as plt
import math

import read_windtunnel_edit as windtunnel

from dataclasses import dataclass

@dataclass
class Corrections:

    h = ??? # windtunnel height, m
    c = ??? # chord length, m


    def get_sigma(self):
        sigma = ((math.pi**2) / 48) * ((self.c / self.h)**2)
        return sigma

    def get_tau(self):
        tau = (1/4) * (self.c/self.h)
        return tau

    def get_Lambda(self):
        for x in x_pos:


