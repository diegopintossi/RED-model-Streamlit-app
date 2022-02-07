# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 11:34:35 2019
@author: dpin

TCPC model, based on Ge et al. 2007 (DOI: 10.1021/je060451k)
Three characteristic parameters correlation

Empirical fitting based on experimental data, valid in wide range of ionic strengths
Similar behavior to Pitzer equations

DPIN, 2019-06-21
DPIN, 2019-07-11 (added osmotic coeff.)

"""

import numpy as np


class Salt:
    def __init__(self, zplus, zmin, vplus, vmin, b, S, n):
        self.zplus = zplus
        self.zmin = zmin
        self.vplus = vplus
        self.vmin = vmin
        self.b = b
        self.S = S
        self.n = n
    
    Aphi = 0.392    # valid at 298.15 K
    T = 298.15
    
    def TCPC(self, I): # calculate activity coefficient
        fgamma = - self.Aphi * (
            (I**0.5) / (1 + self.b*I**0.5) + (2/self.b) * np.log(1 + self.b*I**0.5)
        )
        gammaSV = (self.S/self.T)*I**(2*self.n)/(self.vplus + self.vmin)
        return np.exp((self.zplus * self.zmin) * fgamma + gammaSV)
    
    def osmo(self, I): # calculate osmotic coefficient
        return (
            1 - (self.zplus * self.zmin)
            * self.Aphi * ((I ** 0.5) / (1 + self.b * I ** 0.5))
            + (self.S / (self.T * (self.vplus + self.vmin)))
            * (2 * self.n / (2 * self.n + 1)) * I ** (2 * self.n)
        )
