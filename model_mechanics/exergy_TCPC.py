# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 18:02:58 2019

@author: dpin

DPIN, 2019-06-21

Script to calculate the exergy of a RED system
Based on the Excel spreadsheet of CSIM

This version is limited to NaCl solutions

"""

import numpy as np
from model_mechanics.TCPC import Salt

# Salt parameters for the TCPC model
NaCl = Salt(1, 1, 1, 1, 3.0210, 36.1573, 0.8998)


# Function definition
def moles(C,f):
    # calculation of molar amounts [mol s-1] for the various species in solution
    MW_NaCl, MW_H2O = 58.44/1000, 18.01/1000
    dens = (0.0375*(C/1000) + 0.9987)*1000   # [kg m-3]
    # calculated by linear interpolation of CRC handbook values
    # the factor 1000 is needed because the linear correlation
    # is valid for dens[kg L-1] and C[mol m-3]
    nNa = C*f                         # [mol s-1]
    nCl = C*f                         # [mol s-1]
    nH2O = f*(dens-C*MW_NaCl)/MW_H2O  # [mol s-1]
    nTot = nNa + nCl + nH2O           # [mol s-1]
    return [nNa, nCl, nH2O, nTot]


def fractions(n):
    # calculation of molar fraction [-] for the various species in solution
    # the input n needs to be a list in the form [nNa, nCl, nH2O, nTot]
    xNa = n[0]/n[3]
    xCl = n[1]/n[3]
    xH2O = n[2]/n[3]
    return [xNa, xCl, xH2O]


def activity(C):
    # TCPC model of Ge et al
    return NaCl.TCPC(C/1000)


def entropy(C,f):
    # calculation of entropy rate [W K-1] for a NaCl solution of concentration C
    n = moles(C,f)
    x = fractions(n)
    y = activity(C)
    R = 8.314
    entropy = - R * n[3] * (
            x[0]*np.log(y*x[0]) + x[1]*np.log(y*x[1]) + x[2]*np.log(x[2])
    )
    return [n, x, y, entropy]


def exergy(T,cSW,cRW,cMix,fSW,fRW,fMix):
    # exergy rate [W] calculation
    SW = entropy(cSW,fSW)
    RW = entropy(cRW,fRW)
    Mix = entropy(cMix,fMix)
    S = [SW[3], RW[3], Mix[3]]
    exergy = T * (Mix[3] - SW[3] - RW[3])
    return [S, exergy]
