# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 16:35:33 2020

@author: diegopintossi

RED model (J. Veerman thesis, chapter 8)

All units converted to SI (e.g., cm --> m)


CO-FLOW configuration

         SW -----x---->
         RW -----y---->
         
            |----L----|

  @ x = 0, y = 0
  @ x = L, y = L
  
"""

from scipy.integrate import simps
from model_mechanics.exergy_TCPC import *


def power_density(Uload, E, Rcell,
                  W, L, x):
    """
    Calculation of the power and power density.

    Parameters
    ----------
    Uload : external load voltage [V].
    E : electromotive force (matrix) [V].
    Rcell : internal cell resistance (matrix) [Ohm m-2].
    W : width of the active area [m].
    L : length of the active area [m].
    x : discretized length [m].

    Returns
    -------
    Power_tot : power [W].
    Pd_avg : power density [W m-2].

    """
    
    U = Uload
    
    # Power density (gross) [W m-2]
    Pd = (E*U - U**2)/Rcell
    # Average power density (gross) [W m-2] (double Simpson rule for integration)
    Pd_avg = simps(Pd,x) / (2 * L)
    Power_tot = Pd_avg * (2 * W * L)
    
    return Power_tot, Pd_avg


def efficiencies(cRW, cSW,
                 fRW, fSW,
                 Temp, Power):
    """
    

    Parameters
    ----------
    cRW : NaCl concentration in RW (matrix) [mol m-3].
    cSW : NaCl concentration in SW (matrix) [mol m-3].
    fRW : RW flow rate [m3 s-1].
    fSW : SW flow rate [m3 s-1].
    Power : (gross) power [W].

    Returns
    -------
    dG_in : Gibbs free energy per sencond available in thegradient at the
        inlet [W].
    dG_out : Gibbs free energy per sencond available in thegradient at the
        outlet [W].
    deltaC_RW : change in NaCl concentration between inlet and outlet in the
        RW [g L-1].
    deltaC_SW : change in NaCl concentration between inlet and outlet in the
        SW [g L-1].
    eta_energy : energy efficiency [%].
    eta_thermo : thermodynamic efficiency [%].
    eta_net_energy : net energy efficiency [%].

    """
    
    cSWout = cSW[-1]
    cRWout = cRW[-1]
    
    fMix = fRW + fSW
    
    cMix = (fSW * cSW + fRW * cRW) / fMix
    cMixout = (fSW * cSWout + fRW * cRWout) / fMix
    
    dG_in = exergy(Temp, cSW[0], cRW[0], cMix[0], fSW, fRW, fMix)[1]
    dG_out = exergy(Temp, cSWout, cRWout, cMixout, fSW, fRW, fMix)[1]
    
    eta_energy = 100 * Power / dG_in
    eta_thermo = 100 * Power / (dG_in - dG_out)
    
    return eta_energy, eta_thermo


def current(J, W, x):
    """
    Calculation of the current produced by the stack.

    Parameters
    ----------
    J : current density [A m-2].
    W : width of the active area [m].
    x : discretized length [m].

    Returns
    -------
    I : current [A].

    """
    
    return W * simps(J, x)  # It would be (W*L) * (1/L) * simps(J,x)


def E_average(E, L, x):
    """
    Calculation of the current produced by the stack.

    Parameters
    ----------
    E : emf [V].
    L : length of the active area [m].
    x : discretized length [m].

    Returns
    -------
    average emf [V].

    """
    
    return simps(E, x) / L
