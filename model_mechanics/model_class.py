# -*- coding: utf-8 -*-
"""
RED model (J. Veerman thesis, chapter 8)

All units converted to SI (e.g., cm --> m)

CO-FLOW configuration

         SW -----x---->
         RW -----y---->

            |----L----|
            
  @ x = 0, y = 0
  @ x = L, y = L

"""

from scipy.integrate import odeint
from model_mechanics.coFlow_NaCl_outputs import *


class CoFlow:
    def __init__(
            self,
            tRW, tSW,
            W, L, d,
            membranes,
            Uload,
            cRW0, cSW0,
            intervals,
    ):
        """Class to model RED in co-flow configuration.

        :param tRW: float
            Residence time of river water [s]
        :param tSW: float
            Residence time of seawater [s]
        :param W: float
            Width of the active area [m]
        :param L: float
            Length of the flow path [m]
        :param d: float
            Thickness of the water channels [m]
        :param membranes: string
            Type of membranes used
            ('ideal membranes' OR 'Fujifilm CEM-AEM type 10')
        :param Uload: float
            External load voltage per single cell pair [V]
        :param cRW0: float
            Inlet concentration of river water [mM]
        :param cSW0: float
            Inlet concentrations of seawater [mM]
        :param intervals: integer
            Number of intervals for the discretization of L [-]
        """
        # TCPC parameters
        self.NaCl = Salt(1, 1, 1, 1, 3.0210, 36.1573, 0.8998)

        # Constants
        self.R = 8.3143  # [J mol-1 K -1], universal gas constant
        self.Temp = 298.0  # [K], absolute temperature
        self.Fad = 96485.0  # [C mol-1], Farad constant
        self.Vw = 18E-6  # [m3 mol-1], molar volume of water
        self.obstr = 1.65  # [-], obstruction factor
        self.Lambda = 0.01287  # [m2] Ohm-1 mol-1], molar conductivity

        # Residence times
        self.tRW = tRW  # [s], residence time RW
        self.tSW = tSW  # [s], residence time SW

        # Stack parameters
        self.W = W  # [m], width of the active area
        self.L = L  # [m], length of the active area
        self.d = d  # [m], thickness of the water compartments

        # Flow rates [m3 s-1] CO-FLOW
        self.fSW = L * W * d / tSW  # flow rate of SW
        self.fRW = L * W * d / tRW  # flow rate of RW

        # Membrane data (const.)
        #  Raem = [Ohm m2], AEM area resistance
        #  Rcem = [Ohm m2], CEM area resistance
        #     a = [-], average permselectivity of the membranes
        # Dnacl = [m2 s-1], diffusion coefficient of NaCl in the membrane
        #    Dw = [m2 s-1], diffusion coefficient of water in the membrane
        #    Hm = [m], membrane thickness

        # Perfect membranes
        if membranes == 'Ideal CEM/AEM':
            self.Raem = 1.0E-4
            self.Rcem = 1.0E-4
            self.a = 0.90
            self.Dnacl = 0
            self.Dw = 0
            self.Hm = 80E-6
            self.Keosm = 0

        # Fujifilm CEM/AEM type 10
        if membranes == 'Fujifilm CEM/AEM type 10':
            self.Raem = 1.77E-4
            self.Rcem = 2.69E-4
            self.a = 0.946
            self.Dnacl = 1.5E-12
            self.Dw = 4.5E-9
            self.Hm = 145E-6
            self.Keosm = 6

        # Load
        self.Uload = Uload  # [V], electrode potential (to be adjusted for Pmax)

        # Initial/boundary conditions
        self.cSW0 = cSW0  # [mol m-3], initial concentration of NaCl in SW
        self.cRW0 = cRW0  # [mol m-3], initial concentration of NaCl in RW
        self.Qzero = [cSW0, cRW0]

        # Domain
        self.x = np.linspace(0, L, intervals)

        # Solved or not
        self.solved = False

        # Placeholders for the solution
        self.cRW, self.cSW = None, None
        self.E = None
        self.Rcell = None
        self.J = None
        self.power, self.pd_avg = None, None
        self.eta_energy, self.eta_thermo = None, None

    def water_ohmic_drop(self, C):
        """Method to calculate the resistance of a water channel.

        :param C: list
            Vector of concentrations [mol m-3 = mM]
        :return: list
            Vector of resistances [Ohm]
        """
        return self.obstr * self.d / (self.Lambda * C)

    def activity(self, C):
        """Method to calculate the activity along the water channel.
        Based on the TCPC model.

        :param C: list
            Vector of concentrations [mol m-3 = mM]
        :return: list
            Vector of activities [-]
        """
        return C * NaCl.TCPC(C / 1000)

    # RED Model
    # Definition of the system of ODEs
    def dQdx(self, Q, x):
        """
        Definition of the system of ODEs for scipy.integrate.odeint
        """
        cSW, cRW = Q

        # resistance of the river water compartment [Ohm m2]
        rRW = self.water_ohmic_drop(cRW)
        # resistance of the sea water compartment [Ohm m2]
        rSW = self.water_ohmic_drop(cSW)
        # activity of NaCl in seawater [-]
        aSW = self.activity(cSW)
        # activity of NaCl in seawater [-]
        aRW = self.activity(cRW)
        # Nernst equation defining the electromotive force [V]
        preF = self.a * 2 * self.R * self.Temp / self.Fad
        E = preF * np.log(aSW / aRW)
        # cell resistance [Ohm m2]
        Rcell = self.Raem + self.Rcem + rRW + rSW
        # unsegmented
        U = self.Uload  # [V], electrode potential
        J = (E - U) / Rcell  # [A m-2], current density

        # Salt flux [mol s-1 m-2]
        Tnacl = J / self.Fad + (cSW - cRW) * 2 * self.Dnacl / self.Hm
        # Water flux [m s-1]
        Tw = - (cSW - cRW) * 2 * self.Dw * self.Vw / self.Hm + self.Keosm * self.Vw * J / self.Fad

        # ODEs
        dcSWdx = - self.W * Tnacl / self.fSW + self.W * cSW * Tw / self.fSW
        dcRWdx = + self.W * Tnacl / self.fRW - self.W * cRW * Tw / self.fRW

        return [dcSWdx, dcRWdx]

    def solve_system(self):
        """
        Solve the system of ODEs using scipy.integrate.odeint
        """
        # Solution of the ODEs in the domain x
        solution = odeint(self.dQdx, self.Qzero, self.x)
        # Concentrations [mM]
        self.cRW = solution[:, 1]
        self.cSW = solution[:, 0]
        # emf [V]
        aRW = self.activity(self.cRW)
        aSW = self.activity(self.cSW)
        preF = self.a * 2 * self.R * self.Temp / self.Fad
        self.E = preF * np.log(aSW / aRW)
        # Cell resistance [Ohm m2]
        rRW = self.water_ohmic_drop(self.cRW)
        rSW = self.water_ohmic_drop(self.cSW)
        self.Rcell = self.Raem + self.Rcem + rRW + rSW
        # Current density [A m-2]
        self.J = (self.E - self.Uload) / self.Rcell
        # Power [W] and power density [W m-2] calculations
        self.power, self.pd_avg = power_density(
            self.Uload, self.E, self.Rcell, self.W, self.L, self.x
        )
        # Efficiency calculations [%]
        self.eta_energy, self.eta_thermo = efficiencies(
            self.cRW, self.cSW, self.fRW, self.fSW, self.Temp, self.power
        )
        # Set the system as solved
        self.solved = True



