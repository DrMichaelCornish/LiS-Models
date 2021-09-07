#
# Standard parameters for lithium-sulfur battery models
#
import pybamm


class Cornish2021Parameters:
    """
    Standard parameters for lithium-sulfur battery models
    """

    def __init__(self):

        # Physical constants
        self.R = pybamm.constants.R
        self.F = pybamm.constants.F
        self.N = pybamm.constants.R/pybamm.constants.k_b
        self.T_ref = pybamm.Parameter("Reference Temperature [K]")
        self.Ta = pybamm.Parameter("Ambient Temperature [K]")

        # Known Parameters for any model
        self.Ms = pybamm.Parameter("Molar mass of S8 [g.mol-1]")
        self.ns = pybamm.Parameter("Number of S atoms in S [atoms]")
        self.ns2 = pybamm.Parameter("Number of S atoms in S2 [atoms]")
        self.ns4 = pybamm.Parameter("Number of S atoms in S4 [atoms]")
        self.ns8 = pybamm.Parameter("Number of S atoms in S8 [atoms]")
        self.ne = pybamm.Parameter("Electron number per reaction [electrons]")
        self.n_cells = pybamm.Parameter(
            "Number of cells connected in series to make a battery"
        )
        self.voltage_low_cut = pybamm.Parameter("Lower voltage cut-off [V]")
        self.voltage_high_cut = pybamm.Parameter("Upper voltage cut-off [V]")
        self.V_initial = pybamm.Parameter("Initial Condition for Terminal Voltage [V]")
        self.Tc_initial = pybamm.Parameter("Initial Cell Temperature [K]")
        self.timescale = 1
        self.dimensional_current_with_time = pybamm.FunctionParameter(
            "Current function [A]", {"Time[s]": pybamm.t * self.timescale}
        )

        
        # Standard Parameters
        self.il0 = pybamm.Parameter("Exchange current density L [A.m-2]")
        self.im0 = pybamm.Parameter("Exchange current density M [A.m-2]")
        self.ih0 = pybamm.Parameter("Exchange current density H [A.m-2]")
        self.m_s = pybamm.Parameter("Mass of active sulfur per cell [g]")
        self.rho_s = pybamm.Parameter("Density of precipitated Sulfur [g.L-1]")
        self.EH0 = pybamm.Parameter("Standard Potential H [V]")
        self.EM0 = pybamm.Parameter("Standard Potential M [V]")
        self.EL0 = pybamm.Parameter("Standard Potential L [V]")
        self.v = pybamm.Parameter("Electrolyte volume per cell [L]")
        self.ar = pybamm.Parameter("Active reaction area per cell [m2]")
        self.k_p = pybamm.Parameter("Precipitation rate [s-1]")
        self.S_star = pybamm.Parameter("S saturation mass [g]")
        self.k_s_charge = pybamm.Parameter(
            "Shuttle rate coefficient during charge [s-1]"
        )
        self.k_s_discharge = pybamm.Parameter(
            "Shuttle rate coefficient during discharge [s-1]"
        )
        self.f_s = pybamm.Parameter("Loss rate due to shuttle [s-1]")
        self.c_h = pybamm.Parameter("Cell heat capacity [J.g-1.K-1]")
        self.A = pybamm.Parameter("Pre-Exponential factor in Arrhenius Equation [J.mol-1]")
        self.h = pybamm.Parameter("Cell heat transfer coefficient [W.K-1]")
        self.m_c = pybamm.Parameter("Cell mass [g]")

        
        # Standard Initial Conditions
        self.S8_initial = pybamm.Parameter("Initial Condition for S8 ion [g]")
        self.S4_initial = pybamm.Parameter("Initial Condition for S4 ion [g]")
        self.S2_initial = pybamm.Parameter("Initial Condition for S2 ion [g]")
        self.S_initial = pybamm.Parameter("Initial Condition for S ion [g]")
        self.Sp_initial = pybamm.Parameter(
            "Initial Condition for Precipitated Sulfur [g]"
        )
        
        
        