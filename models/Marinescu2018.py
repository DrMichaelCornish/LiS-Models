import pybamm
from .base_lithium_sulfur_model import BaseModel



class Marinescu2018(BaseModel):
    """
    Zero Dimensional model with the following electrochemistry
    
    S_{8}^{0} + 4e^{-} => 2 S_{4}^{2-}
    S_{4}^{2-} + 4e^{-} => S_{2}^{2-} + 2S^{2-}

    Parameters
    ----------
    options : dict, optional
        A dictionary of options to be passed to the model.
    name : str, optional
        The name of the model.

    References
    ----------
    .. [1]  Marinescu, M., O’Neill, L., Zhang, T., Walus, S., Wilson, T. E., &
            Offer, G. J. (2018).
            Irreversible vs reversible capacity fade of lithium-sulfur batteries
            during cycling: the effects of precipitation and shuttle. Journal of
            The Electrochemical Society, 165(1), A6107-A6118.
    """

    def __init__(self, options=None, name="Marinescu et al. (2018) Zero Dimensional Model"):
        super().__init__(options, name)
            
        # set external variables
        self.set_external_circuit_submodel()
        V = self.variables["Terminal voltage [V]"]
        I = self.variables["Current [A]"]

        # set internal variables
        S8 = pybamm.Variable("S8 [g]")
        S4 = pybamm.Variable("S4 [g]")
        S2 = pybamm.Variable("S2 [g]")
        S = pybamm.Variable("S [g]")
        Sp = pybamm.Variable("Precipitated Sulfur [g]")
        Ss = pybamm.Variable("Shuttled Sulfur [g]")
        

        #######################################
        # Model parameters as defined in table (1) in [1]. Parameters with 'H' or
        # 'L' in the name represent the high and low plateau parameter, respectively.
        #######################################
        from .parameters.marinescu2018_parameters import Marinescu2018Parameters
        self.param = Marinescu2018Parameters()
        param = self.param
        param = self.param

        # standard parameters
        R = param.R
        F = param.F
        Ta = param.Ta
        N = param.N
        T = param.T_ref    
        Tc = T

        # model-specific known parameters
        Ms = param.Ms
        ns = param.ns
        ns2 = param.ns2
        ns4 = param.ns4
        ns8 = param.ns8
        ne = param.ne
        ih0 = param.ih0
        il0 = param.il0
        rho_s = param.rho_s
        EH0 = param.EH0
        EL0 = param.EL0

        # model-specific unknown parameters
        v = param.v
        ar = param.ar
        k_p = param.k_p
        S_star = param.S_star
        k_s_charge = param.k_s_charge
        k_s_discharge = param.k_s_discharge
        f_s = param.f_s
        c_h = param.c_h
        m_c = param.m_c
        A = param.A
        h = param.h

        i_coef = ne * F / (2 * R * T)
        E_H_coef = R * T / (4 * F)
        f_h = (ns4 ** 2) * Ms * v / ns8
        f_l = (ns ** 2) * ns2 * Ms ** 2 * (v ** 2) / ns4

        #######################################################
        # Non-dynamic model functions
        #######################################################

        # High plateau potenital [V] as defined by equation (3a) in [1]
        E_H = EH0 + E_H_coef * ( pybamm.log(f_h) + pybamm.log(S8) - 2*pybamm.log(S4) )

        # Low plateau potenital [V] as defined by equation (3b) in [1]
        E_L = EL0 + E_H_coef * ( pybamm.log(f_l) + pybamm.log(S4) - pybamm.log(S2) - 2*pybamm.log(S) ) 

        # High plateau over-potenital [V] as defined by equation (7a) in [1]
        eta_H = V - E_H

        # Low plateau over-potenital [V] as defined by equation (7b) in [1]
        eta_L = V - E_L

        # High plateau current [A] as defined by equation (6a) in [1]
        i_H = -2 * ih0 * ar * pybamm.sinh(i_coef * eta_H)

        # Low plateau current [A] as defined by equation (6b) in [1]
        i_L = -2 * il0 * ar * pybamm.sinh(i_coef * eta_L)

        # Shuttle coefficient
        k_s =  k_s_discharge * (I >= 0) + k_s_charge *  (I < 0) * pybamm.exp(-A*N*((1/Tc)-(1/T))/R)

        ###################################
        # Dynamic model functions
        ###################################

        # Algebraic constraint on currents as defined by equation (5) in [1]
        algebraic_condition = i_H + i_L - I
        self.algebraic.update({V: algebraic_condition})

        # Differential equation (8a) in [2]
        dS8dt = -(ns8 * Ms * i_H / (ne * F)) - k_s * S8
        
        # Differential equation (8b) in [2] 
        dS4dt = (ns8 * Ms * i_H / (ne * F)) + (1- (f_s/Ms)*Ss)* k_s * S8 - (ns4 * Ms * i_L / (ne * F))
        
        # Differential equation (8c) in [2]
        dS2dt = ns2 * Ms * i_L / (ne * F)

        # Differential equation (8d) in [2]
        dSdt = (2 * ns * Ms * i_L / (ne * F)) - k_p * Sp * (S - S_star) / (v * rho_s)
    
        # Differential equation (8e) in [2]
        dSpdt = k_p * Sp * (S - S_star) / (v * rho_s)
        
        # Differential equation (8f) in [2]
        dSsdt = k_s * S8
        
        self.rhs.update({S8: dS8dt, S4: dS4dt, S2: dS2dt, S: dSdt, Sp: dSpdt,Ss : dSsdt})
        
        ##############################
        # Model variables
        #############################
        
        self.variables.update(
            {
                "Time [s]": pybamm.t * self.timescale,
                "Capacity [Ah]": pybamm.t * self.timescale * pybamm.AbsoluteValue(I) / 3600,
                "S8 [g]": S8,
                "S8 [g]": S8,
                "S4 [g]": S4,
                "S2 [g]": S2,
                "S [g]": S,
                "Precipitated Sulfur [g]": Sp,
                "Shuttle coefficient [s-1]": k_s,
                "Shuttle rate [g-1.s-1]": k_s * S8,
                "High plateau potential [V]": E_H,
                "Low plateau potential [V]": E_L,
                "High plateau over-potential [V]": eta_H,
                "Low plateau over-potential [V]": eta_L,
                "High plateau current [A]": i_H,
                "Low plateau current [A]": i_L,
                "Algebraic condition": algebraic_condition,
                "Shuttled Sulfur [g]" : Ss
            }
        )

        ######################################
        # Discharge initial condition
        # The values are found by considering the zero-current
        # state of the battery. Set S8, S4, and Sp as written
        # below. Then, solve eta_H = V, eta_L = V, the algebraic
        # condition, and mass conservation for the remaining values.
        ######################################

        self.initial_conditions.update(
            {
                self.variables["S8 [g]"]: param.S8_initial,
                self.variables["S4 [g]"]: param.S4_initial,
                self.variables["S2 [g]"]: param.S2_initial,
                self.variables["S [g]"]: param.S_initial,
                self.variables["Precipitated Sulfur [g]"]: param.Sp_initial,
                self.variables["Terminal voltage [V]"]: param.V_initial,
                self.variables["Shuttled Sulfur [g]"] : param.Ss_initial
            }
        )

        ######################################
        # Model events
        ######################################
        tol = 1e-4
        self.events.append(
            pybamm.Event(
                "Minimum voltage",
                V - self.param.voltage_low_cut,
                pybamm.EventType.TERMINATION,
            )
        )
        self.events.append(
            pybamm.Event(
                "Maximum voltage",
                V - self.param.voltage_high_cut,
                pybamm.EventType.TERMINATION,
            )
        )
        #self.events.append(
            #pybamm.Event(
                #"Zero theoretical capacity", cth - tol, pybamm.EventType.TERMINATION
            #)
        #)
    @property
    def default_parameter_values(self):
        # TODO: separate parameters out by component and create a parameter set
        # that can be called (see pybamm/parameters/parameter_sets.py)
        file = "models/inputs/parameters/lithium-sulfur/marinescu2018_parameters.csv"
        values_path = pybamm.get_parameters_filepath(file)
        return pybamm.ParameterValues(values=values_path)
    
    
    def print_bibtex(self):
        bibtex = '''
                @article{Marinescu2018,
author = {Marinescu, Monica and O'Neill, Laura and Zhang, Teng and Walus, Sylwia and Wilson, Timothy E. and Offer, Gregory J.},
doi = {10.1149/2.0171801jes},
issn = {0013-4651},
journal = {Journal of The Electrochemical Society},
pages = {A6107--A6118},
title = {{Irreversible vs Reversible Capacity Fade of Lithium-Sulfur Batteries during Cycling: The Effects of Precipitation and Shuttle}},
volume = {165},
year = {2018}
}
                '''
        print(bibtex)
        
    def print_citations(self):
        citation = '''Marinescu, M. et al. Irreversible vs Reversible Capacity Fade of Lithium-Sulfur Batteries during Cycling: The Effects of Precipitation and Shuttle. J. Electrochem. Soc. 165, A6107–A6118 (2018).'''
        
        print(citation)  