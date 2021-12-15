import pybamm
from .base_lithium_sulfur_model import BaseModel


class Cornish2021(BaseModel):
    """
    Zero Dimensional with the following electrochemistry
    
    S_{8}^{0} + 4e^{-} => 2 S_{4}^{2-}
    S_{4}^{2-} + 2e^{-} => 2S_{2}^{2-}
    S_{2}^{2-} + 2e^{-} => 2S^{2-}

    Parameters
    ----------
    options : dict, optional
        A dictionary of options to be passed to the model.
    name : str, optional
        The name of the model.
    """

<<<<<<< HEAD:models/Cornish2021.py
    def __init__(self, options = None, temp_control = 30, name="Cornish & Marinescu (2021) Zero Dimensional Model"):
=======
    def __init__(self, options=None, name="Cornish & Marinescu (2021) Zero Dimensional Model"):
>>>>>>> master:models/ZeroD_Cornish_type.py
        super().__init__(options, name)
        
        # citations
        #pybamm.citations.register("Cornish2021")
        
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

        #######################################
        # Model parameters
        #######################################
        from .parameters.cornish2021_parameters import Cornish2021Parameters
        self.param = Cornish2021Parameters()
        param = self.param

        # standard parameters
        R = param.R
        F = param.F
        T = param.T_ref

        # parameters
        Ms = param.Ms
        ns = param.ns
        ns2 = param.ns2
        ns4 = param.ns4
        ns8 = param.ns8
        ne = param.ne
        ih0 = param.ih0
        im0 = param.im0
        il0 = param.il0
        rho_s = param.rho_s
        EH0 = param.EH0
        EM0 = param.EM0
        EL0 = param.EL0
        v = param.v
        ar = param.ar
        k_p = param.k_p
        k_d = param.k_d
        S_star = param.S_star
        k_s_charge = param.k_s_charge
        k_s_discharge = param.k_s_discharge
        
        
        # parameters derived from other parameters
        nH = 1
        nM = 1
        nL = 1
        iH_coef = nH * F / (2 * R * T)
        iM_coef = nM * F / (2 * R * T)
        iL_coef = nL * F / (2 * R * T)
        E_H_coef = R * T / (nH * F)
        E_M_coef = R * T / (nM * F)
        E_L_coef = R * T / (nL * F)
        f_h = ((ns4 ** 2) * Ms * v / ns8)**(1/4)
        f_m = ((ns2 ** 2) * Ms * v / ns4)**(1/2)
        f_l = ((ns ** 2) * Ms * v / ns2)**(1/2)

        #######################################################
        # Non-dynamic model functions
        #######################################################

        # High plateau potenital [V] as defined by equation (2a) in [1]
        E_H = EH0 + E_H_coef * (pybamm.log(f_h) + 0.25*pybamm.log(S8) - 0.5*pybamm.log(S4) )
        
        E_M = EM0 + E_M_coef * (pybamm.log(f_m) + 0.5*pybamm.log(S4) - pybamm.log(S2) ) 

        # Low plateau potenital [V] as defined by equation (2b) in [1]
        E_L = EL0 + E_L_coef * (pybamm.log(f_l) + 0.5*pybamm.log(S2) - pybamm.log(S) ) 

        # High plateau over-potenital [V] as defined by equation (6a) in [1]
        eta_H = V - E_H
        
        eta_M = V - E_M

        # Low plateau over-potenital [V] as defined by equation (6b) in [1]
        eta_L = V - E_L
        
        # High plateau current [A] as defined by equation (5a) in [1]
        i_H = -2 * ih0 * ar * pybamm.sinh(iH_coef * eta_H)
        
        i_M = -2 * im0 * ar * pybamm.sinh(iM_coef * eta_M)

        # Low plateau current [A] as defined by equation (5b) in [1]
        i_L = -2 * il0 * ar * pybamm.sinh(iL_coef * eta_L)

        # Theoretical capacity [Ah] of the cell as defined by equation (2) in [2]
        cth = (3 * ne * F * S8 / (ns8 * Ms) + ne * F * S4 / (ns4 * Ms)) / 3600

        # Shuttle coefficient
        k_s = ( k_s_charge * (I < 0) ) + ( k_s_discharge * (I >= 0) )
        
        # Precipitation/Dissolution rate
        k_pd = ( k_d * (I < 0) ) + ( k_p * (I >= 0) )
        
        ###################################
        # Dynamic model functions
        ###################################

        # Algebraic constraint on currents as defined by equation (7) in [1]
        algebraic_condition = i_H + i_M + i_L - I
        self.algebraic.update({V: algebraic_condition})

        # Differential equation (8a) in [1]
        dS8dt = -(ns8 * Ms * i_H / (nH * F)) - k_s * S8

        # Differential equation (8b) in [1]
        dS4dt = (ns8 * Ms * i_H / (nH * F)) + k_s * S8 - (ns4 * Ms * i_M / (nM * F))

        # Differential equation (8c) in [1]
        dS2dt = ns4 * Ms * i_M / (nM * F) - (ns2 * Ms * i_L / (nM * F)) 

        # Differential equation (8d) in [1]
        dSdt = (ns2 * Ms * i_L / (nM * F)) - k_pd * Sp * (S - S_star) / (v * rho_s)

        # Differential equation (8e) in [1]
        dSpdt = k_pd * Sp * (S - S_star) / (v * rho_s)

        self.rhs.update({S8: dS8dt, S4: dS4dt, S2: dS2dt, S: dSdt, Sp: dSpdt})

        ##############################
        # Model variables
        #############################

        self.variables.update(
            {
                "Time [s]": pybamm.t * self.timescale,
                "Capacity [Ah]": pybamm.t * self.timescale * pybamm.AbsoluteValue(I) / 3600,
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
                "Middle plateau over-potential [V]": eta_M,
                "Low plateau over-potential [V]": eta_L,
                "High plateau current [A]": i_H,
                "Middle plateau current [A]": i_M,
                "Low plateau current [A]": i_L,
                "Theoretical capacity [Ah]": cth,
                "Algebraic condition": algebraic_condition,
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
                self.variables['Terminal voltage [V]']: param.V_initial
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
<<<<<<< HEAD:models/Cornish2021.py
        
        
        self.temp_control = temp_control
    @property
    def default_parameter_values(self):
        
        if self.temp_control == 20:
            file = "models/inputs/parameters/lithium-sulfur/cornish2021_parameters_20.csv"
        elif self.temp_control == 30:
            file = "models/inputs/parameters/lithium-sulfur/cornish2021_parameters_30.csv"
        elif self.temp_control == 40:
            file = "models/inputs/parameters/lithium-sulfur/cornish2021_parameters_40.csv"
        else:
            print('No defined standard parameter set for temperature control value of {} [C]'.format(temp_control))
            
        #file = "models/inputs/parameters/lithium-sulfur/cornish2021_parameters.csv"
=======
        #self.events.append(
        #    pybamm.Event(
        #        "Zero theoretical capacity", cth - tol, pybamm.EventType.TERMINATION
        #    )
        #)
        
    @property
    def default_parameter_values(self):
        # TODO: separate parameters out by component and create a parameter set
        # that can be called (see pybamm/parameters/parameter_sets.py)
        file = "models/inputs/parameters/lithium-sulfur/cornish2021_parameters.csv"
>>>>>>> master:models/ZeroD_Cornish_type.py
        values_path = pybamm.get_parameters_filepath(file)
        return pybamm.ParameterValues(values=values_path)