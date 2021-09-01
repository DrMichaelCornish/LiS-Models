import pybamm
from .base_lithium_sulfur_model import BaseModel



class ZeroD_Marinescu_type(BaseModel):
    """
    Zero Dimensional model with Chemistry 1
    
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
    .. [1]  Marinescu, M., Zhang, T. & Offer, G. J. (2016).
            A zero dimensional model of lithium-sulfur batteries during charge
            and discharge. Physical Chemistry Chemical Physics, 18, 584-593.

    .. [2]  Marinescu, M., O’Neill, L., Zhang, T., Walus, S., Wilson, T. E., &
            Offer, G. J. (2018).
            Irreversible vs reversible capacity fade of lithium-sulfur batteries
            during cycling: the effects of precipitation and shuttle. Journal of
            The Electrochemical Society, 165(1), A6107-A6118.
            
    .. [3]  Hua, X., Zhang, T., Offer, G., & Marinescu, M. (2019).
            Towards online tracking of the shuttle effect in lithium sulfur 
            batteries using differential thermal voltemetry. Journal of Energy 
            Storage, 21 (2019), 765-772.
    """

    def __init__(self, options=None, name="Zero Dimensional Marinescu Type"):
        super().__init__(options, name)
        
        # citations
        #pybamm.citations.register("Marinescu2016")
        #if 'degradation' in options:
        #    pybamm.citations.register("Marinescu2018")
        #if 'thermal' in options:
        #    pybamm.citations.register("Hua2019")
        
        # set base options
        if 'degradation' not in self.options.keys():
            self.options['degradation'] = 'none'
            
        if 'thermal' not in self.options.keys():
            self.options['thermal'] = 'none'
            
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
        if self.options['degradation'] == 'Marinescu2018':
            Ss = pybamm.Variable("Shuttled Sulfur [g]")
        if self.options['thermal'] == 'Hua2019':
            Tc = pybamm.Variable("Cell Temperature [K]")

        #######################################
        # Model parameters as defined in table (1) in [1]. Parameters with 'H' or
        # 'L' in the name represent the high and low plateau parameter, respectively.
        #######################################
        param = self.param

        # standard parameters
        R = param.R
        F = param.F
        T = param.T_ref
        Ta = param.Ta
        if self.options['thermal'] == 'none':
            Tc = T
        N = param.N

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

        # High plateau potenital [V] as defined by equation (2a) in [1]
        E_H = EH0 + E_H_coef * ( pybamm.log(f_h) + pybamm.log(S8) - 2*pybamm.log(S4) )

        # Low plateau potenital [V] as defined by equation (2b) in [1]
        E_L = EL0 + E_H_coef * ( pybamm.log(f_l) + pybamm.log(S4) - pybamm.log(S2) - 2*pybamm.log(S) ) 

        # High plateau over-potenital [V] as defined by equation (6a) in [1]
        eta_H = V - E_H

        # Low plateau over-potenital [V] as defined by equation (6b) in [1]
        eta_L = V - E_L

        # High plateau current [A] as defined by equation (5a) in [1]
        i_H = -2 * ih0 * ar * pybamm.sinh(i_coef * eta_H)

        # Low plateau current [A] as defined by equation (5b) in [1]
        i_L = -2 * il0 * ar * pybamm.sinh(i_coef * eta_L)

        # Theoretical capacity [Ah] of the cell as defined by equation (2) in [2]
        cth = (3 * ne * F * S8 / (ns8 * Ms) + ne * F * S4 / (ns4 * Ms)) / 3600

        # Shuttle coefficient (set Tc = T and H = 0 to retrieve Marinescu et al. (2016,2018)
        k_s =  k_s_discharge * (I >= 0) + k_s_charge *  (I < 0) * pybamm.exp(-A*N*((1/Tc)-(1/T))/R)

        ###################################
        # Dynamic model functions
        ###################################

        # Algebraic constraint on currents as defined by equation (7) in [1]
        algebraic_condition = i_H + i_L - I
        self.algebraic.update({V: algebraic_condition})

        # Differential equation (8a) in [2]
        dS8dt = -(ns8 * Ms * i_H / (ne * F)) - k_s * S8

        if self.options["degradation"] == 'Marinescu2018': 
            # Differential equation (8b) in [2] 
            dS4dt = (ns8 * Ms * i_H / (ne * F)) + (1- (f_s/Ms)*Ss)* k_s * S8 - (ns4 * Ms * i_L / (ne * F))
            
            # Differential equation (8f) in [2]
            dSsdt = k_s * S8
        else:
            # Differential equation (8b) in [1] 
            dS4dt = (ns8 * Ms * i_H / (ne * F)) + k_s * S8 - (ns4 * Ms * i_L / (ne * F))

        # Differential equation (8c) in [2]
        dS2dt = ns2 * Ms * i_L / (ne * F)

        # Differential equation (8d) in [2]
        dSdt = (2 * ns * Ms * i_L / (ne * F)) - k_p * Sp * (S - S_star) / (v * rho_s)
    
        # Differential equation (8e) in [2]
        dSpdt = k_p * Sp * (S - S_star) / (v * rho_s)
        
        self.rhs.update({S8: dS8dt, S4: dS4dt, S2: dS2dt, S: dSdt, Sp: dSpdt})#
        
        if self.options["degradation"] == 'Marinescu2018':
            # Differential equation (8f) in [2]
            dSsdt = k_s * S8
            
            self.rhs.update({Ss : dSsdt})
            
        if self.options["thermal"] == 'Hua2019':
            # Differential equation (8g) in Supplementary material A of [3]
            dTcdt = ((k_s*ne*F*S8*V/(ns8*Ms)) - h*(Tc-Ta))/(m_c*c_h)
            
            self.rhs.update({Tc : dTcdt})
        
        ##############################
        # Model variables
        #############################
        
        self.variables.update(
            {
                "Time [s]": pybamm.t * self.timescale,
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
                "Theoretical capacity [Ah]": cth,
                "Algebraic condition": algebraic_condition,
            }
        )
        
        if self.options['degradation'] == 'Marinescu2018':
            self.variables.update({"Shuttled Sulfur [g]" : Ss})
        if self.options['thermal'] == 'Hua2019':
            self.variables.update({"Cell Temperature [K]" : Tc})

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
                self.variables["Terminal voltage [V]"]: param.V_initial
            }
        )
        
        if self.options['degradation'] == 'Marinescu2018':
            self.initial_conditions.update({self.variables["Shuttled Sulfur [g]"] : param.Ss_initial})
        if self.options['thermal'] == 'Hua2019':
            self.initial_conditions.update({self.variables["Cell Temperature [K]"] : param.Tc_initial})

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
        self.events.append(
            pybamm.Event(
                "Zero theoretical capacity", cth - tol, pybamm.EventType.TERMINATION
            )
        )
