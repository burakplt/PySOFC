#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""utility.py: This module contains methods to calculate some thermodynamic, electrochemical and diffusion properties."""
__author__ = "Burak Polat"
__email__ = "burakplt@outlook.com"
__license__ = "MIT"

import matplotlib.pyplot as plt
import numpy as np
from math import log, pi, exp

F = 96485   # C/mol
R = 8.3144  # J/mol.K
molecular_weights = {"H2":2.016, "H2O": 18.015, "O2":31.999, "N2":28.013}     # g/mol
diffusion_volumes_list = {"H2": 6.12, "N2": 18.5, "O2": 16.3, "H2O": 13.1}    # VDI-Wärmeatlas, 2002.

class BasicThermo:
    """
    Thermodynamic calculations for SOFC system using polynomial gas phase heat capacity correlations.
    Cp(T) = A + B.T + C.T^2 + D.T^3 in J/mol.K
    H(T) = H_ref + A.T + B/2.T^2 + C/3.T^3 + D/4.T^4 in J/mol
    S(T) = S_ref + A.ln(T) + B.T + C/2.T^2 + D/3.T^3 in J/mol.K
    G(T) = H(T) - T.S(T)
    E°(T) = -dG/(2.F) where dG is the gibbs free energy of reaction
    E_Nernst = E° + RT/2F ln(p_H2.p_O2^0.5 / p_H2O)
    """

    # gas_heat_capacity_coefficients(dict) coefficients are given as "Specie" : [A, B, C, D]
    gas_heat_capacity_coefficients = {"H2":[27.124, 9.267e-3, -13.799e-6, 7.64e-9],
                                      "H2O":[32.22, 1.923e-3, 10.548e-6, -3.594e-9],
                                      "O2":[28.087, -0.004e-3, 17.447e-6, -10.644e-9],
                                      "CO":[30.848, -12.84e-3, 27.87e-6, -12.7e-9],
                                      "CO2":[19.78, 73.39e-3, -55.98e-6, 17.14e-9],
                                      "CH4":[19.238, 52.09e-3, 11.966e-6, -11.309e-9],
                                      "N2":[31.128, -13.556e-3, 26.777e-6, -11.673e-9]}

    # Correction factors to match standard state enthalpy and entropy for consistency
    gas_phase_correction_enthalpy = {"H2":8392, "H2O": 9778+241824, "O2":8507, "CO":8847.7+110530, "CO2":8693.08+393510,
                                     "CH4":8129.08+74870, "N2":8887}
    gas_phase_correction_entropy = {"H2": -26.05, "H2O": 4.25, "O2": 44.4, "CO":24.61, "CO2":81.57, "CH4":60.64,
                                    "N2":17.22}

    def __init__(self):
        pass

    @staticmethod
    def enthalpy_gas(specie: str, T: float) -> float:
        """
        Gas state enthalpy of gas at given temperature.
        :param specie: Options are "H2", "H2O", "O2", "N2", "CO", "CO2", "CH4"
        :param T: Temperature [K]
        :return: Enthalpy [J/mol]
        """
        k = BasicThermo.gas_heat_capacity_coefficients[specie]
        A = k[0]
        B = k[1]
        C = k[2]
        D = k[3]
        h_ref = BasicThermo.gas_phase_correction_enthalpy[specie]   # J/mol at T_ref = 298.15 K

        h = T*A + T**2 *B/2 + T**3 *C/3 + T**4 *D/4    # J/mol
        return h - h_ref

    @staticmethod
    def entropy_gas(specie: str, T: float) -> float:
        """
        Gas state entropy of gas at given temperature.
        :param specie: Options are "H2", "H2O", "O2", "N2", "CO", "CO2", "CH4"
        :param T: Temperature [K]
        :return: Entropy [J/mol.K]
        """
        k = BasicThermo.gas_heat_capacity_coefficients[specie]
        A = k[0]
        B = k[1]
        C = k[2]
        D = k[3]
        s_correction = BasicThermo.gas_phase_correction_entropy[specie]   # J/mol.K at T_ref = 298.15 K

        s = A*log(T) + T*B + T**2 *C/2 + T**3 *D/3    # J/mol.K
        return s + s_correction

    @staticmethod
    def gibbs_gas(specie: str, T: float) -> float:
        """
        Gas state gibbs free energy of gas at given temperature.
        :param specie: Options are "H2", "H2O", "O2", "N2", "CO", "CO2", "CH4"
        :param T: Temperature [K]
        :return: Gibbs free energy [J/mol]
        """
        g = BasicThermo.enthalpy_gas(specie, T) - T * BasicThermo.entropy_gas(specie, T)
        return g

    @staticmethod
    def equilibrium_constant_wgs(T: float) -> float:
        """
        Calculates equilibrium constant of water-gas shift reaction.
        :param T: Temperature [K]
        :return: Equilibrium constant [-]
        """
        G = BasicThermo.gibbs_gas
        dG = G("CO2", T) + G("H2", T) - G("CO", T) - G("H2O", T)
        k = exp(-dG/R/T)
        return k

    @staticmethod
    def std_potential_hydrogen(T: float) -> float:
        """
        Calculates standard potential for reaction H2 + 1/2 O2 --> H2O at P=1 bar and pure H2, O2 feeds.
        Generated via: BasicThermo.generate_std_potential_curve_hydrogen(380, 1300, 200, 2, show_plot=True)
        :param T:Temperature [K]
        :return: Standard potential E°(T) [Volt]
        """
        E = -2.52296740e-8* T**2 - 2.35159319e-4*T + 1.25893994
        return E

    @staticmethod
    def std_potential_CO(T: float) -> float:
        """
        Calculates standard potential for reaction CO + 1/2 O2 --> CO2 at P=1 bar and pure CO, O2 feeds.
        Generated via: BasicThermo.generate_std_potential_curve_CO(380, 1300, 200, 2, show_plot=True)
        :param T:Temperature [K]
        :return: Standard potential E°(T) [Volt]
        """
        E = 7.80893469e-9 * T **2 - 4.66575523e-4 * T + 1.47245339
        return E

    @staticmethod
    def nernst_potential_hydrogen(T: float, P: float, x_H2: float, x_O2: float, x_H2O: float) -> float:
        """
        Nernst potential of Hydrogen/Water fueled SOFC. E_Nernst = E° + RT/2F ln(p_H2.p_O2^0.5 / p_H2O)

        :param T: Temperature [K]
        :param x_H2: Mole fraction of Hydrogen
        :param x_O2: Mole fraction of Oxygen
        :param x_H2O: Mole fraction of Water
        :return: Nernst potential E_Nernst(T) [Volt]
        """
        if x_H2O == 0:
            E_nernst = BasicThermo.std_potential_hydrogen(T) + R * T / 4 / F * log(P/1)
        else:
            E_nernst = BasicThermo.std_potential_hydrogen(T) + R * T / 2 / F * \
                       log(x_H2 * x_O2 ** 0.5 / x_H2O * P**0.5) + R * T / 4 / F * log(P/1)
        return E_nernst

    @staticmethod
    def nernst_potential_CO(T: float, P: float, x_CO: float, x_O2: float, x_CO2: float) -> float:
        """
        Nernst potential of Carbon monoxide fueled SOFC. E_Nernst = E° + RT/2F ln(p_CO.p_O2^0.5 / p_CO2)

        :param T: Temperature [K]
        :param x_CO: Mole fraction of Carbon monoxide
        :param x_O2: Mole fraction of Oxygen
        :param x_CO2: Mole fraction of Carbon dioxide
        :return: Nernst potential E_Nernst(T) [Volt]
        """
        if x_CO2 == 0:
            E_nernst = BasicThermo.std_potential_CO(T) + R * T / 4 / F * log(P / 1)
        else:
            E_nernst = BasicThermo.std_potential_CO(T) + R * T / 2 / F * \
                       log(x_CO * x_O2 ** 0.5 / x_CO2 * P ** 0.5) + R * T / 4 / F * log(P / 1)
        return E_nernst

    @staticmethod
    def generate_std_potential_curve_hydrogen(T_min: float, T_max: float, sample_size: int,
                                              polynom_order: int = 2, show_plot: bool = False):
        """
        Generates a temperature dependent standard potential E°(T) for a given temperature range, sample size and
        the order of polynomial to be fitted. An already fitted equation using this function is employed in this class.
        It is valid only for the reaction H2 + 1/2 O2 --> H2O at P=1 bar and produced water is in gas phase.

        :param T_min: Temperature lower limit in K
        :param T_max: Temperature upper limit in K
        :param sample_size: Number of points to divide temperature range
        :param polynom_order: Order of polynomial to fit Temperature-OCV data
        :param show_plot: If True, a plot with calculated data points and fitted curve is shown.
        :return: Coefficients of the fitted polynomial as numpy array, the highest power is first element of array.
        """

        T_range = np.linspace(T_min, T_max, sample_size)    # Temperature interval
        potential = []
        for T in T_range:
            # Calculate the change in gibbs free energy *dG* of reaction H2 + 1/2 O2 --> H2O
            dG = BasicThermo.gibbs_gas("H2O", T) - BasicThermo.gibbs_gas("H2", T) - BasicThermo.gibbs_gas("O2", T) / 2
            potential.append(-dG / 2 / F)     # Potential in Volt

        coefficients = np.polyfit(T_range, potential, polynom_order)
        if show_plot:
            f = np.poly1d(coefficients)
            plt.plot(T_range, potential,"*", label="E°(T)")
            plt.plot(T_range, f(T_range), label="Fitted curve")
            plt.legend()
            plt.xlabel("Temperature [K]")
            plt.ylabel("Standard potential E° [Volt]")
            plt.show()

        return coefficients

    @staticmethod
    def generate_std_potential_curve_CO(T_min: float, T_max: float, sample_size: int,
                                              polynom_order: int = 2, show_plot: bool = False):
        """
        Generates a temperature dependent standard potential E°(T) for a given temperature range, sample size and
        the order of polynomial to be fitted. An already fitted equation using this function is employed in this class.
        It is valid only for the reaction CO + 1/2 O2 --> CO2 at P=1 bar.

        :param T_min: Temperature lower limit in K
        :param T_max: Temperature upper limit in K
        :param sample_size: Number of points to divide temperature range
        :param polynom_order: Order of polynomial to fit Temperature-OCV data
        :param show_plot: If True, a plot with calculated data points and fitted curve is shown.
        :return: Coefficients of the fitted polynomial as numpy array, the highest power is first element of array.
        """

        T_range = np.linspace(T_min, T_max, sample_size)    # Temperature interval
        potential = []
        for T in T_range:
            # Calculate the change in gibbs free energy *dG* of reaction H2 + 1/2 O2 --> H2O
            dG = BasicThermo.gibbs_gas("CO2", T) - BasicThermo.gibbs_gas("CO", T) - BasicThermo.gibbs_gas("O2", T) / 2
            potential.append(-dG / 2 / F)     # Potential in Volt

        coefficients = np.polyfit(T_range, potential, polynom_order)
        if show_plot:
            f = np.poly1d(coefficients)
            plt.plot(T_range, potential,"*", label="E°(T)")
            plt.plot(T_range, f(T_range), label="Fitted curve")
            plt.legend()
            plt.xlabel("Temperature [K]")
            plt.ylabel("Standard potential E° [Volt]")
            plt.show()

        return coefficients


class PolarizationFunctions:
    """
    This class contains some static methods for calculation of polarization functions (voltage loss).
    """

    def __init__(self):
        pass

    @staticmethod
    def knudsen_diffusivity_explicit(specie: str, pore_diameter: float, temperature: float) -> float:
        """
        Calculates Knudsen diffusivity coefficient D_Kn explicitly for given specie and temperature.
        :param specie: Options are "H2", "O2", "H2O", "N2"
        :param pore_diameter: Pore diameter of the electrode [m]
        :param temperature: Temperature [K]
        :return: Knudsen diffusivity D_Kn [m^2/s]
        """
        MW = molecular_weights[specie]   # g/mol
        D_Kn = pore_diameter/3 * ( 8*R*temperature/pi/MW )**0.5   # D_Kn in [m^2/s]. ? 1000 is mass conversion factor
        return D_Kn

    @staticmethod
    def calculate_knudsen_factors(pore_diameter: float):
        """
        Calculates coefficients *C* to express Knudsen diffusivity in form D_Kn(T) = C. T^0.5
        It reduces computational repetition, and is called whenever a sofc.SolidOxidFuelCell object is instantiated.
        :param pore_diameter: Electrode pore diameter [m]
        :return: dict object containing coefficients {"H2":C1, "O2":C2, "N2":C3, "H2O":C4}
        """
        f = lambda specie: pore_diameter / 3 * (8 * 1000 * R / pi / molecular_weights[specie]) ** 0.5
        constants = {"H2": f("H2"), "O2": f("O2"), "N2": f("N2"), "H2O": f("H2O")}  # m^2/s K^0.5
        return constants

    @staticmethod
    def binary_diffusion_explicit(specie_1: str, specie_2:str, temperature: float, pressure: float,
                                  diffusion_volume=None) -> float:
        """
        Binary diffusion coefficient for given species at given temperature and pressure.
        :param specie_1: Options are "H2", "O2", "N2", "H2O"
        :param specie_2: Options are "H2", "O2", "N2", "H2O"
        :param temperature: Temperature [K]
        :param pressure: Total pressure [bar]
        :param diffusion_volume: Diffusion volumes, by default: {"H2": 6.12, "N2": 18.5, "O2": 16.3, "H2O": 13.1}
        :return: Binary diffusion coefficient D_ij in [m^2/s]
        """
        if not diffusion_volume:
            diffusion_volume = diffusion_volumes_list

        D = 1.43e-7 * temperature ** (1.75) / (
                    2 / (1 / molecular_weights[specie_1] + 1 / molecular_weights[specie_2])) ** 0.5 / \
            (diffusion_volume[specie_1] ** (1 / 3) + diffusion_volume[specie_2] ** (1 / 3)) ** 2 / pressure
        return D

    @staticmethod
    def calculate_binary_diffusion_factors(diffusion_volumes=None):
        """
        Calculates coefficients *C* to express binary Diffusion coefficient in form D_ij(T) = C. T^1.75 / P
        :param diffusion_volumes: Diffusion volumes, by default: {"H2": 6.12, "N2": 18.5, "O2": 16.3, "H2O": 13.1}
        :return: dict object containing coefficients {"H2-H2O":C1, "O2-N2":C2}
        """
        if not diffusion_volumes:
            diffusion_volume = diffusion_volumes_list

        D_H2_H2O = 1.43e-7 / ( 2 / (1/molecular_weights["H2"] + 1/molecular_weights["H2O"]) )** 0.5 / \
            (diffusion_volume["H2"] ** (1 / 3) + diffusion_volume["H2O"] ** (1 / 3)) ** 2

        D_O2_N2 = 1.43e-7 / ( 2 / (1/molecular_weights["O2"] + 1/molecular_weights["N2"]) )** 0.5 / \
                   ( diffusion_volume["O2"]** (1 / 3) + diffusion_volume["N2"]** (1 / 3) ) ** 2

        return {"H2-H2O": D_H2_H2O, "O2-N2": D_O2_N2}
