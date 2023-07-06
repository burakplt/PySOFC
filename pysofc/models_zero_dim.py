#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""models_zero_dim.py: This module contains zero-dimensional Solid Oxide Fuel Cell simulation."""
__author__ = "Burak Polat"
__email__ = "burakplt@outlook.com"
__license__ = "MIT"

from .utility import *
import numpy as np
import matplotlib.pyplot as plt


class LumpedModelHydrogenSOFC:

    def __init__(self, fuel_cell):
        """
        Modeling Solid Oxide Fuel and Electrolysis Cell by assuming lumped parameter (0-Dimensional) model.
        :param fuel_cell: Related SOC object. Required for cell parameters and properties.
        """
        self.FC = fuel_cell
        self.input = {"T_fuel_in": None, "T_oxygen_in": None, "pressure": None, "T_operation": None,
                     "T_fuel_out": None, "T_oxygen_out": None, "pressure_drop": None,
                     "X_fuel_in": None, "X_oxygen_in": None, "X_fuel_out": None, "X_oxygen_out": None,
                     "fuel_utilization": None, "oxygen_utilization": None, "fuel_to_air_ratio": None,
                     "molar_flowrate_fuel": None, "molar_flowrate_oxygen": None,
                     "operation_voltage": None, "mean_current_density": None}

        self.input_units = {"T_fuel_in": "K", "T_oxygen_in": "K", "pressure": "bar", "T_operation": "[K]",
                     "T_fuel_out": "K", "T_oxygen_out": "K", "pressure_drop": "bar",
                     "X_fuel_in": "dict", "X_oxygen_in": "dict", "X_fuel_out": "dict", "X_oxygen_out": "dict",
                     "fuel_utilization": "-", "oxygen_utilization": "-", "fuel_to_air_ratio": "-",
                     "molar_flowrate_fuel": "mol/s", "molar_flowrate_oxygen": "mol/s",
                     "operation_voltage": "Volt", "mean_current_density": "A/m^2"}

        self.temperature_rise_fuel = 50     # [K] is used to set operation temperature if T_op not specified
        self.T_op = None
        self.N_oxygen_in = None
        self.max_current = None

    def print_parameter_list(self) -> None:
        """
        Prints current operation parameters (self.input) alongside with their units.
        :return: None
        """
        for param in list(self.input.keys()):
            print(f"{param}: in [{self.input_units[param]}]")

    def reset_parameter_list(self) -> None:
        """
        Sets all entry in input parameters (self.input) as None
        :return: None
        """
        for param in list(self.input.keys()):
            self.input[param] = None

    def evaluate_input_data(self) -> None:
        """
        Pre-evaluation of input data. Input list is extensive, means all operation conditions are not required for a
        particular case. For instance operation temperature can be given directly as "T_operation", or can be calculated
        using "T_fuel_in" and "T_fuel_out".
        :return: None
        """

        # Operation temperature T_op is given and constant throughout SOFC operation
        if self.input["T_operation"]:
            self.T_op = self.input["T_operation"]

        # If inlet and outlet temperatures of fuel side are given
        elif self.input["T_fuel_in"] and self.input["T_fuel_out"]:
            self.T_op = (self.input["T_fuel_in"] + self.input["T_fuel_out"]) / 2

        # If T_op and T_fuel_out not specified: T_op = T_fuel_in + temperature_rise_fuel
        else:
            self.T_op = self.input["T_fuel_in"] + self.temperature_rise_fuel

        if self.input["molar_flowrate_oxygen"]:
            self.N_oxygen_in = self.input["molar_flowrate_oxygen"]  # Inlet molar flow rate of oxidant
        elif self.input["fuel_to_air_ratio"]:
            self.N_oxygen_in = self.input["molar_flowrate_fuel"] / self.input["fuel_to_air_ratio"]
        else:
            x_H2_in = self.input["X_fuel_in"]["H2"]         # Inlet mole fraction of Hydrogen
            x_O2_in = self.input["X_oxygen_in"]["O2"]       # Inlet mole fraction of Oxygen
            self.N_oxygen_in = self.input["molar_flowrate_fuel"] * x_H2_in / 2 / x_O2_in
            print(f"Inlet flow rate of oxidant was not given. It is set stoichiometric flow rate: {self.N_oxygen_in} mol/s")

        self.max_current = self.input["molar_flowrate_fuel"] * self.input["X_fuel_in"]["H2"] * 2 * F

    def calculate_operation_voltage(self, current_density=None, temperature=None, fuel_utilisation=None) -> float:
        """
        Calculates operational voltage of SOC at an operation point.
        :param current_density: Optional. If not given either fuel utilisation or mean_current_density required
        :param temperature: Optional. If not given, self.T_op is used
        :param fuel_utilisation: Optional. Fuel utilisation is not used if current density is known.
        :return: Operating Voltage [V]
        """
        self.evaluate_input_data()
        x_H2_in = self.input["X_fuel_in"]["H2"]                     # Inlet mole fraction of Hydrogen
        x_O2_in = self.input["X_oxygen_in"]["O2"]                   # Inlet mole fraction of Oxygen
        x_H2O_in = self.input["X_fuel_in"]["H2O"]                   # Inlet mole fraction of Water
        N_fuel_in = self.input["molar_flowrate_fuel"]               # Inlet molar flow rate of fuel
        N_oxygen_in = self.N_oxygen_in                              # Inlet molar flow rate of oxidant
        A_cell = self.FC.dimensions["total_electroactive_area"]     # Total active cell area
        T_op = temperature if temperature else self.T_op            # Operation temperature [K]
        P = self.input["pressure"]                                  # Pressure [bar]
        is_pure_oxygen = True if x_O2_in == 1 else False

        if current_density is not None:
            J = current_density * 1e4
            N_reacted = J * A_cell / (2 * F)
        elif fuel_utilisation is not None:
            N_reacted =  N_fuel_in * x_H2_in * fuel_utilisation
            J = N_reacted * 2 * F / A_cell
        else:
            J = self.input["mean_current_density"]
            N_reacted = J * 1e4 * A_cell / (2 * F)

        if J == 0:
            OCV = BasicThermo.nernst_potential_hydrogen(T_op, P, x_H2_in, x_O2_in, x_H2O_in)  # OCV at j = 0
            return OCV

        # Outlet molar fractions
        x_H2_out = (x_H2_in * N_fuel_in - N_reacted) / N_fuel_in
        x_O2_out = (x_O2_in * N_oxygen_in - N_reacted / 2) / (N_oxygen_in - N_reacted / 2)
        x_H2O_out = (x_H2O_in * N_fuel_in + N_reacted) / N_fuel_in

        # Average molar fractions
        x_H2 = (x_H2_in + x_H2_out) / 2
        x_H2O = (x_H2O_in + x_H2O_out) / 2
        x_O2 = (x_O2_in + x_O2_out) / 2
        V_Nernst = BasicThermo.nernst_potential_hydrogen(T_op, P, x_H2, x_O2, x_H2O)    # Nernst voltage at given J

        activation_polarization = self.FC.activation_polarization(J, T_op)              # Activation polarization
        ohmic_polarization = self.FC.ohmic_polarization(J, T_op)                        # Ohmic polarization
        dif_f = self.FC.diffusion_polarization_anode(J, T_op, P, x_H2, x_H2O)
        dif_a = self.FC.diffusion_polarization_cathode(J, T_op, P, x_O2, pure_oxygen=is_pure_oxygen)
        diffusion_polarization = dif_a + dif_f  # Diffusion polarization

        V_op = V_Nernst - activation_polarization - ohmic_polarization - diffusion_polarization  # Cell voltage
        return V_op

    def calculate_current_density(self, operation_voltage=None):
        """
        Calculates current density at a given operating voltage.
        :param operation_voltage: Optional. If not given, operation_voltage input is used.
        :return: Current density [A/cm^2] or None if numeric method not converged
        """
        self.evaluate_input_data()
        V_op = operation_voltage if operation_voltage else self.input["operation_voltage"]

        f = lambda J : V_op - self.calculate_operation_voltage(current_density = J/1e4)
        j0 = 0
        j1 = self.max_current / self.FC.dimensions["total_electroactive_area"]  # Max current density

        for i in range(20):
            j = j0 - (j1-j0) * f(j0) / ( f(j1) - f(j0) )
            j0 = j1
            j1 = j
            if abs( f(j) ) < 1e-5:
                return j/1e4

        print("Secant method couldn't converge in 20 steps. Return None!")
        return None

    def isothermal_voltage_current_curve(self, J_list, experiment_J=None, experiment_V=None, experiment_name=None, show_plot=True):
        """
        Generates Current density-Voltage curve at given operational conditions. Current density is independent variable.
        :param J_list: (list) Current density list [A/m^2]
        :param experiment_J: (list) Optional. Current density list of Experimental V-I data [A/m^2]
        :param experiment_V: (list) Optional. Operating voltage list of Experimental V-I data [Volt]
        :param experiment_name: (str) Name or identifier of experimental data
        :param show_plot: Result plots is showed immediately if True
        :return: {"J":current_density, "V":voltage, "ASR":asr, "FU":fu, "POWER":power, "ERROR":relative_error}
        """
        self.evaluate_input_data()
        x_H2_in = self.input["X_fuel_in"]["H2"]                     # Inlet mole fraction of Hydrogen
        x_O2_in = self.input["X_oxygen_in"]["O2"]                   # Inlet mole fraction of Oxygen
        x_H2O_in = self.input["X_fuel_in"]["H2O"]                   # Inlet mole fraction of Water
        N_fuel_in = self.input["molar_flowrate_fuel"]               # Inlet molar flow rate of fuel
        N_oxygen_in = self.N_oxygen_in                              # Inlet molar flow rate of oxidant
        A_cell = self.FC.dimensions["total_electroactive_area"]     # Total active cell area
        T_op = self.T_op                                            # Operation temperature [K]
        P = self.input["pressure"]                                  # Pressure [bar]

        is_pure_oxygen = True if x_O2_in == 1 else False
        OCV = BasicThermo.nernst_potential_hydrogen(T_op, P, x_H2_in, x_O2_in, x_H2O_in)     # OCV at j = 0

        # Define parameter arrays
        voltage, current_density = [], []
        act, dif_fuel, dif_oxy, ohm, power, fu, asr = [], [], [], [], [], [], []
        
        for J in J_list:
            J *= 1e+4
            I = J * A_cell                          # Current [A]
            N_reacted = I / (2*F)                   # Reaction extend

            # Outlet molar fractions
            x_H2_out = (x_H2_in*N_fuel_in - N_reacted) / N_fuel_in
            x_O2_out = (x_O2_in*N_oxygen_in - N_reacted/2) / (N_oxygen_in - N_reacted/2)
            x_H2O_out = (x_H2O_in*N_fuel_in + N_reacted) / N_fuel_in

            # Average molar fractions
            x_H2 = (x_H2_in + x_H2_out) / 2
            x_H2O = (x_H2O_in + x_H2O_out) / 2
            x_O2 = (x_O2_in + x_O2_out) / 2

            V_Nernst = BasicThermo.nernst_potential_hydrogen(T_op, P, x_H2, x_O2, x_H2O) # Nernst voltage at given J

            activation_polarization = self.FC.activation_polarization(J, T_op)          # Activation polarization
            ohmic_polarization = self.FC.ohmic_polarization(J, T_op)                    # Ohmic polarization
            dif_f = self.FC.diffusion_polarization_anode(J, T_op, P, x_H2, x_H2O)
            dif_a = self.FC.diffusion_polarization_cathode(J, T_op, P, x_O2, pure_oxygen=is_pure_oxygen)
            diffusion_polarization = dif_a + dif_f                                     # Diffusion polarization

            V_op = V_Nernst - activation_polarization - ohmic_polarization - diffusion_polarization     # Cell voltage
            voltage.append(V_op)
            current_density.append(J/1e4)                   # Current density in A/cm2
            act.append(activation_polarization)
            ohm.append(ohmic_polarization)
            dif_fuel.append(dif_f)
            dif_oxy.append(dif_a)
            _power = V_op * J/1e4
            power.append(_power)
            _asr = (OCV-V_op) / J * 1e4
            asr.append(_asr)                              # ASR
            _fu = N_reacted / N_fuel_in / x_H2_in
            fu.append(_fu)                                  # Fuel utilization

        j_sim, j_exp = np.array(current_density), np.array(experiment_J)
        relative_error = None
        if j_sim.size == j_exp.size and np.abs(j_sim - j_exp).sum() < 1 :
            relative_error = (np.array(voltage) - np.array(experiment_V)) / np.array(experiment_V)
            relative_error = np.abs(relative_error) * 100

        fig, ax = plt.subplots(2, 2, sharex=True)
        ax[0,0].plot(current_density, voltage, "*--", label="Lumped Model")
        if experiment_V and experiment_J:
            ax[0,0].plot(experiment_J, experiment_V, "+--", label=f"Experiment:{experiment_name}")
        ax[0, 0].set_xlabel(r"Current Density [A/$cm^2$]")
        ax[0, 0].set_ylabel("Cell Potential [V]")

        ax[0,1].plot(current_density, act, "*-", label=r"$\eta_{act}$")
        ax[0,1].plot(current_density, ohm, "+-", label=r"$\eta_{ohm}$")
        ax[0,1].plot(current_density, dif_fuel, "o-", label=r"$\eta_{dif,fuel}$")
        ax[0,1].plot(current_density, dif_oxy, "x-", label=r"$\eta_{dif,oxy}$")
        ax[0,1].set_xlabel("Current Density [A/$cm^2$]")
        ax[0,1].set_ylabel("Polarization (Voltage drop) [V]")

        ax[1,0].plot(current_density, power, "o-", label="Power Output")
        ax[1,0].set_xlabel("Current Density [A/$cm^2$]")
        ax[1,0].set_ylabel("Power density [W/$cm^2$]")

        ax[1, 1].plot(current_density, asr, "o-", label="ASR")
        ax[1, 1].set_xlabel("Current Density [A/$cm^2$]")
        ax[1, 1].set_ylabel(r"Area Specific Resistance [$\Omega.cm^2$]")

        ax[0, 0].legend()
        ax[0, 1].legend()
        ax[1, 0].legend()
        ax[1, 1].legend()
        fig.suptitle(f"Parameters {self.FC.cell_id}, {T_op} K, {P} bar \n" + \
                     f"Inlet: $H_2$: %{x_H2_in*100}, $H_2O$: %{x_H2O_in*100}, $O_2$: %{x_O2_in*100}, "+ r"$N_{fuel,in}$: "+ f"{N_fuel_in*60:.5f} mol/min" )
        if show_plot:
            plt.show()

        return {"J":current_density, "V":voltage, "ASR":asr, "FU":fu, "POWER":power, "ERROR":relative_error}
