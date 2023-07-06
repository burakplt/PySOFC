#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""solide_oxide_cell.py: This module contains HydrogenSOFC class that represents a physical model of a Solid Oxide Fuel Cell"""
__author__ = "Burak Polat"
__email__ = "burakplt@outlook.com"
__license__ = "MIT"

from .utility import *
import json
from math import exp, asinh


class HydrogenSOFC:

    def __init__(self, cell_id: str, load_from_json=None):
        """
        Physical model of a Solid Oxide Cell.
        :param cell_id: Name or ID of a HydrogenSOFC instance
        :param load_from_json: None or path to cell parameter .json file
        """
        self.cell_id = cell_id
        self.binary_diffusion_factors = None  # Factors to express binary Diffusion as D_ij(T) = C. T^1.75 / P
        self.knudsen_factors = None           # Pore diameter-specific Knudsen factors

        self.dimensions = {
            "channel_height": None,
            "channel_width": None,
            "channel_length": None,
            "total_electroactive_area": None,
            "thickness_cathode_electrode": None,
            "thickness_anode_electrode": None,
            "thickness_electrolyte": None,
        }

        self.activation_parameters = {
            "activation_energy_fuel": None,
            "activation_energy_oxygen": None,
            "pre_exp_factor_fuel": None,
            "pre_exp_factor_oxygen": None
        }

        self.diffusion_parameters = {
            "porosity": None,
            "tortuosity": None,
            "pore_diameter": None,
            "diffusion_volumes": {"H2": 6.12, "N2": 18.5, "O2": 16.3, "H2O": 13.1}  # VDI-Wärmeatlas, 2002.
        }

        self.ohmic_parameters = {
            "activation_energy_electrolyte": None,
            "pre_exp_factor_electrolyte": None,
            "other_ohmic_resistance": None
        }

        self.parameter_units = {
            "channel_height": "m",
            "channel_width": "m",
            "channel_length": "m",
            "total_electroactive_area": "m^2",
            "thickness_cathode_electrode": "m",
            "thickness_anode_electrode": "m",
            "thickness_electrolyte": "m",
            "activation_energy_fuel": "J/mol",
            "activation_energy_oxygen": "J/mol",
            "pre_exp_factor_fuel": "A/m^2",
            "pre_exp_factor_oxygen": "A/m^2",
            "porosity": "dimensionless",
            "tortuosity": "dimensionless",
            "pore_diameter": "m",
            "diffusion_volumes": "dimensionless",
            "activation_energy_electrolyte": "J/mol",
            "pre_exp_factor_electrolyte": "1/ohm.m",
            "other_ohmic_resistance": "ohm.m^2"
        }

        # Load parameters from json file
        if load_from_json:
            self.import_parameters_json(load_from_json)
            self.pre_calculations()

    def pre_calculations(self):
        # Calculate factors to express binary Diffusion coefficient in form D_ij(T) = C. T^1.75 / P
        self.binary_diffusion_factors = PolarizationFunctions.calculate_binary_diffusion_factors()

        # Calculate pore diameter-specific Knudsen factors
        self.knudsen_factors = PolarizationFunctions.calculate_knudsen_factors(self.diffusion_parameters["pore_diameter"])

    def import_parameters_json(self, file_path: str) -> None:
        """
        Reads parameters from .json file. Structure must match the output of self.export_parameters_json function.
        :param file_path: Full path of cell parameter .json file
        :return: None
        """
        try:
            with open(file_path, 'r') as f:
                data = json.load(f)
                self.cell_id = data["cell_id"]
                self.dimensions = data["dimensions"]
                self.activation_parameters = data["activation_parameters"]
                self.diffusion_parameters = data["diffusion_parameters"]
                self.ohmic_parameters = data["ohmic_parameters"]
        except:
            raise Exception("Cell parameters file structure is not valid!")

    def export_parameters_json(self, name:str) -> None:
        """
        Writes the current cell parameters to .json file
        :param name: File name without any extensions
        :return: None
        """
        parameters = {"cell_id":name, "dimensions": self.dimensions, "activation_parameters": self.activation_parameters,
                      "diffusion_parameters": self.diffusion_parameters, "ohmic_parameters": self.ohmic_parameters}

        with open(f"Parameter Files\cell_parameters_{name}.json", 'w') as json_file:
            json.dump(parameters, json_file, indent=4)

    def knudsen_diffusion(self, specie: str, temperature: float) -> float:
        """
        Calculates Knudsen diffusion of given specie as a function of temperature.
        :param specie: Options are "H2", "O2", "H2O"
        :param temperature: Temperature
        :return: (float) Knudsen diffusion [m^2/s]
        """
        return self.knudsen_factors[specie] * temperature ** 0.5

    def binary_diffusion(self, pair: str, temperature: float, pressure: float) -> float:
        """
        Calculates Binary diffusion of given pair as a function of temperature and pressure.
        :param pair: Options are "H2-H2O" and "O2-N2"
        :param temperature: Temperature [K]
        :param pressure: Total pressure [bar]
        :return: (float) Binary diffusion coefficient [m^2/s]
        """
        D_ij = self.binary_diffusion_factors[pair] * temperature ** 1.75 / pressure
        return D_ij

    def diffusion_polarization_anode(self, current_density: float, temperature: float, pressure: float,
                                     x_H2: float, x_H2O: float) -> float:
        """
        Calculates voltage drop (polarization) due to concentration losses at anode side.
        :param current_density: Current density [A/m^2]
        :param temperature: Temperature [K]
        :param pressure: Total pressure of anode [bar]
        :param x_H2: Molar fraction of H2 in bulk phase
        :param x_H2O: Molar fraction of H2O in bulk phase
        :return: Concentration polarization anode side [Volt]
        """
        RT2F = R*temperature/2/F        # R.T/2.F
        P_H2_bulk = pressure * x_H2     # Partial pressure of H2 in bulk flow
        P_H2O_bulk = pressure * x_H2O   # Partial pressure of H2O in bulk flow

        por_tort = self.diffusion_parameters["porosity"] / self.diffusion_parameters["tortuosity"]    # e/tau
        D_binary = self.binary_diffusion("H2-H2O", temperature, pressure)   # Binary diffusion coefficient D_H2-H2O

        # Effective diffusion coefficients
        D_eff_H2 = por_tort * self.knudsen_diffusion("H2", temperature) * D_binary / (self.knudsen_diffusion("H2", temperature) + D_binary)
        D_eff_H2O = por_tort * self.knudsen_diffusion("H2O", temperature) * D_binary / (self.knudsen_diffusion("H2", temperature) + D_binary)

        # Partial pressures of H2 and H2O at triple phase boundary (TPB)
        P_H2_tpb = P_H2_bulk - (RT2F*1e-5 * current_density * self.dimensions["thickness_anode_electrode"] / D_eff_H2)
        P_H2O_tpb = P_H2O_bulk + (RT2F*1e-5 * current_density * self.dimensions["thickness_anode_electrode"] / D_eff_H2O)

        polarization_anode = RT2F * log(P_H2O_tpb * P_H2_bulk / P_H2O_bulk / P_H2_tpb)
        return polarization_anode

    def diffusion_polarization_cathode(self, current_density: float, temperature: float, pressure: float,
                                       x_O2: float, pure_oxygen: bool = True) -> float:
        """
        Calculates voltage drop (polarization) due to concentration losses at cathode side.
        :param current_density: Current density [A/m^2]
        :param temperature: Temperature [K]
        :param pressure: Total pressure of cathode [bar]
        :param x_O2: Molar fraction of O2 in bulk phase
        :param pure_oxygen: True if pure Oxygen feed, False if Nitrogen is also present in cathode.
        :return: Concentration polarization cathode side [Volt]
        """
        RT2F = R * temperature / 2 / F  # R.T/2.F
        P_O2_bulk = pressure * x_O2  # Partial pressure of O2 in bulk flow
        por_tort = self.diffusion_parameters["porosity"] / self.diffusion_parameters["tortuosity"]  # e/tau

        # If Nitrogen is present in cathode gas, binary diffusion of O2-N2 is also considered.
        if not pure_oxygen:
            D_binary = self.binary_diffusion("O2-N2", temperature, pressure)  # Binary diffusion coefficient D_O2-N2
            # Effective diffusion coefficient
            D_eff_O2 = por_tort * self.knudsen_diffusion("O2", temperature) * D_binary / (self.knudsen_diffusion("O2", temperature) + D_binary)
        else:
            D_eff_O2 = por_tort * self.knudsen_diffusion("O2", temperature)  / (self.knudsen_diffusion("O2", temperature) )

        # Partial pressure of O2 at triple phase boundary (TPB)
        P_O2_tpb = P_O2_bulk - (RT2F*1e-5 * current_density * self.dimensions["thickness_cathode_electrode"] / D_eff_O2)

        polarization_cathode = RT2F * log( P_O2_bulk / P_O2_tpb )
        return polarization_cathode

    def activation_polarization(self, current_density: float, temperature: float) -> float:
        """
        Calculates activation polarization at anode side and cathode side.
        :param current_density: Current density [A/m^2]
        :param temperature: Temperature [K]
        :return: Activation polarization [Volt]
        """
        # Calculate exchange current densities j0
        j0_anode = self.activation_parameters["pre_exp_factor_fuel"]* \
                  exp(-self.activation_parameters["activation_energy_fuel"] / R / temperature)

        j0_cathode = self.activation_parameters["pre_exp_factor_oxygen"] * \
                  exp(-self.activation_parameters["activation_energy_oxygen"] / R / temperature)

        polarization_anode = R * temperature / F * asinh(current_density / j0_anode / 2)
        polarization_cathode = R * temperature / F * asinh(current_density / j0_cathode / 2)

        return polarization_cathode + polarization_anode

    def ohmic_polarization(self, current_density: float, temperature: float) -> float:
        """
        Calculates ohmic polarization. Constant resistance (lumped resistance for multiple parts) is included.
        :param current_density: Current density [A/m^2]
        :param temperature: Temperature [K]
        :return: Ohmic polarization [Volt]
        """
        ionic_conductivity = self.ohmic_parameters["pre_exp_factor_electrolyte"] * \
                             exp(-self.ohmic_parameters["activation_energy_electrolyte"] / R / temperature)

        if self.ohmic_parameters["other_ohmic_resistance"]:
            R_other = self.ohmic_parameters["other_ohmic_resistance"]
        else:
            R_other = 0
        polarization_ohmic = current_density * (self.dimensions["thickness_electrolyte"] / ionic_conductivity + R_other )
        return polarization_ohmic
