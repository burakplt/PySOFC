#!/usr/bin/env python
# -*- coding: utf-8 -*-

from PyQt5.QtWidgets import (QMainWindow, QApplication, QFileDialog, QTableWidgetItem)
from .GUI.nullDimWindow import Ui_NullDimWindow
from .GUI.resultsWindow import Ui_ResultWindow
from .GUI.parametersWindow import Ui_ParameterManager
from pathlib import Path
from .solid_oxide_cell import HydrogenSOFC
from .models_zero_dim import LumpedModelHydrogenSOFC
import numpy as np
import matplotlib.pyplot as plt
import json
from datetime import datetime


class NullDimWindow(QMainWindow, Ui_NullDimWindow):

    def __init__(self, *args, **kwargs):
        super(NullDimWindow, self).__init__(*args, **kwargs)
        self.setupUi(self)
        self.rsoc = None        # HydrogenSOFC object
        self.model = None       # LumpedHydrogen Model object

        self.browseFile.clicked.connect(self.file_browser)                  # Browse file button
        self.generateCurve.clicked.connect(self.generate_VI_curve)          # Generate V-I curve button
        self.cellManager.clicked.connect(self.cell_manager)                 # Cell Manager button
        self.showCurrentCell.clicked.connect(self.show_cell_parameters)     # Show Current Parameters button

        self.res_window = None              # Simulation results window
        self.parameter_window = None        # Cell parameter manager window
        self.sim_param_window = None        # Cell parameters used in simulation window

        self.cell_file = None   # Path to json parameter cell file
        self.J_list, self.mean_temperature = None, None
        self.experimental_data = {"name":None, "J":None, "V":None}
        self.is_input_valid = [False, False]    # [Essential Input , Current-Temperature Input]
        self.show()

    def timestamp(self):
        # Return current timestamp
        t = datetime.now()
        return t.strftime("%H:%M:%S")

    def file_browser(self):
        """
        File browser to select json cell parameter file.
        Assigns: self.cell_file, self.rsoc, self.model, self.parameter_window, self.sim_param_window, self.res_window
        :return: None
        """
        home_dir = str(Path(__file__).parent)
        pth = QFileDialog.getOpenFileName(self, 'Open file', home_dir)
        if pth[0]:
            self.cell_file = pth[0]
            self.labelFileName.setText(f"Cell File: {pth[0].split('/')[-1]}")
            try:
                self.rsoc = HydrogenSOFC("SOFC", load_from_json=self.cell_file)
                self.model = LumpedModelHydrogenSOFC(self.rsoc)
                self.parameter_window = ParameterWindow(self)
                self.sim_param_window = SimParameterWindow(self)
                self.res_window = ResultWindow(self)

                self.logText.append(f"{self.timestamp()} [OK] Parameters from {pth[0].split('/')[-1]} loaded.\n")
            except:
                self.logText.append(f"{self.timestamp()} [ERROR] Parameters from {pth[0].split('/')[-1]} could not be loaded. Check file!\n")

    def get_essential_input(self):
        """
        Read essential input: Pressure, Molar flow rates, Molar fractions
        :return: None
        """
        try:
            self.model.input["pressure"] = float(self.input_pressure.text())
            self.model.input["molar_flowrate_fuel"] = float(self.input_molar_flowrate_fuel.text())
            self.model.input["molar_flowrate_oxygen"] = float(self.input_molar_flowrate_oxygen.text())
            self.model.input["X_fuel_in"] = {"H2":float(self.input_x_H2.text()), "H2O":float(self.input_x_H2O.text())}
            x_N2 = float(self.input_x_N2.text()) if self.input_x_N2.text() else 0
            self.model.input["X_oxygen_in"] = {"O2": float(self.input_x_O2.text()), "N2":x_N2}
            self.is_input_valid[0] = True
        except:
            self.logText.append(f"{self.timestamp()} [ERROR] Check essential parameters!\n")
            self.is_input_valid[0] = False

    def get_voltage_current_input(self):
        """
        Read current density list, temperature, experimental J and V lists
        :return: None
        """
        try:
            self.mean_temperature = float(self.input_meanTemperature.text())
            J_list = self.input_jList.text()
            if len(J_list.split(":")) == 3:
                bounds = J_list.split(":")
                self.J_list = np.linspace(float(bounds[0]), float(bounds[1]), int(bounds[2])) * 1e+4
            else:
                self.J_list = np.array(J_list.split(","), dtype=float) * 1e+4
            self.logText.append(f"{self.timestamp()} [OK] Parameters for V-I curve are valid.\n")
            self.is_input_valid[1] = True

            try:
                self.experimental_data["name"] = self.input_experimentName.text()
                self.experimental_data["J"] = np.array(self.input_jExperiment.text().split(","), dtype=float)
                self.experimental_data["V"] = np.array(self.input_vExperiment.text().split(","), dtype=float)
                self.logText.append(f"{self.timestamp()} [OK] Input for experimental data are valid.\n")
                if self.experimental_data["J"].size != self.experimental_data["V"].size:
                    self.logText.append(f"{self.timestamp()} [WARNING] Length of current density and voltage not matching!\n")
                    raise Exception
            except:
                self.logText.append(
                    f"{self.timestamp()} [WARNING] Check experimental data input! Comparison with experiment is ignored.\n")
        except:
            self.logText.append(f"{self.timestamp()} [ERROR] Check current density / temperature input!\n")
            self.is_input_valid[1] = False

    def generate_VI_curve(self):
        """
        Generate Voltage-Current curve and plots etc.
        :return: None
        """
        try:
            self.get_essential_input()
            self.get_voltage_current_input()
            if sum(self.is_input_valid) < 2:
                raise Exception

            self.model.input["T_operation"] = self.mean_temperature
            sim = self.model.isothermal_voltage_current_curve(self.J_list, list(self.experimental_data["J"]),
                  list(self.experimental_data["V"]), self.input_experimentName.text(), show_plot=False)

            self.logText.append(f"{self.timestamp()} [OK] Inputs to generate V_I curve are valid.\n")
            self.res_window.write_results_table(sim)
            self.res_window.write_parameters_table()
            self.res_window.show()
            plt.show()
        except:
            self.logText.append(f"{self.timestamp()} [ERROR] Input error, see messages above.\n")

    def cell_manager(self):
        """
        Handles cell manager button click. Connects to **ParameterWindow**.
        :return: None
        """
        self.parameter_window.show()

    def show_cell_parameters(self):
        # Handles show cell parameters button click. Connects to ParameterWindow
        is_parameter_window = self.parameter_window.show_current_parameters()
        if is_parameter_window:
            self.parameter_window.show()
        else:
            self.logText.append(f"{self.timestamp()} [ERROR] No cell file is selected!\n")


class ResultWindow(QMainWindow, Ui_ResultWindow):

    def __init__(self, main_window:NullDimWindow, *args, **kwargs):
        super(ResultWindow, self).__init__(*args, **kwargs)
        self.setupUi(self)
        self.main_window = main_window  # Main Window NullDimWindow
        self.showCellParameters.clicked.connect(self.handle_show_parameters)    # Show parameters used in simulation

    def write_results_table(self, result:dict):
        """
        Creates table of simulation output: self.tableResult table on the right column
        :param result: {"J", "V", "FU", "POWER", "ERROR"}
        :return: None
        """
        J, V, FU, POWER, ERROR = result["J"], result["V"], result["FU"], result["POWER"], result["ERROR"]
        header = ["Current Density (A/cm^2)", "Voltage (V)", "Fuel Utilisation", "Power Density (W/cm^2)", "Relative Error %"]
        col_size = 5 if ERROR is not None else 4

        rows = len(J)
        self.tableResult.setRowCount(rows)
        self.tableResult.setColumnCount(col_size)
        self.tableResult.setHorizontalHeaderLabels(header)
        for row in range(rows):
            self.tableResult.setItem(row, 0, QTableWidgetItem(f"{J[row]:.4f}"))
            self.tableResult.setItem(row, 1, QTableWidgetItem(f"{V[row]:.4f}"))
            self.tableResult.setItem(row, 2, QTableWidgetItem(f"{FU[row]:.4f}"))
            self.tableResult.setItem(row, 3, QTableWidgetItem(f"{POWER[row]:.4f}"))
            if col_size == 5:
                self.tableResult.setItem(row, 4, QTableWidgetItem(f"{ERROR[row]:.4f}"))

        self.tableResult.resizeColumnsToContents()

    def write_parameters_table(self):
        """
        Creates table of simulation operational parameters: self.tableParameters table on the left column
        :return: None
        """
        model = self.main_window.model
        self.tableParameters.setRowCount(7)
        self.tableParameters.setColumnCount(2)
        self.tableParameters.setHorizontalHeaderLabels(["Parameter", "Value"])

        self.tableParameters.setItem(0, 0, QTableWidgetItem( "Cell Name" ) )
        self.tableParameters.setItem(0, 1, QTableWidgetItem(str(model.FC.cell_id)))
        self.tableParameters.setItem(1, 0, QTableWidgetItem("Temperature (K)"))
        self.tableParameters.setItem(1, 1, QTableWidgetItem(str(model.T_op)))
        self.tableParameters.setItem(2, 0, QTableWidgetItem("Pressure (bar)"))
        self.tableParameters.setItem(2, 1, QTableWidgetItem(str(model.input["pressure"])))
        self.tableParameters.setItem(3, 0, QTableWidgetItem("Fuel flow rate (mol/s)"))
        self.tableParameters.setItem(3, 1, QTableWidgetItem(str(model.input["molar_flowrate_fuel"])))
        self.tableParameters.setItem(4, 0, QTableWidgetItem("Air flow rate (mol/s)"))
        self.tableParameters.setItem(4, 1, QTableWidgetItem(str(model.input["molar_flowrate_oxygen"])))
        self.tableParameters.setItem(5, 0, QTableWidgetItem("Fuel molar composition"))
        self.tableParameters.setItem(5, 1, QTableWidgetItem(
            f"H2: {model.input['X_fuel_in']['H2']}, H2O: {model.input['X_fuel_in']['H2O']}"))
        self.tableParameters.setItem(6, 0, QTableWidgetItem("Air molar composition"))
        self.tableParameters.setItem(6, 1, QTableWidgetItem(
            f"O2: {model.input['X_oxygen_in']['O2']}, N2: {model.input['X_oxygen_in']['N2']}"))

        self.tableParameters.resizeColumnsToContents()

    def handle_show_parameters(self):
        """
        Handles show parameters button click. Connects to **SimParameterWindow**.
        :return: None
        """
        win = self.main_window.sim_param_window
        win.parameter_view()
        win.show()


class ParameterWindow(QMainWindow, Ui_ParameterManager):

    def __init__(self, main_window:NullDimWindow, *args, **kwargs):
        super(ParameterWindow, self).__init__(*args, **kwargs)
        self.setupUi(self)
        self.main_window = main_window          # Main window object NullDimWindow
        self.cell_object = main_window.rsoc     # HydrogenSOFC object
        self.applyChanges.hide()                # Hide apply changes button

        self.dimensions = {
            "channel_height": self.channel_height,
            "channel_width": self.channel_width,
            "channel_length": self.channel_length,
            "total_electroactive_area": self.total_electroactive_area,
            "thickness_cathode_electrode": self.thickness_cathode_electrode,
            "thickness_anode_electrode": self.thickness_anode_electrode,
            "thickness_electrolyte": self.thickness_electrolyte,
        }
        self.activation_parameters = {
            "activation_energy_fuel": self.activation_energy_fuel,
            "activation_energy_oxygen": self.activation_energy_oxygen,
            "pre_exp_factor_fuel": self.pre_exp_factor_fuel,
            "pre_exp_factor_oxygen": self.pre_exp_factor_oxygen
        }
        self.diffusion_parameters = {
            "porosity": self.porosity,
            "tortuosity": self.tortuosity,
            "pore_diameter": self.pore_diameter,
        }
        self.ohmic_parameters = {
            "activation_energy_electrolyte": self.activation_energy_electrolyte,
            "pre_exp_factor_electrolyte": self.pre_exp_factor_electrolyte,
            "other_ohmic_resistance": self.other_ohmic_resistance
        }

        self.saveFile.clicked.connect(self.handle_save)                 # Save as file button
        self.loadFile.clicked.connect(self.handle_load)                 # Load from file button
        self.applyChanges.clicked.connect(self.handle_apply_changes)    # Apply changes and return simulation button

    def read_input(self):
        # Reads input text boxes and convert to numeric values
        dimension, activation, diffusion, ohmic = {}, {}, {}, {}
        for key, value in self.dimensions.items():
            try:
                dimension[key] = float(value.text())
            except:
                dimension[key] = None

        for key, value in self.activation_parameters.items():
            try:
                activation[key] = float(value.text())
            except:
                dimension[key] = None

        for key, value in self.diffusion_parameters.items():
            try:
                diffusion[key] = float(value.text())
            except:
                dimension[key] = None

        for key, value in self.ohmic_parameters.items():
            try:
                ohmic[key] = float(value.text())
            except:
                dimension[key] = None
        return dimension, activation, diffusion, ohmic

    def write_input_from_json(self, json_data):
        """
        Populate input text boxes with date read from json file
        :param json_data: {"cell_id", "dimensions", "activation_parameters", "diffusion_parameters" "ohmic_parameters"}
        :return: None
        """
        for key, value in self.dimensions.items():
            value.setText(str(json_data["dimensions"][key]))

        for key, value in self.activation_parameters.items():
            value.setText(str(json_data["activation_parameters"][key]))

        for key, value in self.diffusion_parameters.items():
            value.setText(str(json_data["diffusion_parameters"][key]))

        for key, value in self.ohmic_parameters.items():
            value.setText(str(json_data["ohmic_parameters"][key]))

    def handle_save(self):
        # Save input data as cell parameters json file. File name read from text box 'Label'
        name = self.cell_id.text()
        if name:
            input_data = self.read_input()
            parameters = {"cell_id": name, "dimensions": input_data[0], "activation_parameters": input_data[1],
                          "diffusion_parameters": input_data[2], "ohmic_parameters": input_data[3]}
            parameters["diffusion_parameters"]["diffusion_volumes"] = {"H2": 6.12, "N2": 18.5, "O2": 16.3, "H2O": 13.1}
            with open(f"cell_parameters_{name}.json", 'w') as json_file:
                json.dump(parameters, json_file, indent=4)

            self.textBrowser.append(f"Parameters were saved as file: cell_parameters_{name}.json\n")
        else:
            self.textBrowser.append("Name/ID is not given! Please enter a cell identifier.\n")

    def handle_load(self):
        # Handle load from file button
        home_dir = str(Path(__file__).parent)
        pth = QFileDialog.getOpenFileName(self, 'Open file', home_dir)
        if pth[0]:
            with open(pth[0], 'r') as f:
                data = json.load(f)
                self.write_input_from_json(data)
        self.textBrowser.append(f"[OK] Parameters were loaded from: {pth[0]}\n")

    def show_current_parameters(self) -> bool:
        """
        Handle show current parameters button on the **Main Window**.
        :return: True if a cell file is selected, False if not.
        """
        if self.cell_object:
            fc = self.cell_object
            data = {"dimensions":fc.dimensions, "activation_parameters": fc.activation_parameters,
                              "diffusion_parameters": fc.diffusion_parameters, "ohmic_parameters": fc.ohmic_parameters }

            self.write_input_from_json(data)
            self.cell_id.setText(self.cell_object.cell_id)
            self.applyChanges.show()
            return True
        else:
            return False

    def handle_apply_changes(self):
        """
        Handle apply changes button click. Read all input data > Update HydrogenSOFC > Recalculate diffusion parameters
        :return: None
        """
        fc = self.cell_object
        input_data = self.read_input()
        input_data[2]["diffusion_volumes"] = {"H2": 6.12, "N2": 18.5, "O2": 16.3, "H2O": 13.1}
        fc.dimensions = input_data[0]
        fc.activation_parameters = input_data[1]
        fc.diffusion_parameters = input_data[2]
        fc.ohmic_parameters = input_data[3]
        fc.pre_calculations()
        self.close()


class SimParameterWindow(QMainWindow, Ui_ParameterManager):
    def __init__(self, main_window, *args, **kwargs):
        super(SimParameterWindow, self).__init__(*args, **kwargs)
        self.setupUi(self)
        self.main_window = main_window              # Main window object NullDimWindow
        self.cell_object = self.main_window.rsoc    # HydrogenSOFC object

        self.dimensions = {
            "channel_height": self.channel_height,
            "channel_width": self.channel_width,
            "channel_length": self.channel_length,
            "total_electroactive_area": self.total_electroactive_area,
            "thickness_cathode_electrode": self.thickness_cathode_electrode,
            "thickness_anode_electrode": self.thickness_anode_electrode,
            "thickness_electrolyte": self.thickness_electrolyte,
        }
        self.activation_parameters = {
            "activation_energy_fuel": self.activation_energy_fuel,
            "activation_energy_oxygen": self.activation_energy_oxygen,
            "pre_exp_factor_fuel": self.pre_exp_factor_fuel,
            "pre_exp_factor_oxygen": self.pre_exp_factor_oxygen
        }
        self.diffusion_parameters = {
            "porosity": self.porosity,
            "tortuosity": self.tortuosity,
            "pore_diameter": self.pore_diameter,
        }
        self.ohmic_parameters = {
            "activation_energy_electrolyte": self.activation_energy_electrolyte,
            "pre_exp_factor_electrolyte": self.pre_exp_factor_electrolyte,
            "other_ohmic_resistance": self.other_ohmic_resistance
        }

    def write_input_from_json(self, json_data):
        """
        Populate input text boxes with date read from json file
        :param json_data: {"cell_id", "dimensions", "activation_parameters", "diffusion_parameters" "ohmic_parameters"}
        :return: None
        """
        for key, value in self.dimensions.items():
            value.setText(str(json_data["dimensions"][key]))

        for key, value in self.activation_parameters.items():
            value.setText(str(json_data["activation_parameters"][key]))

        for key, value in self.diffusion_parameters.items():
            value.setText(str(json_data["diffusion_parameters"][key]))

        for key, value in self.ohmic_parameters.items():
            value.setText(str(json_data["ohmic_parameters"][key]))

    def parameter_view(self):
        """
        Handle show parameters button click on **ResultWindow** to show parameters used in simulation
        :return: None
        """
        fc = self.cell_object
        data = {"dimensions": fc.dimensions, "activation_parameters": fc.activation_parameters,
                "diffusion_parameters": fc.diffusion_parameters, "ohmic_parameters": fc.ohmic_parameters}

        self.write_input_from_json(data)
        self.setWindowTitle(f"Cell Parameters {self.cell_object.cell_id}")
        self.rightBar.hide()
        self.setFixedSize(680, 610)

def start_GUI():
    import sys
    app = QApplication(sys.argv)
    window = NullDimWindow()
    window.setWindowTitle("Lumped Hydrogen rSOC")
    sys.exit(app.exec_())

if __name__ == "__main__":
    import sys
    app = QApplication(sys.argv)
    window = NullDimWindow()
    window.setWindowTitle("Lumped Hydrogen rSOC")
    sys.exit(app.exec_())