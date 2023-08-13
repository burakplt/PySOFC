
# PySOFC: Solid Oxide Fuel Cell Simulations in Python

Main objective of this project is to build an open source tool for modeling Solid Oxide Fuel Cells (SOFC). I wanted to upload the MVP (minimum liable product) for the part  **zero-dimensional models with hydrogen fuel** while I work on other modules.

Code contains lots of comments and explanations. My purpose is to provide a working example for other researchers to build/modify their own code on similar topics.  

## Usage/Examples
A zero-dimensional (lumped parameters) model of a Solid Oxide Fuel Cell can be simulated by writing a script, or using graphical user interface.

### Using graphical user interface
Just write a python script with following lines and run. GUI window will show up.
```python
from pysofc import startGUI

startGUI()
```
![App Screenshot](https://github.com/burakplt/PySOFC/assets/22001926/2c85fe0f-88d8-4094-b7e7-b6f90ca4e5a9)

Watch the video of an example usage case of graphical user interface. It provides much more information than pages long documentation. 


https://github.com/burakplt/PySOFC/assets/22001926/8793145b-03a4-42f9-9ead-f80e24179da9


### Writing a script 
Second option is to write Python script.

```python
from pysofc import HydrogenSOFC                 # Physical representation of an SOFC
from pysofc import LumpedModelHydrogenSOFC      # Zero-Dimensional (Lumped parameters) simulation 

cell_file = r"Parameter Files\cell_parameters_jung2005.json"    # Journal of Power Sources 155 (2006) 145–151

fuel_cell = HydrogenSOFC("SOFC-1", load_from_json= cell_file)
model = LumpedModelHydrogenSOFC(fuel_cell)

model.input["X_fuel_in"] = {"H2": 0.97, "H2O": 0.03}            # Define fuel gas inlet molar composition
model.input["X_oxygen_in"] = {"O2": 0.21, "N2": 0.79}           # Define air/sweep gas inlet molar composition
model.input["molar_flowrate_fuel"] = 0.0001859                  # Fuel gas feed rate [mol/s]
model.input["molar_flowrate_oxygen"] = 0.000495733              # Air/sweep gas feed rate [mol/s]
model.input["T_operation"] = 923                                # Operating temperature [K]
model.input["pressure"] = 1                                     # Operating pressure [bar]
```
For most of the cases these operational inputs would be sufficient. If you want to implement a modified version, please see the all available operation parameters of class  *LumpedModelHydrogenSOFC*  in source file `pysofc/models_zero_dim`.

To generate a **voltage-current (V-I) curve** at given operation conditions you should supply a list of current densities. If related experimental data are present, you can compare simulated V-I curve with experiments.

```python
...
# Define current density points to be evaluated as list with unit [A/cm^2]
current_list = [0.064, 0.13, 0.188, 0.248, 0.375, 0.502, 0.626, 0.753, 0.874, 0.998, 1.125, 1.246, 1.373, 1.496, 1.623, 1.744, 1.868]

# If experimental data are present, give them as two seperate lists for current density and voltage
experiment_current = [0.064, 0.13, 0.188, 0.248, 0.375, 0.502, 0.626, 0.753, 0.874, 0.998, 1.125, 1.246, 1.373, 1.496, 1.623, 1.744, 1.868]
experiment_voltage = [1.017, 0.931, 0.872, 0.826, 0.753, 0.695, 0.644, 0.597, 0.553, 0.517, 0.479, 0.439, 0.404, 0.369, 0.331, 0.296, 0.257]

sim = model.isothermal_voltage_current_curve(J_list=current_list, experiment_J=experiment_current, experiment_V=experiment_voltage, experiment_name='Jung-2005')
```
It returns the following dictionary of results and shows plots by default. If you don't want to have plots to be shown, set function parameter `show_plot=False`.

```python
>>> sim
{"J":Current Density, "V":Voltage, "ASR":Area Specific Resistance, 
 "FU":Fuel Utilisation, "POWER":Electrical Power, "ERROR":Error relative to experiments}
# Each item contains results as list in the same order of current density. 
# ERROR is calculated only if experimental data are provided, otherwise retturns None.
``` 
***Screenshot of generated plots***
![App Screenshot](https://github.com/burakplt/PySOFC/assets/22001926/44ddd6c8-2e18-4546-8eef-aa1c4268d523)

You can also calculate operation voltage for a single current density point or vice versa.
```python
...
model.calculate_operation_voltage(current_density= 0.5)
>>> 0.71975723

model.calculate_current_density(operation_voltage= 0.71975723)
>>> 0.500009
```

In example above all required parameters of a Solide Oxide Cell are read from the file `cell_parameters_jung2005.json`. Please see the structure of a cell parameter file `pysofc/Parameter Files/cell_parameters_jung2005.json`.
Alternatively, these parameters can be defined or altered by changing corresponding variables.

```python
...

fuel_cell = HydrogenSOFC("SOFC-1")                # Note that load_from_json is removed
model = LumpedModelHydrogenSOFC(fuel_cell)

# Define or alter parameters one by one
fuel_cell.dimensions["cell_width"] = 0.005

# Or replace with dictionary as whole. 
# If it is not relevant for the simulation, you can define it as None.
fuel_cell.dimensions = {
            "channel_height": None,
            "channel_width": 0.0005,
            "channel_length": 0.01,
            "total_electroactive_area": 0.0016,
            "thickness_cathode_electrode": 30e-6,
            "thickness_anode_electrode": 20e-6,
            "thickness_electrolyte": 8e-6,
        }

...
```
Please see the exact structure of parameters of *HydrogenSOFC* class in `pysofc/solid_oxide_cell.py`. Using GUI program makes this process easier, and is the recommended method to create/modify parameter .json files.

Functions in `pysofc/utility.py` module can be used stand-alone (static methods). 

```python
from pysofc import BasicThermo
from pysofc import PolarizationFunctions

# Available species: "H2", "H2O", "O2", "CO", "CO2", "CH4", "N2"

BasicThermo.enthalpy_gas("H2", T= 600)   # Calculates the gas phase enthalpy of Hydrogen at 600 K
>>> 8804.468                            # J/mol

BasicThermo.entropy_gas("H2O", T= 600)   # Calculates the gas phase entropy of Hydrogen at 600 K
>>> 151.086                              # J/mol.K 

# Binary diffusion coefficient [m^2/s] of Water-Hydrogen at 800 K and 1 bar
PolarizationFunctions.binary_diffusion_explicit(specie_1="H2", specie_2="H2O", temperature=800, pressure=1)
>>> 0.0005156019
```
See `pysofc/utility.py` for other possible usages.  
## Roadmap

Ongiong works:
- Simulation of SOFC running on reformate fuels containing Methane, Carbon monoxide, Carbon dioxode 
- One-Dimensional SOFC modeling both hydrogen and reformate fuels.
- Integration of surface reaction kinetics for reformate fuels using [Cantera](https://github.com/Cantera/cantera)
- New GUI that covers other chemical components and models

Possible improvements:
- Implementing a One-Dimensional model using [SIMPLE](https://en.wikipedia.org/wiki/SIMPLE_algorithm) algorithm.
- Publishing all as Python package

