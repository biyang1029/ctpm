# Patient Temperature Simulation Project

## Overview
This project is focused on simulating the temperature distribution within a patient's body during various conditions, such as perioperative scenarios. The main goal is to model the dynamics of body temperature in a realistic way, taking into account various influencing factors such as convective heat transfer, radiation, and metabolic heat generation. The simulation uses a finite difference approach to calculate temperature distributions within a three-dimensional grid that represents the patient.

The project is implemented in C++ and makes use of classes and functions to represent patients, initialize simulation parameters, calculate boundary conditions, and iteratively solve the temperature field.

## Features
- **Three-dimensional Temperature Distribution**: The temperature field is represented by a 3D matrix that models the heat transfer within the patient's body.
- **Boundary Condition Calculations**: Includes complex boundary conditions like convective heat transfer coefficients and radiative heat losses, which vary depending on different regions of the body.
- **Iterative Heat Transfer Calculation**: Simulates the heat transfer process iteratively until convergence is achieved, providing an accurate simulation of real-world conditions.
- **Variable Thermal Conductivity**: Adjusts the thermal conductivity dynamically depending on core temperature conditions to simulate various physiological responses.
- **File Output**: Outputs temperature data to a file for further analysis.

## Files
### 1. `main.cpp`
The entry point of the program. Initializes the patient and simulation data and starts the simulation.

### 2. `Initialization.h` & `Initialization.cpp`
Contains the `InitializationData` struct and the `initializeSimulation` function. These handle setting up the grid, initial conditions, and other simulation parameters.

### 3. `PatientData.h`
Defines the `Patient` structure that holds basic patient information like height, weight, core temperature, and environmental temperature.

### 4. `Simulation.h` & `Simulation.cpp`
Implements the main simulation function (`runSimulation`). Handles the temperature calculations using finite difference techniques and calculates boundary conditions for each node.

## Compilation Instructions
This project uses C++ and requires a modern compiler that supports the C++11 standard or newer. You can compile the project using the following commands in your terminal:

```sh
# Using g++ compiler
g++ main.cpp Initialization.cpp Simulation.cpp -o patient_simulation
```
Alternatively, if you are using Visual Studio, you can open the solution file and build the project directly from the IDE.

## How to Run
After compiling the code, you can run the executable as follows:

```sh
./patient_simulation
```
Upon execution, the simulation data, including the temperature field at each grid point, will be saved in a file named `temperature_data.txt` (or `temperature_data.csv` if you use the CSV output method).

## Usage
The project can be used for various purposes, including:
1. **Medical Training**: To visualize and understand body temperature distribution during surgeries or different clinical conditions.
2. **Research**: To validate new thermal models for patients in operating rooms or develop new temperature regulation techniques.
3. **Simulation-Based Testing**: Testing how different factors, such as room temperature or changes in convective heat transfer, can impact core and skin temperatures.

## Output
The simulation produces output in the form of a text or CSV file that contains temperature values for each grid point, indexed by their coordinates `(i, j, k)`. You can analyze this data further using tools like Excel or Python for visualization and post-processing.

## Key Variables
- `temperatureField`: A 3D vector storing the temperature at each grid point.
- `ConvectiveCoefficient`: Varies by the different parts of the body, used in calculating convective heat transfer.
- `BoundaryConditions`: Struct that holds the boundary condition parameters for each grid point.
- `resi`: Residual used in the convergence check to determine if the temperature field has reached stability.

## Known Issues
1. **NaN Values in Temperature**: Sometimes during the iteration, `NaN` (Not a Number) values can appear due to issues such as division by zero or incorrect initial conditions. Double-check the boundary conditions and ensure proper initialization of all variables.
2. **Infinite Loop in Iterations**: The main iteration loop may not converge in some cases, leading to an infinite loop. Make sure the convergence criteria are properly set and that numerical instability is minimized by using appropriate time step values (`DeltaTime`).

## Future Work
- **Refine Thermal Conductivity Modeling**: Implement a more accurate model for dynamically adjusting thermal conductivity based on physiological feedback.
- **Add More Patient Variability**: Introduce more diverse patient models, including age, BMI, or health condition-related adjustments.
- **Parallel Computation**: The iteration loop could be parallelized using OpenMP to speed up the simulation for large grid sizes.


## Contact
For more information, questions, or suggestions, please contact the developer at:

- Name: Yang Bi
- Email: yang.bi@sjtu.edu.cn
