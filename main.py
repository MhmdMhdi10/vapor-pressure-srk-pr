import numpy as np
import matplotlib.pyplot as plt
from Antoine_equation.antoine import antoine_equation
from srk_eos.srk_eos import srk_eos, srk_fugacity_coefficient
from newthon_raphson.newton_raphson import newton_raphson
from Utility.fv_calculator import fv_calculator
from Utility.fl_calculator import fl_calculator
from peng_robinson.peng_robinson import peng_robinson, peng_robinson_fugacity_coefficient

# J·mol^−1·K^−1
R = 8.3144598

antoine_constants = {
    "C1": {"A": 3.9895, "B": 443.028, "C": -0.49, "T_c": 190.6, "P_c": 46.1, "omega": 0.011},
    "C3": {"A": 4.01158, "B": 834.26, "C": -22.763, "T_c": 369.9, "P_c": 42.5, "omega": 0.152},
    "C5": {"A": 3.9892, "B": 1070.617, "C": -40.454, "T_c": 469.8, "P_c": 33.6, "omega": 0.249},
    "C7": {"A": 4.81803, "B": 1635.409, "C": -27.338, "T_c": 540, "P_c": 27.4, "omega": 0.302},
    "H2O": {"A": 4.6543, "B": 1435.264, "C": -64.848, "T_c": 647, "P_c": 220.64, "omega": 0.344},
}

def calculate_vapor_pressure(component, T_min, T_max):
    temperatures_f = []
    vapor_pressures_psi = []

    for t in range(T_min, T_max):
        # Calculate SRK parameters and vapor pressure
        P_sat = antoine_equation(antoine_constants[component]["A"],
                                 antoine_constants[component]["B"],
                                 antoine_constants[component]["C"],
                                 t) * 1e5  # Convert to Pascals

        Z_factors, A, B = srk_eos(t, P_sat, antoine_constants[component])

        # Vapor and saturated liquid Z factors
        Z_vapor = max(Z_factors)  # Choose the root closest to 1 for vapor phase
        Z_sat = max(Z_factors)  # Choose the root closest to 0 for liquid phase (saturated)

        # Fugacity coefficients for vapor and saturated liquid
        phi_vapor = srk_fugacity_coefficient(Z_vapor, A, B)
        phi_sat = srk_fugacity_coefficient(Z_sat, A, B)

        fv = fv_calculator(phi_vapor, P_sat)
        fl = fl_calculator(Z_sat, P_sat, P_sat, phi_sat, R, t)

        # Root finding for SRK fugacity correction
        root_P = (newton_raphson(srk_eos, srk_fugacity_coefficient, R, P_sat, t, antoine_constants[component], fv, fl)
                  * 1e-5)  # Convert to bar

        # Convert bar to psi
        vapor_pressure_psi = root_P * 14.5038

        # Convert Kelvin to Fahrenheit
        temperature_f = (t - 273.15) * 9 / 5 + 32

        # Append converted values
        temperatures_f.append(temperature_f)
        vapor_pressures_psi.append(vapor_pressure_psi)

    return temperatures_f, vapor_pressures_psi


def calculate_vapor_pressure1(component, T_min, T_max):
    temperatures_f = []
    vapor_pressures_psi = []

    for t in range(T_min, T_max):
        # Calculate SRK parameters and vapor pressure
        P_sat = antoine_equation(antoine_constants[component]["A"],
                                 antoine_constants[component]["B"],
                                 antoine_constants[component]["C"],
                                 t) * 1e5  # Convert to Pascals

        Z_factors = peng_robinson(t, P_sat, antoine_constants[component])
        print(Z_factors)

        # Vapor and saturated liquid Z factors
        Z_vapor = Z_factors['Z_vapor']  # Choose the root closest to 1 for vapor phase
        Z_sat = Z_factors['Z_vapor']  # Choose the root closest to 0 for liquid phase (saturated)

        # Fugacity coefficients for vapor and saturated liquid
        phi_vapor = peng_robinson_fugacity_coefficient(t, P_sat, Z_vapor, antoine_constants[component])
        phi_sat = peng_robinson_fugacity_coefficient(t, P_sat, Z_sat, antoine_constants[component])

        fv = fv_calculator(phi_vapor, P_sat)
        fl = fl_calculator(Z_sat, P_sat, P_sat, phi_sat, R, t)

        # Root finding for SRK fugacity correction
        root_P = (newton_raphson(srk_eos, srk_fugacity_coefficient, R, P_sat, t, antoine_constants[component], fv, fl)
                  * 1e-5)  # Convert to bar

        # Convert bar to psi
        vapor_pressure_psi = root_P * 14.5038

        # Convert Kelvin to Fahrenheit
        temperature_f = (t - 273.15) * 9 / 5 + 32

        # Append converted values
        temperatures_f.append(temperature_f)
        vapor_pressures_psi.append(vapor_pressure_psi)

    return temperatures_f, vapor_pressures_psi


number = int(input("enter 1 for srk and 2 for peng robin-robinson: "))

if number == 1:

    # Plot all components
    plt.figure(figsize=(10, 8))
    for component in antoine_constants:
        T_min = 70  # Starting temperature in Kelvin
        T_max = int(antoine_constants[component]["T_c"])  # Critical temperature in Kelvin

        # Calculate vapor pressures for the component
        temperatures_f, vapor_pressures_psi = calculate_vapor_pressure(component, T_min, T_max)

        # Plot each component
        plt.plot(temperatures_f, vapor_pressures_psi, label=component)

    # Add labels, legend, and title
    plt.xlabel("Temperature (°F)", fontsize=12)
    plt.ylabel("Vapor Pressure (psi)", fontsize=12)
    plt.title("Vapor Pressure vs. Temperature for All Components SRK", fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend(fontsize=10)
    plt.show()


elif number == 2:

    # Plot all components
    plt.figure(figsize=(10, 8))
    for component in antoine_constants:
        T_min = 70  # Starting temperature in Kelvin
        T_max = int(antoine_constants[component]["T_c"])  # Critical temperature in Kelvin

        # Calculate vapor pressures for the component
        temperatures_f, vapor_pressures_psi = calculate_vapor_pressure1(component, T_min, T_max)

        # Plot each component
        plt.plot(temperatures_f, vapor_pressures_psi, label=component)

    # Add labels, legend, and title
    plt.xlabel("Temperature (°F)", fontsize=12)
    plt.ylabel("Vapor Pressure (psi)", fontsize=12)
    plt.title("Vapor Pressure vs. Temperature for All Components Peng-Robinson", fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend(fontsize=10)
    plt.show()

else:
    raise ValueError("choose 1 or 2")
