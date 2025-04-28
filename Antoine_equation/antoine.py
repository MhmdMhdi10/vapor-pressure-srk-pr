import numpy as np


antoine_constants = {
    "C1": {"A": 3.9895, "B": 443.028, "C": -0.49, "T_c": 190.6, "P_c": 46.1, "omega": 0.011},
    "C3": {"A": 4.01158, "B": 834.26, "C": -22.763, "T_c": 369.9, "P_c": 42.5, "omega": 0.152},
    "C5": {"A": 3.9892, "B": 1070.617, "C": -40.454, "T_c": 469.8, "P_c": 33.6, "omega": 0.249},
    "C7": {"A": 4.81803, "B": 1635.409, "C": -27.338, "T_c": 540, "P_c": 27.4, "omega": 0.302},
    "H2O": {"A": 4.6543, "B": 1435.264, "C": -64.848, "T_c": 647, "P_c": 220.64, "omega": 0.344},
}


def antoine_equation(A, B, C, T):
    """
    Calculate the vapor pressure using Antoine's equation.

    Parameters:
        A (float): Antoine constant A
        B (float): Antoine constant B
        C (float): Antoine constant C
        T (float): Temperature in Kelvin

    Returns:
        float: Vapor pressure in the unit consistent with Antoine constants
    """
    # Ensure T > 0 to avoid division by zero or log of negative
    if T <= 0:
        raise ValueError("Temperature must be greater than 0 K")
    # Antoine equation
    P_sat = 10 ** (A - (B / (T + C)))
    return P_sat


# # Test the function for methane (C1) in the range 70K to 190K
# component = "C1"
# T_range = np.linspace(70, antoine_constants[component]["T_c"], 70)
#
# # Calculate vapor pressures
# P_sat_values = [antoine_equation(antoine_constants[component]["A"],
#                                  antoine_constants[component]["B"],
#                                  antoine_constants[component]["C"],
#                                  T) for T in T_range]
#
# # Print results for verification
# for T, P_sat in zip(T_range, P_sat_values):
#     print(f"T = {T:.2f} K, P_sat = {P_sat:.4f} (consistent with Antoine units)")
#
