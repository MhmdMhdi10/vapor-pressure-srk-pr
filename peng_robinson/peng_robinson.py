import numpy as np

R = 0.08314

antoine_constants = {
    "C1": {"A": 3.9895, "B": 443.028, "C": -0.49, "T_c": 190.6, "P_c": 46.1, "omega": 0.011},
    "C3": {"A": 4.01158, "B": 834.26, "C": -22.763, "T_c": 369.9, "P_c": 42.5, "omega": 0.152},
    "C5": {"A": 3.9892, "B": 1070.617, "C": -40.454, "T_c": 469.8, "P_c": 33.6, "omega": 0.249},
    "C7": {"A": 4.81803, "B": 1635.409, "C": -27.338, "T_c": 540, "P_c": 27.4, "omega": 0.302},
    "H2O": {"A": 4.6543, "B": 1435.264, "C": -64.848, "T_c": 647, "P_c": 220.64, "omega": 0.344},
}


def peng_robinson(T, P, component):
    """
    Calculate compressibility factor (Z) using the Peng-Robinson EOS.

    Parameters:
        R (float): Universal gas constant (J/mol.K)
        P (float): Pressure (bar)
        T (float): Temperature (K)
        Tc (float): Critical temperature (K)
        Pc (float): Critical pressure (bar)
        omega (float): Acentric factor

    Returns:
        dict: Compressibility factors (Z_vapor, Z_liquid), or None if no roots exist.
    """

    Tc = component["T_c"]
    Pc = component["P_c"] * 1e5  # Convert P_c from atm to Pa
    omega = component["omega"]

    # Convert P and Pc to Pa
    P = P * 1e5
    Pc = Pc * 1e5

    # Reduced temperature and alpha term
    Tr = T / Tc
    alpha = (1 + (0.37464 + 1.54226 * omega - 0.26992 * omega ** 2) * (1 - np.sqrt(Tr))) ** 2

    # Calculate 'a' and 'b'
    a = 0.45724 * (R * Tc) ** 2 / Pc * alpha
    b = 0.07780 * R * Tc / Pc

    # A and B terms for the cubic equation
    A = a * P / (R ** 2 * T ** 2)
    B = b * P / (R * T)

    # Define the cubic equation coefficients
    coeffs = [1, -(1 - B), A - (3 * B ** 2) - (2 * B), -(A * B - B ** 2 - B ** 3)]

    # Solve the cubic equation
    Z_roots = np.roots(coeffs)

    # Filter real roots
    Z_roots_real = Z_roots[np.isreal(Z_roots)].real

    # Sort roots in ascending order
    Z_roots_real.sort()

    # Return vapor and liquid phase Z-factors
    if len(Z_roots_real) == 1:  # Single root (e.g., supercritical state)
        return {"Z_vapor": float(Z_roots_real[0]), "Z_liquid": None}
    elif len(Z_roots_real) > 1:  # Two or three roots
        return {"Z_vapor": float(Z_roots_real[-1]), "Z_liquid": float(Z_roots_real[0])}
    else:  # No real roots (should not happen in practical conditions)
        return {"Z_vapor": None, "Z_liquid": None}


def peng_robinson_fugacity_coefficient(T, P, Z, component):
    """
    Calculate the fugacity coefficient for the vapor phase using Peng-Robinson EOS.

    Parameters:
        T (float): Temperature in Kelvin.
        P (float): Pressure in Pa.
        Z (float): Compressibility factor.
        component (dict): Dictionary containing critical properties (T_c, P_c, omega) for the component.

    Returns:
        float: Fugacity coefficient (phi) for the vapor phase.
    """
    # Retrieve critical properties
    Tc = component["T_c"]
    Pc = component["P_c"] * 1e5  # Convert P_c from atm to Pa
    omega = component["omega"]

    # Calculate a and b parameters
    Tr = T / Tc
    alpha = (1 + (0.37464 + 1.54226 * omega - 0.26992 * omega ** 2) * (1 - np.sqrt(Tr))) ** 2
    a = 0.45724 * (R * Tc) ** 2 * alpha / Pc
    b = 0.07780 * R * Tc / Pc

    # Calculate A and B
    A = a * P / (R ** 2 * T ** 2)
    B = b * P / (R * T)

    # Fugacity coefficient (phi) calculation
    print(Z)
    term1 = Z - 1

    term2 = np.log(Z - B)
    term3 = A / (2 * np.sqrt(2) * B) * np.log((Z + (1 + np.sqrt(2)) * B) / (Z + (1 - np.sqrt(2)) * B))

    ln_phi = term1 - term2 - term3
    phi = np.exp(ln_phi)
    return phi

#
# # Example Usage
# P = 10  # Pressure in bar
# T = 300  # Temperature in K
# Tc = 369.9  # Critical temperature of component in K (example: propane)
# Pc = 42.5  # Critical pressure of component in bar (example: propane)
# omega = 0.1521  # Acentric factor for propane
#
# Z_factors = peng_robinson(R, P, T, Tc, Pc, omega)
# print(f"Compressibility Factors (Z): {Z_factors}")
#
# # Example calculation for C1 at T = 120 K and P = 1 atm
# component = antoine_constants["C1"]
#
# # Take the vapor root (highest Z value for the vapor phase)
# Z_vapor = Z_factors['Z_vapor']
# Z_liquid = Z_factors['Z_liquid']
#
# # Calculate fugacity coefficient for the vapor phase
# phi_vapor = peng_robinson_fugacity_coefficient(T, P, Z_vapor, component)
# phi_liquid = peng_robinson_fugacity_coefficient(T, P, Z_liquid, component)
#
# print(f"Compressibility Factor (Z_vapor): {Z_vapor} , Compressibility Factor (Z_liquid): {Z_liquid}")
# print(f"Fugacity Coefficient (phi_vapor): {phi_vapor} , Fugacity Coefficient (phi_liquid): {phi_liquid}")
#
