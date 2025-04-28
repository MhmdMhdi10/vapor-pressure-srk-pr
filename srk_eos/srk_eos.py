import numpy as np

# Gas constant in J·mol^−1·K^−1
R = 8.3144598

# Antoine constants for example components (C1, C3, C5, C7, H2O)
antoine_constants = {
    "C1": {"A": 3.9895, "B": 443.028, "C": -0.49, "T_c": 190.6, "P_c": 46.1, "omega": 0.011},
    "C3": {"A": 4.01158, "B": 834.26, "C": -22.763, "T_c": 369.9, "P_c": 42.5, "omega": 0.152},
    "C5": {"A": 3.9892, "B": 1070.617, "C": -40.454, "T_c": 469.8, "P_c": 33.6, "omega": 0.249},
    "C7": {"A": 4.81803, "B": 1635.409, "C": -27.338, "T_c": 540, "P_c": 27.4, "omega": 0.302},
    "H2O": {"A": 4.6543, "B": 1435.264, "C": -64.848, "T_c": 647, "P_c": 220.64, "omega": 0.344},
}

def srk_parameters(T, component):
    Tc = component["T_c"]
    Pc = component["P_c"] * 1e5  # Convert P_c from Bar to Pa
    omega = component["omega"]

    Tr = T / Tc  # Reduced temperature
    alpha = (1 + (0.48 + 1.574 * omega - 0.176 * omega ** 2) * (1 - np.sqrt(Tr))) ** 2
    a = 0.42748 * (R * Tc) ** 2 * alpha / Pc
    b = 0.08664 * R * Tc / Pc
    return a, b


def srk_eos(T, P, component):

    a, b = srk_parameters(T, component)
    A = (a * P) / (R ** 2 * T ** 2)
    B = (b * P) / (R * T)

    # SRK cubic equation coefficients
    coeffs = [1, -(1 - B), A - 2 * B - 3 * B ** 2, -(A * B - B ** 2 - B ** 3)]

    # Solve cubic equation for Z
    Z_roots = np.roots(coeffs)
    Z_real = [z.real for z in Z_roots if np.isreal(z) and z.real > 0]
    return Z_real, A, B


def srk_fugacity_coefficient(Z, A, B):
    ln_phi = Z - 1 - np.log(Z - B) - (A / (2 * np.sqrt(2) * B)) * np.log(
        (Z + (1 + np.sqrt(2)) * B) / (Z + (1 - np.sqrt(2)) * B)
    )

    if ln_phi > 700:  # np.exp(700) is the maximum safe value
        ln_phi = 700
    elif ln_phi < -700:  # Prevent underflow
        ln_phi = -700

    phi = np.exp(ln_phi)
    return phi

def srk_straight_fugacity_coefficient(T, P, component):
    a, b = srk_parameters(T, component)
    A = (a * P) / (R ** 2 * T ** 2)
    B = (b * P) / (R * T)

    # SRK cubic equation coefficients
    coeffs = [1, -(1 - B), A - 2 * B - 3 * B ** 2, -(A * B - B ** 2 - B ** 3)]

    # Solve cubic equation for Z
    Z_roots = np.roots(coeffs)
    Z_real = [z.real for z in Z_roots if np.isreal(z) and z.real > 0]


    Z = max(Z_real)

    ln_phi = Z - 1 - np.log(Z - B) - (A / (2 * np.sqrt(2) * B)) * np.log(
        (Z + (1 + np.sqrt(2)) * B) / (Z + (1 - np.sqrt(2)) * B)
    )

    if ln_phi > 700:  # np.exp(700) is the maximum safe value
        ln_phi = 700
    elif ln_phi < -700:  # Prevent underflow
        ln_phi = -700

    phi = np.exp(ln_phi)
    return phi
# # Example calculation for C1 at T = 120 K and P = 1 atm
# T = 120  # Temperature in Kelvin
# P = 101325  # Pressure in Pa (1 atm)
# component = antoine_constants["C1"]
#
# # Calculate SRK parameters, Z-factor, and fugacity coefficient
# Z_factors, A, B = srk_eos(T, P, component)
#
# # Vapor and saturated liquid Z factors
# Z_vapor = Z_factors[0]  # Choose the root closest to 1 for vapor phase
# Z_sat = Z_factors[1]  # Choose the root closest to 0 for liquid phase (saturated)
#
# # Fugacity coefficients for vapor and saturated liquid
# phi_vapor = srk_fugacity_coefficient(Z_vapor, A, B)
# phi_sat = srk_fugacity_coefficient(Z_sat, A, B)
#
# print(f"SRK EOS Parameters for C1 at {T} K and {P / 101325:.2f} atm:")
# print(f"A = {A:.4e}, B = {B:.4e}")
# print(f"Vapor Compressibility Factor (Z_v): {Z_vapor}")
# print(f"Saturated Liquid Compressibility Factor (Z_sat): {Z_sat}")
# print(f"Vapor Fugacity Coefficient (phi_v): {phi_vapor}")
# print(f"Saturated Liquid Fugacity Coefficient (phi_sat): {phi_sat}")
