import numpy as np

def fl_calculator(Z_sat, P, P_sat, phi_sat, R, T):

    v_sat = Z_sat * R * T / P_sat

    fl = abs(phi_sat * P_sat * np.exp((v_sat * (P - P_sat) / R * T)))
    return fl

