import numpy as np
from Utility.fl_calculator import fl_calculator
from Utility.fv_calculator import fv_calculator



def fugacity_error(fv, fl):
    """
    Calculate the absolute error between the vapor and liquid fugacity coefficients.
    """
    error = abs(fv - fl)
    return error


def fugacity_derivative(Z_func, phi_func, R, P, T, component, delta=0.01):
    """
    Numerically calculate the derivative of a function with respect to pressure (P) or temperature (T).

    Parameters:
        func (callable): The function whose derivative we want to calculate (e.g., fugacity coefficient).
        P (float): Current value of pressure.
        T (float): Current value of temperature.
        delta (float): A small perturbation for the finite difference.

    Returns:
        float: The derivative of the function with respect to P or T.
    """

    def fp_calculator(Z_func, phi_func, R, P, T, component):
        P = P + delta
        Z_factors, A, B = Z_func(T, P, component)

        Z_vapor = max(Z_factors)  # Choose the root closest to 1 for vapor phase
        Z_sat = min(Z_factors)  # Choose the root closest to 0 for liquid phase (saturated)
        print(Z_factors)


        phi_vapor = phi_func(Z_vapor, A, B)
        phi_sat = phi_func(Z_sat, A, B)


        print("phi : ", phi_vapor, phi_sat)
        # print(phi_vapor)

        fv = fv_calculator(phi_vapor, P)
        fl = fl_calculator(Z_sat, P, P, phi_sat, R, T)

        print("ff + :  ", fv, fl)

        fp = abs(fv - fl)
        return fp

    def fm_calculator(Z_func, phi_func, R, P, T, component):
        P = P - delta

        Z_factors, A, B = Z_func(T, P, component)

        Z_vapor = max(Z_factors)  # Choose the root closest to 1 for vapor phase
        Z_sat = min(Z_factors)  # Choose the root closest to 0 for liquid phase (saturated)
        print("z - ", Z_factors)

        phi_vapor = phi_func(Z_vapor, A, B)
        phi_sat = phi_func(Z_sat, A, B)

        fv = fv_calculator(phi_vapor, P)
        fl = fl_calculator(Z_sat, P, P, phi_sat, R, T)

        print("ff - :  ", fv, fl)
        fm = abs(fv - fl)
        return fm

    f_plus_delta = fp_calculator(Z_func, phi_func, R, P, T, component)  # Function value at P + delta

    f_minus_delta = fm_calculator(Z_func, phi_func, R, P, T, component)  # Function value at P - delta

    print("f_func", f_plus_delta, f_minus_delta)


    derivative = f_plus_delta - f_minus_delta / (2 * delta)  # Central difference

    return derivative


def newton_raphson(Z_func, phi_func, R, P, T, component, fv, fl, tolerance=1e-5, max_iter=100):
    """
    Solve for the root of the function using the Newton-Raphson method.

    Parameters:
        dfunc (callable): The function whose root we are trying to find (e.g., fugacity error).
        P (float): Current value of pressure.
        T (float): Current value of temperature.
        tolerance (float): Desired tolerance for the solution.
        max_iter (int): Maximum number of iterations.

    Returns:
        float: The solution (root) of the function.
    """

    for _ in range(max_iter):
        f_value = fugacity_error(fv, fl)  # Evaluate the function
        df_value = fugacity_derivative(Z_func, phi_func, R, P, T, component)  # Evaluate the derivative of the function

        if f_value == 0:
            print("Derivative is zero. Cannot proceed with Newton-Raphson.")
            break

        # Newton-Raphson update
        P_new = P - f_value / df_value

        # Check for convergence
        if abs(P_new - P) < tolerance:
            print(f"Converged after {_ + 1} iterations.")
            return P_new

        P = P_new

    print("Max iterations reached without convergence.")
    return P  # Return the last computed value if max_iter is reached
