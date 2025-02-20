import numpy as np

def calculate_c_1_plus_plus_unpolarized_V(
    squared_Q_momentum_transfer: float, 
    x_Bjorken: float, 
    squared_hadronic_momentum_transfer_t: float,
    epsilon: float,
    lepton_energy_fraction_y: float,
    t_prime: float,
    shorthand_k: float,
    verbose: bool = False) -> float:
    """
    """

    try:

        # (1): Calculate the recurrent quantity sqrt(1 + epsilon^2):
        root_one_plus_epsilon_squared = np.sqrt(1. + epsilon**2)

        # (2): Calculate the recurrent quantity t/Q^{2}:
        t_over_Q_squared = squared_hadronic_momentum_transfer_t / squared_Q_momentum_transfer

        # (3): Calculate the first bracket term:
        first_bracket_term = (2. - lepton_energy_fraction_y)**2 * (1. - (1. - 2. * x_Bjorken) * t_over_Q_squared)

        # (4): Compute the first part of the second term in brackets:
        second_bracket_term_first_part = 1. - lepton_energy_fraction_y - epsilon**2 * lepton_energy_fraction_y**2 / 4.

        # (5): Compute the second part of the second term in brackets:
        second_bracket_term_second_part = 0.5 * (1. + root_one_plus_epsilon_squared - 2. * x_Bjorken) * t_prime / squared_Q_momentum_transfer

        # (6): The prefactor in front of the brackets:
        coefficient_prefactor = 16. * shorthand_k * x_Bjorken * t_over_Q_squared / np.power(root_one_plus_epsilon_squared, 5)

        # (7): The entire thing:
        c_1_plus_plus_V_unp = coefficient_prefactor * (first_bracket_term + second_bracket_term_first_part * second_bracket_term_second_part)

        # (7.1): If verbose, log the output:
        if verbose:
            print(f"> Calculated c_1_plus_plus_V_unp to be:\n{c_1_plus_plus_V_unp}")

        # (12): Return the coefficient:
        return c_1_plus_plus_V_unp

    except Exception as ERROR:
        print(f"> Error in calculating c_1_plus_plus_V_unp for Interference Term:\n> {ERROR}")
        return 0.