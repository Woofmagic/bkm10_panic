import numpy as np

def calculate_c_0_plus_plus_unpolarized(
    squared_Q_momentum_transfer: float,
    x_Bjorken: float,
    squared_hadronic_momentum_transfer_t: float,
    epsilon: float,
    lepton_energy_fraction_y: float,
    k_tilde: float,
    verbose: bool = False) -> float:
    """
    """

    try:

        # (1): Calculate the recurrent quantity sqrt(1 + epsilon^2):
        root_one_plus_epsilon_squared = np.sqrt(1. + epsilon**2)

        # (2): Calculate the recurrent quantity t/Q^{2}:
        t_over_Q_squared = squared_hadronic_momentum_transfer_t / squared_Q_momentum_transfer

        # (3): Calculate 1 + sqrt(1 + epsilon^{2}):
        one_plus_root_epsilon_stuff = 1. + root_one_plus_epsilon_squared

        # (4): Calculate 2 - x_{B}:
        two_minus_xb = 2. - x_Bjorken

        # (5): Caluclate 2 - y:
        two_minus_y = 2. - lepton_energy_fraction_y

        # (6): Calculate the first term in the brackets:
        first_term_in_brackets = k_tilde**2 * two_minus_y**2 / (squared_Q_momentum_transfer * root_one_plus_epsilon_squared)

        # (7): Calculate the first part of the second term in brackets:
        second_term_in_brackets_first_part = t_over_Q_squared * two_minus_xb * (1. - lepton_energy_fraction_y - (epsilon**2 * lepton_energy_fraction_y**2 / 4.))
        
        # (8): Calculate the numerator of the second part of the second term in brackets:
        second_term_in_brackets_second_part_numerator = 2. * x_Bjorken * t_over_Q_squared * (two_minus_xb + 0.5 * (root_one_plus_epsilon_squared - 1.) + 0.5 * epsilon**2 / x_Bjorken) + epsilon**2
        
        # (9): Calculate the second part of the second term in brackets:
        second_term_in_brackets_second_part =  1. + second_term_in_brackets_second_part_numerator / (two_minus_xb * one_plus_root_epsilon_stuff)
        
        # (10): Calculate the prefactor:
        prefactor = -4. * two_minus_y * one_plus_root_epsilon_stuff / np.power(root_one_plus_epsilon_squared, 4)

        # (11): Calculate the coefficient
        c_0_plus_plus_unp = prefactor * (first_term_in_brackets + second_term_in_brackets_first_part * second_term_in_brackets_second_part)

        # (11.1): If verbose, log the output:
        if verbose:
            print(f"> Calculated c_0_plus_plus_unp to be:\n{c_0_plus_plus_unp}")

        # (12): Return the coefficient:
        return c_0_plus_plus_unp

    except Exception as ERROR:
        print(f"> Error in calculating c_0_plus_plus_unp for Interference Term:\n> {ERROR}")
        return 0.