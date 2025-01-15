import numpy as np

def calculate_s_2_plus_plus_unpolarized_V(
    lepton_helicity: float,
    squared_Q_momentum_transfer: float, 
    x_Bjorken: float, 
    squared_hadronic_momentum_transfer_t: float,
    epsilon: float,
    lepton_energy_fraction_y: float, 
    verbose: bool = False) -> float:
    """
    """

    try:

        # (1): Calculate the recurrent quantity sqrt(1 + epsilon^2):
        root_one_plus_epsilon_squared = np.sqrt(1. + epsilon**2)

        # (2): Calculate the quantity t/Q^{2}:
        t_over_Q_squared = squared_hadronic_momentum_transfer_t / squared_Q_momentum_transfer

        # (3): Calculate a fancy, annoying quantity:
        fancy_y_stuff = 1. - lepton_energy_fraction_y - epsilon**2 * lepton_energy_fraction_y**2 / 4.

        # (4): Calculate the bracket term:
        one_minus_2xb = 1. - 2. * x_Bjorken

        # (5): Calculate the bracket term:
        bracket_term = root_one_plus_epsilon_squared - 1. + (one_minus_2xb + root_one_plus_epsilon_squared) * t_over_Q_squared

        # (6): Calculate the parentheses term:
        parentheses_term = 1. - one_minus_2xb * t_over_Q_squared

        # (7): Calculate the prefactor:
        prefactor = -4. * lepton_helicity * fancy_y_stuff * lepton_energy_fraction_y * x_Bjorken * t_over_Q_squared / root_one_plus_epsilon_squared**4

        # (8): Calculate the coefficient
        s_2_plus_plus_unp_V = prefactor * parentheses_term * bracket_term
        
        # (8.1): If verbose, log the output:
        if verbose:
            print(f"> Calculated s_2_plus_plus_unp_V to be:\n{s_2_plus_plus_unp_V}")

        # (9): Return the coefficient:
        return s_2_plus_plus_unp_V

    except Exception as ERROR:
        print(f"> Error in calculating s_2_plus_plus_unp_V for Interference Term:\n> {ERROR}")
        return 