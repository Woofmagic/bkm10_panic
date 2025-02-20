from decimal import Decimal

import numpy as np

def calculate_s_1_zero_plus_longitudinally_polarized_V(
    target_polarization: float,
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

        # (1): Calculate the annoying quantity 1 - y - y^{2} epsilon^{2} / 4
        combination_of_y_and_epsilon = 1. - lepton_energy_fraction_y - (lepton_energy_fraction_y**2 * epsilon**2 / 4.)

        # (2): Calculate t/Q^{2}:
        t_over_Q_squared = squared_hadronic_momentum_transfer_t / squared_Q_momentum_transfer

        # (3): Calculate first bracket term:
        first_bracket_term = k_tilde**2 * (2. - lepton_energy_fraction_y)**2 / squared_Q_momentum_transfer

        # (4): Calculate a long contribution to the second bracket term:
        second_bracket_term_long = 4. - 2. * x_Bjorken + 3. * epsilon**2 + t_over_Q_squared * (4. * x_Bjorken * (1. - x_Bjorken) + epsilon**2)

        # (5): Calculate the second bracket term:
        second_bracket_term = (1. + t_over_Q_squared) * combination_of_y_and_epsilon * second_bracket_term_long
        
        # (6): Calculate the prefactor:
        prefactor = -8. * np.sqrt(2.) * target_polarization  * np.sqrt(combination_of_y_and_epsilon) * t_over_Q_squared / np.sqrt((1. + epsilon**2)**5)

        # (7): Calculate everything:
        s_1_zero_plus_V_LP = prefactor * (first_bracket_term + second_bracket_term)

        # (6.1): If verbose, log the output:
        if verbose:
            print(f"> Calculated s_1_zero_plus_V_LP to be:\n{s_1_zero_plus_V_LP}")

        # (7): Return the coefficient:
        return s_1_zero_plus_V_LP

    except Exception as ERROR:
        print(f"> Error in calculating s_1_zero_plus_V_LP for Interference Term:\n> {ERROR}")
        return 0.