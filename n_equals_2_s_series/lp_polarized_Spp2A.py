from decimal import Decimal

import numpy as np

def calculate_s_2_plus_plus_longitudinally_polarized_A(
    target_polarization: float,
    squared_Q_momentum_transfer: float, 
    x_Bjorken: float, 
    squared_hadronic_momentum_transfer_t: float,
    epsilon: float,
    lepton_energy_fraction_y: float,
    t_prime: float,
    k_tilde: float,
    verbose: bool = False) -> float:
    """
    """

    try:

        # (1): Calculate the recurrent quantity sqrt(1 + epsilon^2):
        root_one_plus_epsilon_squared = np.sqrt(1. + epsilon**2)

        # (2): Calculate the first contribution to the bracket term:
        bracket_term_first_term = (1. + root_one_plus_epsilon_squared - 2. * x_Bjorken) * (1. - ((1. - 2. * x_Bjorken) * squared_hadronic_momentum_transfer_t / squared_Q_momentum_transfer)) * t_prime / squared_Q_momentum_transfer

        # (3): Calculate second contribution to the bracket term:
        bracket_term_second_term = 4. * k_tilde**2 / squared_Q_momentum_transfer

        # (4): Calculate the bracket term:
        bracket_term = x_Bjorken * squared_hadronic_momentum_transfer_t * (bracket_term_second_term - bracket_term_first_term) / squared_Q_momentum_transfer

        # (5): Calculate the prefactor:
        prefactor = 4. * target_polarization * (2. - lepton_energy_fraction_y) * (1. - lepton_energy_fraction_y - epsilon**2 * lepton_energy_fraction_y**2 / 4.) / root_one_plus_epsilon_squared**5

        # (6): Calculate the coefficient
        s_2_plus_plus_A_LP = prefactor * bracket_term

        # (6.1): If verbose, log the output:
        if verbose:
            print(f"> Calculated s_2_plus_plus_A_LP to be:\n{s_2_plus_plus_A_LP}")

        # (7): Return the coefficient:
        return s_2_plus_plus_A_LP

    except Exception as ERROR:
        print(f"> Error in calculating s_2_plus_plus_A_LP for Interference Term:\n> {ERROR}")
        return 0