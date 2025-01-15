import numpy as np

def calculate_s_1_zero_plus_unpolarized_V(
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

        # (1): Calculate the quantity (1 + epsilon^2)^{2}:
        one_plus_epsilon_squared_squared = (1. + epsilon**2)**2

        # (2): Calculate the quantity t/Q^{2}:
        t_over_Q_squared = squared_hadronic_momentum_transfer_t / squared_Q_momentum_transfer

        # (3): Calculate a fancy, annoying quantity:
        fancy_y_stuff = 1. - lepton_energy_fraction_y - epsilon**2 * lepton_energy_fraction_y**2 / 4.

        # (4): Calculate the bracket term:
        bracket_term = 4. * (1. - 2. * x_Bjorken) * t_over_Q_squared * (1. + x_Bjorken * t_over_Q_squared) + epsilon**2 * (1. + t_over_Q_squared)**2

        # (5): Calculate the prefactor:
        prefactor = 4. * np.sqrt(2. * fancy_y_stuff) * lepton_helicity * lepton_energy_fraction_y * (2. - lepton_energy_fraction_y) * x_Bjorken * t_over_Q_squared / one_plus_epsilon_squared_squared

        # (6): Calculate the coefficient
        s_1_zero_plus_unp_V = prefactor * bracket_term
        
        # (6.1): If verbose, log the output:
        if verbose:
            print(f"> Calculated s_1_zero_plus_unp_V to be: {s_1_zero_plus_unp_V}")

        # (7): Return the coefficient:
        return s_1_zero_plus_unp_V

    except Exception as ERROR:
        print(f"> Error in calculating s_1_zero_plus_unp_V for Interference Term:\n> {ERROR}")
        return 0.