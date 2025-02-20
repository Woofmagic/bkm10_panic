import numpy as np

def calculate_s_2_plus_plus_unpolarized(
    lepton_helicity: float,
    squared_Q_momentum_transfer: float, 
    x_Bjorken: float, 
    epsilon: float,
    lepton_energy_fraction_y: float,
    t_prime: float,
    verbose: bool = False) -> float:
    """
    """

    try:

        # (1): Calculate the recurrent quantity sqrt(1 + epsilon^2):
        root_one_plus_epsilon_squared = np.sqrt(1. + epsilon**2)

        # (2): Calculate the quantity t'/Q^{2}:
        tPrime_over_Q_squared = t_prime / squared_Q_momentum_transfer

        # (3): Calculate a fancy, annoying quantity:
        fancy_y_stuff = 1. - lepton_energy_fraction_y - epsilon**2 * lepton_energy_fraction_y**2 / 4.

        # (4): Calculate the first bracket term:
        first_bracket_term = (epsilon**2 - x_Bjorken * (root_one_plus_epsilon_squared - 1.)) / (1. + root_one_plus_epsilon_squared - 2. * x_Bjorken)

        # (5): Calculate the second bracket term:
        second_bracket_term = (2. * x_Bjorken + epsilon**2) * tPrime_over_Q_squared / (2. * root_one_plus_epsilon_squared)

        # (6): Calculate the prefactor:
        prefactor = -4. * lepton_helicity * fancy_y_stuff * lepton_energy_fraction_y * (1. + root_one_plus_epsilon_squared - 2. * x_Bjorken) * tPrime_over_Q_squared / root_one_plus_epsilon_squared**3

        # (7): Calculate the coefficient
        s_2_plus_plus_unp = prefactor * (first_bracket_term - second_bracket_term)
        
        # (7.1): If verbose, log the output:
        if verbose:
            print(f"> Calculated s_2_plus_plus_unp to be:\n{s_2_plus_plus_unp}")

        # (6): Return the coefficient:
        return s_2_plus_plus_unp

    except Exception as ERROR:
        print(f"> Error in calculating s_2_plus_plus_unp for Interference Term:\n> {ERROR}")
        return 0.