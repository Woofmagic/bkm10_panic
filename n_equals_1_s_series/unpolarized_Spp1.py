import numpy as np

def calculate_s_1_plus_plus_unpolarized(
    lepton_helicity: float,
    squared_Q_momentum_transfer: float, 
    x_Bjorken: float, 
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

        # (2): Calculate the quantity t'/Q^{2}:
        tPrime_over_Q_squared = t_prime / squared_Q_momentum_transfer

        # (3): Calculate the bracket term:
        bracket_term = 1. + ((1. - x_Bjorken + 0.5 * (root_one_plus_epsilon_squared - 1.)) / root_one_plus_epsilon_squared**2) * tPrime_over_Q_squared
        
        # (4): Calculate the prefactor:
        prefactor = 8. * lepton_helicity * shorthand_k * lepton_energy_fraction_y * (2. - lepton_energy_fraction_y) / root_one_plus_epsilon_squared**2

        # (5): Calculate the coefficient
        s_1_plus_plus_unp = prefactor * bracket_term
        
        # (5.1): If verbose, log the output:
        if verbose:
            print(f"> Calculated s_1_plus_plus_unp to be:\n{s_1_plus_plus_unp}")

        # (6): Return the coefficient:
        return s_1_plus_plus_unp

    except Exception as ERROR:
        print(f"> Error in calculating s_1_plus_plus_unp for Interference Term:\n> {ERROR}")
        return 0.