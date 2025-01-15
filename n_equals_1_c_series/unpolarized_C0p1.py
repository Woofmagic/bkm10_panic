import numpy as np

def calculate_c_1_zero_plus_unpolarized(
    squared_Q_momentum_transfer: float, 
    x_Bjorken: float, 
    squared_hadronic_momentum_transfer_t: float,
    epsilon: float,
    lepton_energy_fraction_y: float, 
    t_prime: float,
    verbose: bool = False) -> float:
    """
    """

    try:

        # (1): Calculate the recurrent quantity sqrt(1 + epsilon^2):
        root_one_plus_epsilon_squared = np.sqrt(1. + epsilon**2)

        # (2): Calculate the recurrent quantity t/Q^{2}:
        t_over_Q_squared = squared_hadronic_momentum_transfer_t / squared_Q_momentum_transfer

        # (3): Calculate t'/Q^{2}
        t_prime_over_Q_squared = t_prime / squared_Q_momentum_transfer

        # (4): Calculate 1 - x_{B}:
        one_minus_xb = 1. - x_Bjorken

        # (5): Calculate the annoying y quantity:
        y_quantity = 1. - lepton_energy_fraction_y - (epsilon**2 * lepton_energy_fraction_y**2 / 4.)

        # (6): Calculate the first term:
        first_bracket_term = (2. - lepton_energy_fraction_y)**2 * t_prime_over_Q_squared * (one_minus_xb + (one_minus_xb * x_Bjorken + (epsilon**2 / 4.)) * t_prime_over_Q_squared / root_one_plus_epsilon_squared)
        
        # (7): Calculate the second term:
        second_bracket_term = y_quantity * (1. - (1. - 2. * x_Bjorken) * t_over_Q_squared) * (epsilon**2 - 2. * (1. + (epsilon**2 / (2. * x_Bjorken))) * x_Bjorken * t_over_Q_squared)
        
        # (8): Calculate part of the prefactor:
        prefactor = 8. * np.sqrt(2. * y_quantity) / root_one_plus_epsilon_squared**4
        
        # (9): Calculate the coefficient:
        c_1_zero_plus_unp = prefactor * (first_bracket_term + second_bracket_term)
        
        # (9.1): If verbose, log the output:
        if verbose:
            print(f"> Calculated c_1_zero_plus_unp to be:\n{c_1_zero_plus_unp}")

        # (9): Return the coefficient:
        return c_1_zero_plus_unp

    except Exception as ERROR:
        print(f"> Error in calculating c_1_zero_plus_unp for Interference Term:\n> {ERROR}")
        return 0.