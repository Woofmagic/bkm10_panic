import numpy as np

def calculate_c_2_plus_plus_unpolarized_A(
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

        # (2): Calculate the recurrent quantity t/Q^{2}:
        t_over_Q_squared = squared_hadronic_momentum_transfer_t / squared_Q_momentum_transfer

        # (3): Calculate t'/Q^{2}
        t_prime_over_Q_squared = t_prime / squared_Q_momentum_transfer

        # (4): Calculate the first bracket term:
        first_bracket_term = 4. * (1. - 2. * x_Bjorken) * k_tilde**2 / (root_one_plus_epsilon_squared * squared_Q_momentum_transfer)

        # (5): Calculate the second bracket term:
        second_bracket_term = (3.  - root_one_plus_epsilon_squared - 2. * x_Bjorken + epsilon**2 / x_Bjorken ) * x_Bjorken * t_prime_over_Q_squared

        # (6): Calculate the prefactor: 
        prefactor = 4. * (2. - lepton_energy_fraction_y) * (1. - lepton_energy_fraction_y - epsilon**2 * lepton_energy_fraction_y**2 / 4.) * t_over_Q_squared / root_one_plus_epsilon_squared**4
        
        # (7): The entire thing:
        c_2_plus_plus_A_unp = prefactor * (first_bracket_term - second_bracket_term)

        # (7.1): If verbose, log the output:
        if verbose:
            print(f"> Calculated c_2_plus_plus_A_unp to be:\n{c_2_plus_plus_A_unp}")

        # (8): Return the coefficient:
        return c_2_plus_plus_A_unp

    except Exception as ERROR:
        print(f"> Error in calculating c_2_plus_plus_A_unp for Interference Term:\n> {ERROR}")
        return 0.