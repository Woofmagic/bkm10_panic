import numpy as np

def calculate_c_1_plus_plus_unpolarized(
    squared_Q_momentum_transfer: float, 
    x_Bjorken: float, 
    squared_hadronic_momentum_transfer_t: float,
    epsilon: float,
    lepton_energy_fraction_y: float, 
    shorthand_k: float,
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

        # (4): Calculate first term in first brackets
        first_bracket_first_term = (1. + (1. - x_Bjorken) * (root_one_plus_epsilon_squared - 1.) / (2. * x_Bjorken) + epsilon**2 / (4. * x_Bjorken)) * x_Bjorken * t_over_Q_squared

        # (5): Calculate the first bracket term:
        first_bracket_term = first_bracket_first_term - 3. * epsilon**2 / 4.

        # (6): Calculate the second bracket term:
        second_bracket_term = 1. - (1. - 3. * x_Bjorken) * t_over_Q_squared + (1. - root_one_plus_epsilon_squared + 3. * epsilon**2) * x_Bjorken * t_over_Q_squared / (one_plus_root_epsilon_stuff - epsilon**2)

        # (7): Calculate the crazy coefficient with all the y's:
        fancy_y_coefficient = 2. - 2. * lepton_energy_fraction_y + lepton_energy_fraction_y**2 + epsilon**2 * lepton_energy_fraction_y**2 / 2.

        # (8): Calculate the entire second term:
        second_term = 4. * shorthand_k * fancy_y_coefficient * (one_plus_root_epsilon_stuff - epsilon**2) * second_bracket_term / root_one_plus_epsilon_squared**5

        # (9): Calculate the first term:
        first_term = 16. * shorthand_k * (1. - lepton_energy_fraction_y - epsilon**2 * lepton_energy_fraction_y**2 / 4.) * first_bracket_term / root_one_plus_epsilon_squared**5

        # (10): Calculate the coefficient
        c_1_plus_plus_unp = first_term + second_term
        
        # (11.1): If verbose, log the output:
        if verbose:
            print(f"> Calculated c_1_plus_plus_unp to be:\n{c_1_plus_plus_unp}")

        # (12): Return the coefficient:
        return c_1_plus_plus_unp

    except Exception as ERROR:
        print(f"> Error in calculating c_1_plus_plus_unp for Interference Term:\n> {ERROR}")
        return 0.