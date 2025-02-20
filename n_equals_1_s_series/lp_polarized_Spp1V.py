import numpy as np

def calculate_s_1_plus_plus_longitudinally_polarized_V(
    target_polarization: float,
    squared_Q_momentum_transfer: float, 
    x_Bjorken: float, 
    squared_hadronic_momentum_transfer_t: float,
    epsilon: float,
    lepton_energy_fraction_y: float, 
    t_prime: float,
    shorthand_k: float,
    verbose: bool = False) -> float:
    """
    """

    try:

        # (1): Calculate epsilon squared:
        ep_squared = epsilon**2

        # (2): Calculate the recurrent quantity sqrt(1 + epsilon^2):
        root_one_plus_epsilon_squared = np.sqrt(1. + ep_squared)

        # (3): Calculate the recurrent quantity t/Q^{2}
        t_over_Q_squared = squared_hadronic_momentum_transfer_t / squared_Q_momentum_transfer

        # (4): Calculate the quantity t'/Q^{2}
        t_prime_over_Q_squared = t_prime / squared_Q_momentum_transfer

        # (5): Calculate epsilon^{2} y^{2} / 4
        epsilon_y_over_2_squared = ep_squared * lepton_energy_fraction_y**2 / 4.

        # (6): Calculate the first bracket term:
        first_bracket_term = 1. - (t_prime_over_Q_squared * ((1. - 2. * x_Bjorken) * (1. - 2. * x_Bjorken + root_one_plus_epsilon_squared)) / (2. * root_one_plus_epsilon_squared**2))

        # (7): Calculate the second multiplicative factor:
        second_term_parentheses_term = t_over_Q_squared * (1. - (x_Bjorken * ((3. + root_one_plus_epsilon_squared) / 4.)) + (5. * ep_squared / 8.))

        # (8): Calculate the numerator of the second term in brackets
        second_bracket_term_numerator = 1. - root_one_plus_epsilon_squared + (ep_squared / 2.) - (2. * x_Bjorken * (3. * (1. - x_Bjorken) - root_one_plus_epsilon_squared))
        
        # (9): Calculate the denominator of the second term in brackets
        second_bracket_term_denominator = 4. - (x_Bjorken * (root_one_plus_epsilon_squared + 3.)) + (5. * ep_squared / 2.)

        # (10): Calculate the second bracket term:
        second_bracket_term = 1. - (t_over_Q_squared * second_bracket_term_numerator / second_bracket_term_denominator)
        
        # (11): Calculate the almost_prefactor:
        almost_prefactor = 8. * target_polarization * shorthand_k / root_one_plus_epsilon_squared**4

        # (12): Calculate the first prefactor:
        prefactor_one = almost_prefactor * (2. - 2. * lepton_energy_fraction_y + lepton_energy_fraction_y**2 + 2. * epsilon_y_over_2_squared) * t_over_Q_squared

        # (13): Calculate the second prefactor:
        prefactor_two = 4. * almost_prefactor * (1. - lepton_energy_fraction_y - epsilon_y_over_2_squared) / root_one_plus_epsilon_squared**2

        # (14): Calculate the coefficient:
        s_1_plus_plus_V_LP = prefactor_one * first_bracket_term + prefactor_two * second_term_parentheses_term * second_bracket_term

        # (14.1): If verbose, log the output:
        if verbose:
            print(f"> Calculated s_1_plus_plus_V_LP to be:\n{s_1_plus_plus_V_LP}")

        # (15): Return the coefficient:
        return s_1_plus_plus_V_LP

    except Exception as ERROR:
        print(f"> Error in calculating s_1_plus_plus_V_LP for Interference Term:\n> {ERROR}")
        return 0.