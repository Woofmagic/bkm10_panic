import numpy as np

def calculate_c_0_plus_plus_longitudinally_polarized(
    lepton_helicity: float,
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

        # (1): Calculate the recurrent quantity sqrt(1 + epsilon^2):
        root_one_plus_epsilon_squared = np.sqrt(1. + epsilon**2)

        # (2): Calculate the recurrent quantity t/Q^{2}:
        t_over_Q_squared = squared_hadronic_momentum_transfer_t / squared_Q_momentum_transfer

        # (3): Calculate the first term in the brackets: 
        first_bracket_term = (2. - lepton_energy_fraction_y)**2 * k_tilde**2 / squared_Q_momentum_transfer

        # (4): Calculate the first part of the second term in brackets:
        second_bracket_term_first_part = 1. - lepton_energy_fraction_y + (epsilon**2 * lepton_energy_fraction_y**2 / 4.)

        # (5): Calculate the second part of the second term in brackets:
        second_bracket_term_second_part = x_Bjorken * t_over_Q_squared - (epsilon**2 * (1. - t_over_Q_squared) / 2.)

        # (6): Calculate the third part of the second term in brackets:
        second_bracket_term_third_part = 1. + t_over_Q_squared * ((root_one_plus_epsilon_squared - 1. + 2. * x_Bjorken) / (1. + root_one_plus_epsilon_squared))

        # (7): Stitch together the second bracket term:
        second_bracket_term = second_bracket_term_first_part * second_bracket_term_second_part * second_bracket_term_third_part

        # (8): Calculate the prefactor:
        prefactor = -4. * lepton_helicity * target_polarization * lepton_energy_fraction_y * (1. + root_one_plus_epsilon_squared) / root_one_plus_epsilon_squared**5

        # (9): Calculate the entire thing:
        c_0_plus_plus_LP = prefactor * (first_bracket_term + second_bracket_term)

        # (9.1): If verbose, log the output:
        if verbose:
            print(f"> Calculated c_0_plus_plus_LP to be:\n{c_0_plus_plus_LP}")

        # (10): Return the coefficient:
        return c_0_plus_plus_LP

    except Exception as ERROR:
        print(f"> Error in calculating c_0_plus_plus_LP for Interference Term:\n> {ERROR}")
        return 0.