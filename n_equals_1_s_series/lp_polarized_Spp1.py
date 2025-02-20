import numpy as np

def calculate_s_1_plus_plus_longitudinally_polarized(
    target_polarization: float,
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

        # (2): Calculate 1 + sqrt(1 + epsilon^2):
        one_plus_root_epsilon_stuff = 1. + root_one_plus_epsilon_squared

        # (3): Calculate the recurrent quantity t/Q^{2}
        t_over_Q_squared = squared_hadronic_momentum_transfer_t / squared_Q_momentum_transfer

        # (4): Calculate epsilon^{2} y^{2} / 4
        epsilon_y_over_2_squared = (epsilon * lepton_energy_fraction_y / 2.)**2

        # (5): Calculate the first bracket term:
        first_bracket_term = 2. * root_one_plus_epsilon_squared - 1. + (t_over_Q_squared * (one_plus_root_epsilon_stuff - 2. * x_Bjorken) / one_plus_root_epsilon_stuff)

        # (6): Calculate the second multiplicative factor:
        second_bracket_term = (3. * epsilon**2 / 2.) + (t_over_Q_squared * (1. - root_one_plus_epsilon_squared - epsilon**2 / 2. - x_Bjorken * (3.  - root_one_plus_epsilon_squared)))

        # (7): Calculate the almost prefactor:
        almost_prefactor = 4. * target_polarization * shorthand_k / root_one_plus_epsilon_squared**6

        # (8): Calculate prefactor one:
        prefactor_one = almost_prefactor * (2. - 2. * lepton_energy_fraction_y + lepton_energy_fraction_y**2 + 2. * epsilon_y_over_2_squared) * one_plus_root_epsilon_stuff

        # (9): Calculate prefactor two:
        prefactor_two = 2. * almost_prefactor * (1. - lepton_energy_fraction_y - epsilon_y_over_2_squared)
    
        # (10): Calculate the coefficient:
        s_1_plus_plus_LP = prefactor_one * first_bracket_term + prefactor_two * second_bracket_term

        # (10.1): If verbose, log the output:
        if verbose:
            print(f"> Calculated s_1_plus_plus_LP to be:\n{s_1_plus_plus_LP}")

        # (11): Return the coefficient:
        return s_1_plus_plus_LP

    except Exception as ERROR:
        print(f"> Error in calculating s_1_plus_plus_LP for Interference Term:\n> {ERROR}")
        return 0.