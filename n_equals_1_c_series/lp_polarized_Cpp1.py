import numpy as np

def calculate_c_1_plus_plus_longitudinally_polarized(
    lepton_helicity: float,
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

        # (2): Calculate the recurrent quantity 1 + sqrt(1 + epsilon^2):
        one_plus_root_epsilon_stuff = 1. + root_one_plus_epsilon_squared

        # (3): Calculate the recurrent quantity 1 + sqrt(1 + epsilon^2) - epsilon^2:
        one_plus_root_epsilon_minus_epsilon_squared = one_plus_root_epsilon_stuff - epsilon**2

        # (4): Calculate the major term:
        major_factor = 1. - ((squared_hadronic_momentum_transfer_t / squared_Q_momentum_transfer) * (1. - 2. * x_Bjorken * (one_plus_root_epsilon_stuff + 1.) / one_plus_root_epsilon_minus_epsilon_squared))

        # (5): Calculate the prefactor:
        prefactor = -4. * lepton_helicity * target_polarization * lepton_energy_fraction_y * shorthand_k * (2. - lepton_energy_fraction_y) / root_one_plus_epsilon_squared**5

        # (6): Calculate the entire thing:
        c_1_plus_plus_LP = prefactor * one_plus_root_epsilon_minus_epsilon_squared * major_factor

        # (6.1): If verbose, log the output:
        if verbose:
            print(f"> Calculated c_1_plus_plus_LP to be:\n{c_1_plus_plus_LP}")

        # (7): Return the coefficient:
        return c_1_plus_plus_LP

    except Exception as ERROR:
        print(f"> Error in calculating c_1_plus_plus_LP for Interference Term:\n> {ERROR}")
        return 0.