import numpy as np

def calculate_s_3_plus_plus_longitudinally_polarized(
    target_polarization: float,
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

        # (2): Calculate 1 + sqrt(1 + epsilon^2):
        one_plus_root_epsilon_stuff = 1. + root_one_plus_epsilon_squared

        # (3): Calculate the coefficient
        prefactor = -4. * target_polarization * shorthand_k * (1. - lepton_energy_fraction_y - lepton_energy_fraction_y**2 * epsilon**2 / 4.) / root_one_plus_epsilon_squared**6

        # (4): Calculate the coefficient:
        s_3_plus_plus_LP = prefactor * (one_plus_root_epsilon_stuff - 2. * x_Bjorken) * epsilon**2 * t_prime / (squared_Q_momentum_transfer * one_plus_root_epsilon_stuff)

        # (4.1): If verbose, log the output:
        if verbose:
            print(f"> Calculated s_3_plus_plus_LP to be:\n{s_3_plus_plus_LP}")

        # (5): Return the coefficient:
        return s_3_plus_plus_LP

    except Exception as ERROR:
        print(f"> Error in calculating s_3_plus_plus_LP for Interference Term:\n> {ERROR}")
        return 0