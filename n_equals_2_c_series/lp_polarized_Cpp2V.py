from decimal import Decimal

import numpy as np

def calculate_c_2_plus_plus_longitudinally_polarized_V(
    lepton_helicity: float,
    target_polarization: float,
    squared_Q_momentum_transfer: float, 
    x_Bjorken: float, 
    squared_hadronic_momentum_transfer_t: float,
    epsilon: float,
    lepton_energy_fraction_y: float, 
    verbose: bool = False) -> float:
    """
    """

    try:

        # (1): Calculate the recurrent quantity sqrt(1 + epsilon^2):
        root_one_plus_epsilon_squared = np.sqrt(1. + epsilon**2)

        # (2): Calculate the recurrent quantity t/Q^{2}
        t_over_Q_squared = squared_hadronic_momentum_transfer_t / squared_Q_momentum_transfer

        # (3): Calculate the recurrent quantity 1 + sqrt(1 + epsilon^2):
        one_plus_root_epsilon_stuff = 1. + root_one_plus_epsilon_squared

        # (4): Calculate one of the multiplicative factors:
        first_multiplicative_factor = (one_plus_root_epsilon_stuff - 2.) + t_over_Q_squared * (one_plus_root_epsilon_stuff - 2. * x_Bjorken)

        # (5): Calculate the second multiplicative factor:
        second_multiplicative_factor = 1. + (t_over_Q_squared * (4. * (1. - x_Bjorken) * x_Bjorken + epsilon**2 ) / (4. - 2. * x_Bjorken + 3. * epsilon**2))

        # (6): Calculate the second multiplicative factor:
        third_multiplicative_factor = t_over_Q_squared * (4. - 2. * x_Bjorken + 3. * epsilon**2)

        # (7): Calculate the prefactor:
        prefactor = -2. * lepton_helicity * target_polarization * lepton_energy_fraction_y * (1. - lepton_energy_fraction_y - (lepton_energy_fraction_y**2 * epsilon**2 / 4.)) / root_one_plus_epsilon_squared**5

        # (8): Calculate the entire thing:
        c_2_plus_plus_V_LP = prefactor * first_multiplicative_factor * second_multiplicative_factor * third_multiplicative_factor

        # (8.1): If verbose, log the output:
        if verbose:
            print(f"> Calculated c_2_plus_plus_V_LP to be:\n{c_2_plus_plus_V_LP}")

        # (9): Return the coefficient:
        return c_2_plus_plus_V_LP

    except Exception as ERROR:
        print(f"> Error in calculating c_2_plus_plus_V_LP for Interference Term:\n> {ERROR}")
        return 0.