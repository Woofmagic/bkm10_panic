from decimal import Decimal

import numpy as np

def calculate_c_0_zero_plus_longitudinally_polarized(
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
    # Title: `calculate_c_0_zero_plus_longitudinally_polarized`

    ## Description: 
    We calculate the coefficient C++(n = 0) for the longitudinally-polarized
    target.

    ## Arguments:
    
    1. lepton_helicity (float)

    The helicity of the lepton beam. The number, while a float, 
    is usually either -1.0 or +1.0.

    ## Returns:
    
    1. c_0_zero_plus_LP (float)

    ## Examples:
    None
    """

    try:

        # (1): Calculate the annoying quantity sqrt(1 - y - y^{2} epsilon^{2} / 2)
        root_combination_of_y_and_epsilon = np.sqrt(1. - lepton_energy_fraction_y - (lepton_energy_fraction_y**2 * epsilon**2 / 4.))

        # (2): Calculate the "prefactor":
        prefactor = 8. * np.sqrt(2.) * lepton_helicity * target_polarization * shorthand_k * (1. - x_Bjorken) * lepton_energy_fraction_y / (1. + epsilon**2)**2

        # (3): Calculate everything:
        c_0_zero_plus_LP = prefactor * root_combination_of_y_and_epsilon * squared_hadronic_momentum_transfer_t / squared_Q_momentum_transfer

        # (3.1): If verbose, log the output:
        if verbose:
            print(f"> Calculated c_0_zero_plus_LP to be:\n{c_0_zero_plus_LP}")

        # (4): Return the coefficient:
        return c_0_zero_plus_LP

    except Exception as ERROR:
        print(f"> Error in calculating c_0_zero_plus_LP for Interference Term:\n> {ERROR}")
        return 0.