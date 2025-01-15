import numpy as np

def calculate_c_3_plus_plus_unpolarized_A(
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

        # (1): Calculate the main term:
        main_term = squared_hadronic_momentum_transfer_t * t_prime * (x_Bjorken * (1. - x_Bjorken) + epsilon**2 / 4.) / squared_Q_momentum_transfer**2

        # (2): Calculate the prefactor: 
        prefactor = 16. * shorthand_k * (1. - lepton_energy_fraction_y - epsilon**2 * lepton_energy_fraction_y**2 / 4.) / (1. + epsilon**2)**2.5
        
        # (3): The entire thing:
        c_3_plus_plus_A_unp = prefactor * main_term

        # (3.1): If verbose, log the output:
        if verbose:
            print(f"> Calculated c_3_plus_plus_A_unp to be:\n{c_3_plus_plus_A_unp}")

        # (4): Return the coefficient:
        return c_3_plus_plus_A_unp

    except Exception as ERROR:
        print(f"> Error in calculating c_2_plus_plus_A_unp for Interference Term:\n> {ERROR}")
        return 0.