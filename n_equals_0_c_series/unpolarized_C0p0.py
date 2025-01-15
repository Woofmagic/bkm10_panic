import numpy as np

def calculate_c_0_zero_plus_unpolarized(
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

        # (1): Calculate the bracket quantity:
        bracket_quantity = epsilon**2 + squared_hadronic_momentum_transfer_t * (2. - 6.* x_Bjorken - epsilon**2) / (3. * squared_Q_momentum_transfer)
        
        # (2): Calculate part of the prefactor:
        prefactor = 12. * np.sqrt(2.) * shorthand_k * (2. - lepton_energy_fraction_y) * np.sqrt(1. - lepton_energy_fraction_y - (epsilon**2 * lepton_energy_fraction_y**2 / 4)) / np.power(1. + epsilon**2, 2.5)
        
        # (3): Calculate the coefficient:
        c_0_zero_plus_unp = prefactor * bracket_quantity
        
        # (3.1): If verbose, log the output:
        if verbose:
            print(f"> Calculated c_0_zero_plus_unp to be:\n{c_0_zero_plus_unp}")

        # (4): Return the coefficient:
        return c_0_zero_plus_unp

    except Exception as ERROR:
        print(f"> Error in calculating c_0_zero_plus_unp for Interference Term:\n> {ERROR}")
        return 0.