import numpy as np

def calculate_s_2_zero_plus_unpolarized(
    lepton_helicity: float,
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

        # (2): Calculate the recurrent quantity epsilon^2/2:
        epsilon_squared_over_2 = epsilon**2 / 2.

        # (3): Calculate the annoying y quantity:
        y_quantity = 1. - lepton_energy_fraction_y - (epsilon**2 * lepton_energy_fraction_y**2 / 4.)

        # (4): Calculate the bracket term:
        bracket_term = 1. + ((1. + epsilon_squared_over_2 / x_Bjorken) / (1. + epsilon_squared_over_2)) * x_Bjorken * squared_hadronic_momentum_transfer_t / squared_Q_momentum_transfer

        # (5): Calculate the prefactor:
        prefactor = 8. * lepton_helicity * np.sqrt(2. * y_quantity) * shorthand_k * lepton_energy_fraction_y / root_one_plus_epsilon_squared**4
        
        # (6): Calculate the coefficient:
        s_2_zero_plus_unp = prefactor * (1. + epsilon_squared_over_2) * bracket_term
        
        # (6.1): If verbose, log the output:
        if verbose:
            print(f"> Calculated s_2_zero_plus_unp to be:\n{s_2_zero_plus_unp}")

        # (7): Return the coefficient:
        return s_2_zero_plus_unp

    except Exception as ERROR:
        print(f"> Error in calculating s_2_zero_plus_unp for Interference Term:\n> {ERROR}")
        return 0.