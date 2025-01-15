import numpy as np

def calculate_s_1_zero_plus_unpolarized(
    lepton_helicity: float,
    squared_Q_momentum_transfer: float, 
    epsilon: float,
    lepton_energy_fraction_y: float,
    k_tilde: float,
    verbose: bool = False) -> float:
    """
    """

    try:

        # (1): Calculate the  quantity (1 + epsilon^2)^{2}:
        root_one_plus_epsilon_squared = (1. + epsilon**2)**2

        # (2): Calculate the huge y quantity:
        y_quantity = np.sqrt(1. - lepton_energy_fraction_y - (epsilon**2 * lepton_energy_fraction_y**2 / 4.))

        # (3): Calculate the coefficient
        s_1_zero_plus_unp = 8. * np.sqrt(2.) * lepton_helicity * (2. - lepton_energy_fraction_y) * lepton_energy_fraction_y * y_quantity * k_tilde**2 / (root_one_plus_epsilon_squared * squared_Q_momentum_transfer)
        
        # (3.1): If verbose, log the output:
        if verbose:
            print(f"> Calculated s_1_zero_plus_unp to be:\n{s_1_zero_plus_unp}")

        # (4): Return the coefficient:
        return s_1_zero_plus_unp

    except Exception as ERROR:
        print(f"> Error in calculating s_1_zero_plus_unp for Interference Term:\n> {ERROR}")
        return 0.