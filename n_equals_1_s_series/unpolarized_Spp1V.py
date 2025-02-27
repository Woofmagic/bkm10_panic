import numpy as np

def calculate_s_1_plus_plus_unpolarized_V(
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

        # (2): Calculate the recurrent quantity t/Q^{2}:
        t_over_Q_squared = squared_hadronic_momentum_transfer_t / squared_Q_momentum_transfer

        # (3): Calculate the bracket term:
        bracket_term = root_one_plus_epsilon_squared - 1. + (1. + root_one_plus_epsilon_squared - 2. * x_Bjorken) * t_over_Q_squared

        # (4): Calculate the prefactor:
        prefactor = -8. * lepton_helicity * shorthand_k * lepton_energy_fraction_y * (2. - lepton_energy_fraction_y) * x_Bjorken * t_over_Q_squared / root_one_plus_epsilon_squared**4

        # (5): Calculate the coefficient
        s_1_plus_plus_unp_V = prefactor * bracket_term
        
        # (5.1): If verbose, log the output:
        if verbose:
            print(f"> Calculated s_1_plus_plus_unp_V to be:\n{s_1_plus_plus_unp_V}")

        # (6): Return the coefficient:
        return s_1_plus_plus_unp_V

    except Exception as ERROR:
        print(f"> Error in calculating s_1_plus_plus_unp_V for Interference Term:\n> {ERROR}")
        return 0.