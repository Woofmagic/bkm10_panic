import numpy as np

def calculate_c_3_plus_plus_unpolarized(
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

        # (3): Calculate the major term:
        major_term = (1. - x_Bjorken) * t_over_Q_squared + 0.5 * (root_one_plus_epsilon_squared - 1.) * (1. + t_over_Q_squared)
    
        # (4): Calculate the "intermediate" term:
        intermediate_term = (root_one_plus_epsilon_squared - 1.) / root_one_plus_epsilon_squared**5

        # (5): Calculate the prefactor:
        prefactor = -8. * shorthand_k * (1. - lepton_energy_fraction_y - epsilon**2 * lepton_energy_fraction_y**2 / 4.)
        
        # (6): Calculate the coefficient
        c_3_plus_plus_unp = prefactor * intermediate_term * major_term
        
        # (6.1): If verbose, log the output:
        if verbose:
            print(f"> Calculated c_3_plus_plus_unp to be:\n{c_3_plus_plus_unp}")

        # (7): Return the coefficient:
        return c_3_plus_plus_unp

    except Exception as ERROR:
        print(f"> Error in calculating c_3_plus_plus_unp for Interference Term:\n> {ERROR}")
        return 0.