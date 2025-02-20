import numpy as np

def calculate_s_2_plus_plus_longitudinally_polarized(
    target_polarization: float,
    squared_Q_momentum_transfer: float, 
    x_Bjorken: float, 
    squared_hadronic_momentum_transfer_t: float,
    epsilon: float,
    lepton_energy_fraction_y: float,
    t_prime: float,
    k_tilde: float,
    verbose: bool = False) -> float:
    """
    """

    try:

        # (1): Calculate the recurrent quantity sqrt(1 + epsilon^2):
        root_one_plus_epsilon_squared = np.sqrt(1. + epsilon**2)

        # (2): Calculate 1 + sqrt(1 + epsilon^2)
        one_plus_root_epsilon_stuff = 1. + root_one_plus_epsilon_squared

        # (3): Calculate 4 * Kt^{2} * (1 + sqrt(1 + e^{2})) * (1 + sqrt(1 + e^{2}) + xb t / Q^{2})t'/Q^{2}
        bracket_term = 4. * k_tilde**2 * (one_plus_root_epsilon_stuff - 2. * x_Bjorken) * (one_plus_root_epsilon_stuff + x_Bjorken * squared_hadronic_momentum_transfer_t / squared_Q_momentum_transfer) * t_prime / (root_one_plus_epsilon_squared * squared_Q_momentum_transfer**2)

        # (4): Calculate the prefactor
        prefactor = -4. * target_polarization * (2. - lepton_energy_fraction_y) * (1. - lepton_energy_fraction_y - (epsilon**2 * lepton_energy_fraction_y**2 / 4.)) / root_one_plus_epsilon_squared**5

        # (5): Calculate the coefficient
        s_2_plus_plus_LP = prefactor * bracket_term

        # (5.1): If verbose, log the output:
        if verbose:
            print(f"> Calculated s_2_plus_plus_LP to be:\n{s_2_plus_plus_LP}")

        # (6): Return the coefficient:
        return s_2_plus_plus_LP

    except Exception as ERROR:
        print(f"> Error in calculating s_2_plus_plus_LP for Interference Term:\n> {ERROR}")
        return 0.