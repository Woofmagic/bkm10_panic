_MASS_OF_PROTON_IN_GEV = .93827208816

def calculate_curly_C_unpolarized_interference(
    squared_Q_momentum_transfer: float,
    x_Bjorken: float,
    squared_hadronic_momentum_transfer_t: float,
    Dirac_form_factor_F1: float,
    Pauli_form_factor_F2: float,
    compton_form_factor_h: float,
    compton_form_factor_h_tilde: float,
    compton_form_factor_e: float,
    verbose: bool = False) -> float:

    try:

        # (1): Calculate the first two terms: weighted CFFs:
        weighted_cffs = (Dirac_form_factor_F1 * compton_form_factor_h) - (squared_hadronic_momentum_transfer_t * Pauli_form_factor_F2 * compton_form_factor_e / (4. * _MASS_OF_PROTON_IN_GEV**2))

        # (2): Calculate the next term:
        second_term = x_Bjorken * (Dirac_form_factor_F1 + Pauli_form_factor_F2) * compton_form_factor_h_tilde / (2. - x_Bjorken + (x_Bjorken * squared_hadronic_momentum_transfer_t / squared_Q_momentum_transfer))

        # (3): Add them together:
        curly_C_unpolarized_interference = weighted_cffs + second_term

        # (4.1): If verbose, print the calculation:
        if verbose:
            print(f"> Calculated Curly C interference unpolarized target to be:\n{curly_C_unpolarized_interference}")

        # (5): Return the output:
        return curly_C_unpolarized_interference

    except Exception as ERROR:
        print(f"> Error in calculating the Curly C interference unpolarized target: \n> {ERROR}")
        return 0.