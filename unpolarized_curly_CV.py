def calculate_curly_C_unpolarized_interference_V(
    squared_Q_momentum_transfer: float, 
    x_Bjorken: float,
    squared_hadronic_momentum_transfer_t: float,
    Dirac_form_factor_F1: float,
    Pauli_form_factor_F2: float,
    compton_form_factor_h_real_part: float,
    compton_form_factor_e_real_part: float,
    verbose: bool = False) -> float:

    try:

        # (1): Calculate the first two terms: weighted CFFs:
        cff_term = compton_form_factor_h_real_part + compton_form_factor_e_real_part

        # (2): Calculate the next term:
        second_term = x_Bjorken * (Dirac_form_factor_F1 + Pauli_form_factor_F2) / (2. - x_Bjorken + (x_Bjorken * squared_hadronic_momentum_transfer_t / squared_Q_momentum_transfer))

        # (3): Add them together:
        curly_C_unpolarized_interference_V = cff_term * second_term

        # (4.1): If verbose, print the calculation:
        if verbose:
            print(f"> Calculated Curly C interference V unpolarized target to be:\n{curly_C_unpolarized_interference_V}")

        # (5): Return the output:
        return curly_C_unpolarized_interference_V

    except Exception as ERROR:
        print(f"> Error in calculating the Curly C interference V unpolarized target: \n> {ERROR}")
        return 0.