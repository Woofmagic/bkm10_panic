_MASS_OF_PROTON_IN_GEV = .93827208816

def calculate_curly_C_longitudinally_polarized_interference(
    squared_Q_momentum_transfer: float, 
    x_Bjorken: float,
    squared_hadronic_momentum_transfer_t: float,
    Dirac_form_factor_F1: float,
    Pauli_form_factor_F2: float,
    compton_form_factor_h_real_part: float,
    compton_form_factor_h_tilde_real_part: float,
    compton_form_factor_e_real_part: float,
    compton_form_factor_e_tilde_real_part: float,
    verbose: bool = False) -> float:

    try:

        # (1): Calculate t/Q^{2}:
        t_over_Q_squared = squared_hadronic_momentum_transfer_t / squared_Q_momentum_transfer

        # (2): Calculate a fancy quantity:
        ratio_of_xb_to_more_xb = x_Bjorken / (2. - x_Bjorken + x_Bjorken * t_over_Q_squared)

        # (3): Calculate another fancy quantity that appears twice:
        x_Bjorken_correction = x_Bjorken * (1. - t_over_Q_squared) / 2.

        # (4): Calculate the first appearance of CFFs:
        first_cff_contribution = ratio_of_xb_to_more_xb * (Dirac_form_factor_F1 + Pauli_form_factor_F2) * (compton_form_factor_h_real_part + x_Bjorken_correction * compton_form_factor_e_real_part)

        # (5): Calculate the second appearance of CFFs:
        second_cff_contribution = (1. + _MASS_OF_PROTON_IN_GEV**2 * x_Bjorken * ratio_of_xb_to_more_xb * (3. + t_over_Q_squared) / squared_Q_momentum_transfer) * Dirac_form_factor_F1 * compton_form_factor_h_tilde_real_part
        
        # (6): Calculate the third appearance of CFFs:
        third_cff_contribution = t_over_Q_squared * 2 * (1. - 2. * x_Bjorken) * ratio_of_xb_to_more_xb * Pauli_form_factor_F2 * compton_form_factor_h_tilde_real_part

        # (7): Calculate the fourth appearance of the CFFs:
        fourth_cff_contribution = ratio_of_xb_to_more_xb * (x_Bjorken_correction * Dirac_form_factor_F1 + squared_hadronic_momentum_transfer_t * Pauli_form_factor_F2 / (4. * _MASS_OF_PROTON_IN_GEV**2)) * compton_form_factor_e_tilde_real_part

        # (8): Add together with the correct signs the entire thing
        curly_C_longitudinally_polarized_interference = first_cff_contribution + second_cff_contribution - third_cff_contribution - fourth_cff_contribution

        # (8.1): If verbose, log the output:
        if verbose:
            print(f"> Calculated the curly C LP for interference to be:\n{curly_C_longitudinally_polarized_interference}")
        
        # (9): Return the output:
        return curly_C_longitudinally_polarized_interference

    except Exception as ERROR:
        print(f"> Error in calculating the curly C LP contribution amplitude squared\n> {ERROR}")
        return 0