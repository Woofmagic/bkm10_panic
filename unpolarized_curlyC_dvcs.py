_MASS_OF_PROTON_IN_GEV = .93827208816

def calculate_curly_c_unpolarized_dvcs(
    squared_Q_momentum_transfer: float,
    x_Bjorken: float,
    squared_hadronic_momentum_transfer_t: float,
    epsilon: float,
    compton_form_factor_h: complex,
    compton_form_factor_h_tilde: complex,
    compton_form_factor_e: complex,
    compton_form_factor_e_tilde: complex,
    compton_form_factor_h_conjugate: complex,
    compton_form_factor_h_tilde_conjugate: complex,
    compton_form_factor_e_conjugate: complex,
    compton_form_factor_e_tilde_conjugate: complex,
    verbose: bool = False) -> float:

    try:
        
        # (1): Calculate the appearance of Q^{2} + x_{B} t:
        sum_Q_squared_xb_t = squared_Q_momentum_transfer + x_Bjorken * squared_hadronic_momentum_transfer_t

        # (2): Calculate (2 - x_{B})Q^{2} + x_{B} t:
        weighted_sum_Q_squared_xb_t = (2. - x_Bjorken) * squared_Q_momentum_transfer + x_Bjorken * squared_hadronic_momentum_transfer_t

        # (3): Calculate Q^{2} (Q^{2} + x_{B} t):
        Q_squared_times_sum = squared_Q_momentum_transfer * sum_Q_squared_xb_t

        # (4): Calculate the first product of CFFs:
        cff_h_h_star_with_prefactor = compton_form_factor_h * compton_form_factor_h_conjugate * 4. * (1. - x_Bjorken)

        # (5): Calculate the second product of CFFs:
        cff_h_tilde_h_tilde_star = compton_form_factor_h_tilde * compton_form_factor_h_tilde_conjugate

        # (6): Calculate the third product of CFFs:
        cff_h_e_star_plus_e_h_star = compton_form_factor_h * compton_form_factor_e_conjugate + compton_form_factor_e * compton_form_factor_h_conjugate

        # (7): Calculate the fourth product of CFFs:
        cff_h_tilde_e_tilde_star_plus_e_tilde_h_tilde_star = compton_form_factor_h_tilde * compton_form_factor_e_tilde_conjugate + compton_form_factor_e_tilde * compton_form_factor_h_tilde_conjugate
        
        # (8): Calculate the fifth product of CFFs:
        cff_e_e_star = compton_form_factor_e * compton_form_factor_e_conjugate
        
        # (9): Calculate the sixth product of CFFs:
        cff_e_tilde_e_tilde_star = compton_form_factor_e_tilde * compton_form_factor_e_tilde_conjugate

        # (10): Calculate the second bracket term:
        second_bracket_term = 4. * (1. - x_Bjorken + ((2. * squared_Q_momentum_transfer + squared_hadronic_momentum_transfer_t) * epsilon**2 / (4. * sum_Q_squared_xb_t))) * cff_h_tilde_h_tilde_star

        # (11): Calculate the third_bracket term's prefactor
        third_bracket_term_prefactor = x_Bjorken**2 * (squared_Q_momentum_transfer + squared_hadronic_momentum_transfer_t)**2 / Q_squared_times_sum

        # (12): Calculate the fourth bracket term (yes, we're skipping the third for a minute):
        fourth_bracket_term = x_Bjorken**2 * squared_Q_momentum_transfer * cff_h_tilde_e_tilde_star_plus_e_tilde_h_tilde_star / sum_Q_squared_xb_t

        # (13): Calculate the fifth bracket term:
        fifth_bracket_term = (weighted_sum_Q_squared_xb_t**2 * squared_hadronic_momentum_transfer_t / (4. * _MASS_OF_PROTON_IN_GEV**2 * Q_squared_times_sum) + third_bracket_term_prefactor) * cff_e_e_star

        # (14): Calculate the third bracket term:
        third_bracket_term  = third_bracket_term_prefactor * cff_h_e_star_plus_e_h_star

        # (15): Calculate the sixth bracket term:
        sixth_bracket_term = x_Bjorken**2 * squared_Q_momentum_transfer * squared_hadronic_momentum_transfer_t * cff_e_tilde_e_tilde_star / (4. * _MASS_OF_PROTON_IN_GEV**2 * sum_Q_squared_xb_t)

        # (16): Return the entire thing:
        curlyC_unp_DVCS = Q_squared_times_sum * (cff_h_h_star_with_prefactor + second_bracket_term - third_bracket_term - fourth_bracket_term - fifth_bracket_term - sixth_bracket_term) / weighted_sum_Q_squared_xb_t**2
        
        # (13.1): If verbose, log the output:
        if verbose:
            print(f"> Calculated curlyC_unp_DVCS to be:\n{curlyC_unp_DVCS}")

        # (14): Return the coefficient:
        return curlyC_unp_DVCS

    except Exception as ERROR:
        print(f"> Error in calculating curlyC_unp_DVCS for DVCS Amplitude Squared:\n> {ERROR}")
        return 0.