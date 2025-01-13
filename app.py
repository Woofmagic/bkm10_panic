import numpy as np
import pandas as pd

_MASS_OF_PROTON_IN_GEV = .93827208816
_MASS_OF_PROTON_SQUARED_IN_GEV_SQUARED = _MASS_OF_PROTON_IN_GEV *_MASS_OF_PROTON_IN_GEV
_ELECTROMAGNETIC_FINE_STRUCTURE_CONSTANT = (1.0 / 137.0359998)
_ELECTRIC_FORM_FACTOR_CONSTANT = 0.710649
_PROTON_MAGNETIC_MOMENT = 2.79284734463

value_of_Q_squared = 1.82
value_of_x_Bjorken = 0.34
value_of_hadron_recoil = -0.17
value_of_beam_energy = 5.75

compton_form_factor_h_real = -0.897
compton_form_factor_h_imaginary = 2.421
compton_form_factor_e_real = -0.541
compton_form_factor_e_imaginary = 0.903

compton_form_factor_h_tilde_real = 2.444
compton_form_factor_h_tilde_imaginary = 1.131
compton_form_factor_e_tilde_real = 2.207
compton_form_factor_e_tilde_imaginary = 5.383

compton_form_factor_h = complex(compton_form_factor_h_real, compton_form_factor_h_imaginary)
compton_form_factor_h_tilde = complex(compton_form_factor_h_tilde_real, compton_form_factor_h_tilde_imaginary)
compton_form_factor_e = complex(compton_form_factor_e_real, compton_form_factor_e_imaginary)
compton_form_factor_e_tilde = complex(compton_form_factor_e_tilde_real, compton_form_factor_e_tilde_imaginary)

# squared_Q_momentum_transfer = np.array([value_of_Q_squared for _ in range(len(np.arange(0, 361, 1.)))])
# x_Bjorken = np.array([value_of_x_Bjorken for _ in range(len(np.arange(0, 361, 1.)))])
# squared_hadronic_momentum_transfer_t = np.array([value_of_hadron_recoil for _ in range(len(np.arange(0, 361, 1.)))])
# lab_kinematics_k = np.array([value_of_beam_energy for _ in range(len(np.arange(0, 361, 1.)))])
# azimuthal_phi = np.array([phi for phi in range(len(np.arange(0., 361., 1.)))]) 
verbose = False

# These numbers are from Mathematica
epsilon = 0.47293561004973345
lepton_energy_fraction_y = 0.49609612355928445
skewness_parameter = 0.19906188837146524
squared_hadronic_momentum_transfer_t_minimum = -0.13551824472915253
t_prime = -0.034481755270847486
k_tilde = 0.1592415651944438
shorthand_k = 0.08492693191323883

# (2.X): Calculate the cross-section prefactor:
cross_section_prefactor = 3.5309544777485675e-10



def compute_cff_effective(
    skewness_parameter: float,
    compton_form_factor: complex,
    use_ww: bool = False,
    verbose: bool = False) -> complex:
    try:
        if use_ww:
            cff_effective = 2. * compton_form_factor / (1. + skewness_parameter)
        else:
            cff_effective = -2. * skewness_parameter * compton_form_factor / (1. + skewness_parameter)
        if verbose:
            print(f"> Computed the CFF effective to be:\n{cff_effective}")
        return cff_effective
    except Exception as error:
        print(f"> Error in calculating F_effective:\n> {error}")
        return 0.

from unpolarized_curlyC_dvcs import calculate_curly_c_unpolarized_dvcs

curly_C_unpolarized_DVCS_without_WW = calculate_curly_c_unpolarized_dvcs(
    value_of_Q_squared,
    value_of_x_Bjorken,
    value_of_hadron_recoil,
    epsilon,
    compute_cff_effective(skewness_parameter, compton_form_factor_h, False),
    compute_cff_effective(skewness_parameter, compton_form_factor_h_tilde, False),
    compute_cff_effective(skewness_parameter, compton_form_factor_e, False),
    compute_cff_effective(skewness_parameter, compton_form_factor_e_tilde, False),
    compton_form_factor_h.conjugate(),
    compton_form_factor_h_tilde.conjugate(),
    compton_form_factor_e.conjugate(),
    compton_form_factor_e_tilde.conjugate(),
    verbose)

print(curly_C_unpolarized_DVCS_without_WW)

curly_C_unpolarized_DVCS_with_WW = calculate_curly_c_unpolarized_dvcs(
    value_of_Q_squared,
    value_of_x_Bjorken,
    value_of_hadron_recoil,
    epsilon,
    compute_cff_effective(skewness_parameter, compton_form_factor_h, True),
    compute_cff_effective(skewness_parameter, compton_form_factor_h_tilde, True),
    compute_cff_effective(skewness_parameter, compton_form_factor_e, True),
    compute_cff_effective(skewness_parameter, compton_form_factor_e_tilde, True),
    compton_form_factor_h.conjugate(),
    compton_form_factor_h_tilde.conjugate(),
    compton_form_factor_e.conjugate(),
    compton_form_factor_e_tilde.conjugate(),
    verbose)

print(curly_C_unpolarized_DVCS_with_WW)