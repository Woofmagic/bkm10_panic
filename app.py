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
value_of_lepton_helicity = 0.5

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

use_WW = True
verbose = True


# These numbers are from Mathematica
# epsilon = 0.47293561004973345
# lepton_energy_fraction_y = 0.49609612355928445
# skewness_parameter = 0.19906188837146524
# squared_hadronic_momentum_transfer_t_minimum = -0.13551824472915253
# t_prime = -0.034481755270847486
# k_tilde = 0.1592415651944438
# shorthand_k = 0.08492693191323883

#  These numbers are from Python:
epsilon = 0.4729356100497334
lepton_energy_fraction_y = 0.4960961235592845
skewness_parameter = 0.1990618883714652
squared_hadronic_momentum_transfer_t_minimum = -0.13551824472915242
t_prime = -0.0344817552708476
k_tilde = 0.15326511787351363
shorthand_k = 0.08173956475763514
Dirac_form_factor_F1 = 0.7049508167585219
Pauli_form_factor_F2 = 1.1137103937669304

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

from unpolarized_curly_C import calculate_curly_C_unpolarized_interference
from unpolarized_curly_CV import calculate_curly_C_unpolarized_interference_V
from unpolarized_curly_CA import calculate_curly_C_unpolarized_interference_A

from n_equals_0_c_series.unpolarized_Cpp0 import calculate_c_0_plus_plus_unpolarized
from n_equals_0_c_series.unpolarized_Cpp0V import calculate_c_0_plus_plus_unpolarized_V
from n_equals_0_c_series.unpolarized_Cpp0A import calculate_c_0_plus_plus_unpolarized_A

from n_equals_0_c_series.unpolarized_C0p0 import calculate_c_0_zero_plus_unpolarized
from n_equals_0_c_series.unpolarized_C0p0V import calculate_c_0_zero_plus_unpolarized_V
from n_equals_0_c_series.unpolarized_C0p0A import calculate_c_0_zero_plus_unpolarized_A

from n_equals_1_c_series.unpolarized_Cpp1 import calculate_c_1_plus_plus_unpolarized
from n_equals_1_c_series.unpolarized_Cpp1V import calculate_c_1_plus_plus_unpolarized_V
from n_equals_1_c_series.unpolarized_Cpp1A import calculate_c_1_plus_plus_unpolarized_A

from n_equals_1_c_series.unpolarized_C0p1 import calculate_c_1_zero_plus_unpolarized
from n_equals_1_c_series.unpolarized_C0p1V import calculate_c_1_zero_plus_unpolarized_V
from n_equals_1_c_series.unpolarized_C0p1A import calculate_c_1_zero_plus_unpolarized_A

from n_equals_2_c_series.unpolarized_Cpp2 import calculate_c_2_plus_plus_unpolarized
from n_equals_2_c_series.unpolarized_Cpp2V import calculate_c_2_plus_plus_unpolarized_V
from n_equals_2_c_series.unpolarized_Cpp2A import calculate_c_2_plus_plus_unpolarized_A

from n_equals_2_c_series.unpolarized_C0p2 import calculate_c_2_zero_plus_unpolarized
from n_equals_2_c_series.unpolarized_C0p2V import calculate_c_2_zero_plus_unpolarized_V
from n_equals_2_c_series.unpolarized_C0p2A import calculate_c_2_zero_plus_unpolarized_A

from n_equals_3_c_series.unpolarized_Cpp3 import calculate_c_3_plus_plus_unpolarized
from n_equals_3_c_series.unpolarized_Cpp3V import calculate_c_3_plus_plus_unpolarized_V
from n_equals_3_c_series.unpolarized_Cpp3A import calculate_c_3_plus_plus_unpolarized_A

from n_equals_1_s_series.unpolarized_Spp1 import calculate_s_1_plus_plus_unpolarized
from n_equals_1_s_series.unpolarized_Spp1V import calculate_s_1_plus_plus_unpolarized_V
from n_equals_1_s_series.unpolarized_Spp1A import calculate_s_1_plus_plus_unpolarized_A

from n_equals_1_s_series.unpolarized_S0p1 import calculate_s_1_zero_plus_unpolarized
from n_equals_1_s_series.unpolarized_S0p1V import calculate_s_1_zero_plus_unpolarized_V
from n_equals_1_s_series.unpolarized_S0p1A import calculate_s_1_zero_plus_unpolarized_A

from n_equals_2_s_series.unpolarized_Spp2 import calculate_s_2_plus_plus_unpolarized
from n_equals_2_s_series.unpolarized_Spp2V import calculate_s_2_plus_plus_unpolarized_V
from n_equals_2_s_series.unpolarized_Spp2A import calculate_s_2_plus_plus_unpolarized_A

from n_equals_2_s_series.unpolarized_S0p2 import calculate_s_2_zero_plus_unpolarized
from n_equals_2_s_series.unpolarized_S0p2V import calculate_s_2_zero_plus_unpolarized_V
from n_equals_2_s_series.unpolarized_S0p2A import calculate_s_2_zero_plus_unpolarized_A

C0_pp_unpolarized = calculate_c_0_plus_plus_unpolarized(
    value_of_Q_squared,
    value_of_x_Bjorken,
    value_of_hadron_recoil,
    epsilon,
    lepton_energy_fraction_y,
    k_tilde,
    verbose)

C0V_pp_unpolarized = calculate_c_0_plus_plus_unpolarized_V(
    value_of_Q_squared,
    value_of_x_Bjorken,
    value_of_hadron_recoil,
    epsilon,
    lepton_energy_fraction_y,
    k_tilde,
    verbose)

C0A_pp_unpolarized = calculate_c_0_plus_plus_unpolarized_A(
    value_of_Q_squared,
    value_of_x_Bjorken,
    value_of_hadron_recoil,
    epsilon,
    lepton_energy_fraction_y,
    k_tilde,
    verbose)

C0_0p_unpolarized = calculate_c_0_zero_plus_unpolarized(
    value_of_Q_squared,
    value_of_x_Bjorken,
    value_of_hadron_recoil,
    epsilon,
    lepton_energy_fraction_y,
    shorthand_k,
    verbose)

C0V_0p_unpolarized = calculate_c_0_zero_plus_unpolarized_V(
    value_of_Q_squared,
    value_of_x_Bjorken,
    value_of_hadron_recoil,
    epsilon,
    lepton_energy_fraction_y,
    shorthand_k,
    verbose)

C0A_0p_unpolarized = calculate_c_0_zero_plus_unpolarized_A(
    value_of_Q_squared,
    value_of_x_Bjorken,
    value_of_hadron_recoil,
    epsilon,
    lepton_energy_fraction_y,
    shorthand_k,
    verbose)

curly_C_unpolarized_interference_for_pp = calculate_curly_C_unpolarized_interference(
    value_of_Q_squared,
    value_of_x_Bjorken,
    value_of_hadron_recoil,
    Dirac_form_factor_F1,
    Pauli_form_factor_F2,
    compton_form_factor_h,
    compton_form_factor_h_tilde,
    compton_form_factor_e,
    verbose)

curly_C_V_unpolarized_interference_for_pp = calculate_curly_C_unpolarized_interference_V(
    value_of_Q_squared,
    value_of_x_Bjorken,
    value_of_hadron_recoil,
    Dirac_form_factor_F1,
    Pauli_form_factor_F2,
    compton_form_factor_h,
    compton_form_factor_e,
    verbose)

curly_C_A_unpolarized_interference_for_pp = calculate_curly_C_unpolarized_interference_A(
    value_of_Q_squared,
    value_of_x_Bjorken,
    value_of_hadron_recoil,
    Dirac_form_factor_F1,
    Pauli_form_factor_F2,
    compton_form_factor_h_tilde,
    verbose)

curly_C_unpolarized_interference_for_0p = calculate_curly_C_unpolarized_interference(
    value_of_Q_squared,
    value_of_x_Bjorken,
    value_of_hadron_recoil,
    Dirac_form_factor_F1,
    Pauli_form_factor_F2,
    compute_cff_effective(skewness_parameter, compton_form_factor_h, use_WW),
    compute_cff_effective(skewness_parameter, compton_form_factor_h_tilde, use_WW),
    compute_cff_effective(skewness_parameter, compton_form_factor_e, use_WW),
    verbose)

curly_C_V_unpolarized_interference_for_0p = calculate_curly_C_unpolarized_interference_V(
    value_of_Q_squared,
    value_of_x_Bjorken,
    value_of_hadron_recoil,
    Dirac_form_factor_F1,
    Pauli_form_factor_F2,
    compute_cff_effective(skewness_parameter, compton_form_factor_h, use_WW),
    compute_cff_effective(skewness_parameter, compton_form_factor_e, use_WW),
    verbose)

curly_C_A_unpolarized_interference_for_0p = calculate_curly_C_unpolarized_interference_A(
    value_of_Q_squared,
    value_of_x_Bjorken,
    value_of_hadron_recoil,
    Dirac_form_factor_F1,
    Pauli_form_factor_F2,
    compute_cff_effective(skewness_parameter, compton_form_factor_h_tilde, use_WW),
    verbose)

curly_C_0_pp_int =  (curly_C_unpolarized_interference_for_pp
            + (C0V_pp_unpolarized * curly_C_V_unpolarized_interference_for_pp / C0_pp_unpolarized)
            + (C0A_pp_unpolarized * curly_C_A_unpolarized_interference_for_pp / C0_pp_unpolarized))

curly_C_0_0p_int =  ((np.sqrt(2. / value_of_Q_squared) * k_tilde / (2. - value_of_x_Bjorken)) * (curly_C_unpolarized_interference_for_0p
            + (C0V_0p_unpolarized * curly_C_V_unpolarized_interference_for_0p / C0_0p_unpolarized)
            + (C0A_0p_unpolarized * curly_C_A_unpolarized_interference_for_0p / C0_0p_unpolarized)))

C1_pp_unpolarized = calculate_c_1_plus_plus_unpolarized(
    value_of_Q_squared,
    value_of_x_Bjorken,
    value_of_hadron_recoil,
    epsilon,
    lepton_energy_fraction_y,
    shorthand_k,
    verbose)

C1V_pp_unpolarized = calculate_c_1_plus_plus_unpolarized_V(
    value_of_Q_squared,
    value_of_x_Bjorken,
    value_of_hadron_recoil,
    epsilon,
    lepton_energy_fraction_y,
    t_prime,
    shorthand_k,
    verbose)

C1A_pp_unpolarized = calculate_c_1_plus_plus_unpolarized_A(
    value_of_Q_squared,
    value_of_x_Bjorken,
    value_of_hadron_recoil,
    epsilon,
    lepton_energy_fraction_y,
    t_prime,
    shorthand_k,
    verbose)

C1_0p_unpolarized = calculate_c_1_zero_plus_unpolarized(
    value_of_Q_squared,
    value_of_x_Bjorken,
    value_of_hadron_recoil,
    epsilon,
    lepton_energy_fraction_y,
    t_prime,
    verbose)

C1V_0p_unpolarized = calculate_c_1_zero_plus_unpolarized_V(
    value_of_Q_squared,
    value_of_x_Bjorken,
    value_of_hadron_recoil,
    epsilon,
    lepton_energy_fraction_y,
    k_tilde,
    verbose)

C1A_0p_unpolarized = calculate_c_1_zero_plus_unpolarized_A(
    value_of_Q_squared,
    value_of_x_Bjorken,
    value_of_hadron_recoil,
    epsilon,
    lepton_energy_fraction_y,
    k_tilde,
    verbose)

curly_C_1_pp_int =  (curly_C_unpolarized_interference_for_pp
            + (C1V_pp_unpolarized * curly_C_V_unpolarized_interference_for_pp / C1_pp_unpolarized)
            + (C1A_pp_unpolarized * curly_C_A_unpolarized_interference_for_pp / C1_pp_unpolarized))

curly_C_1_0p_int =  ((np.sqrt(2. / value_of_Q_squared) * k_tilde / (2. - value_of_x_Bjorken)) * (curly_C_unpolarized_interference_for_0p
            + (C1V_0p_unpolarized * curly_C_V_unpolarized_interference_for_0p / C1_0p_unpolarized)
            + (C1A_0p_unpolarized * curly_C_A_unpolarized_interference_for_0p / C1_0p_unpolarized)))

C2_pp_unpolarized = calculate_c_2_plus_plus_unpolarized(
    value_of_Q_squared,
    value_of_x_Bjorken,
    value_of_hadron_recoil,
    epsilon,
    lepton_energy_fraction_y,
    t_prime,
    k_tilde,
    verbose)

C2V_pp_unpolarized = calculate_c_2_plus_plus_unpolarized_V(
    value_of_Q_squared,
    value_of_x_Bjorken,
    value_of_hadron_recoil,
    epsilon,
    lepton_energy_fraction_y,
    t_prime,
    k_tilde,
    verbose)

C2A_pp_unpolarized = calculate_c_2_plus_plus_unpolarized_A(
    value_of_Q_squared,
    value_of_x_Bjorken,
    value_of_hadron_recoil,
    epsilon,
    lepton_energy_fraction_y,
    t_prime,
    k_tilde,
    verbose)

C2_0p_unpolarized = calculate_c_2_zero_plus_unpolarized(
    value_of_Q_squared,
    value_of_x_Bjorken,
    value_of_hadron_recoil,
    epsilon,
    lepton_energy_fraction_y,
    shorthand_k,
    verbose)

C2V_0p_unpolarized = calculate_c_2_zero_plus_unpolarized_V(
    value_of_Q_squared,
    value_of_x_Bjorken,
    value_of_hadron_recoil,
    epsilon,
    lepton_energy_fraction_y,
    shorthand_k,
    verbose)

C2A_0p_unpolarized = calculate_c_2_zero_plus_unpolarized_A(
    value_of_Q_squared,
    value_of_x_Bjorken,
    value_of_hadron_recoil,
    epsilon,
    lepton_energy_fraction_y,
    t_prime,
    shorthand_k,
    verbose)

curly_C_2_pp_int =  (curly_C_unpolarized_interference_for_pp
            + (C2V_pp_unpolarized * curly_C_V_unpolarized_interference_for_pp / C2_pp_unpolarized)
            + (C2A_pp_unpolarized * curly_C_A_unpolarized_interference_for_pp / C2_pp_unpolarized))

curly_C_2_0p_int =  ((np.sqrt(2. / value_of_Q_squared) * k_tilde / (2. - value_of_x_Bjorken)) * (curly_C_unpolarized_interference_for_0p
            + (C2V_0p_unpolarized * curly_C_V_unpolarized_interference_for_0p / C2_0p_unpolarized)
            + (C2A_0p_unpolarized * curly_C_A_unpolarized_interference_for_0p / C2_0p_unpolarized)))

C3_pp_unpolarized = calculate_c_3_plus_plus_unpolarized(
    value_of_Q_squared,
    value_of_x_Bjorken,
    value_of_hadron_recoil,
    epsilon,
    lepton_energy_fraction_y,
    shorthand_k,
    verbose)

C3V_pp_unpolarized = calculate_c_3_plus_plus_unpolarized_V(
    value_of_Q_squared,
    value_of_x_Bjorken,
    value_of_hadron_recoil,
    epsilon,
    lepton_energy_fraction_y,
    shorthand_k,
    verbose)

C3A_pp_unpolarized = calculate_c_3_plus_plus_unpolarized_A(
    value_of_Q_squared,
    value_of_x_Bjorken,
    value_of_hadron_recoil,
    epsilon,
    lepton_energy_fraction_y,
    t_prime,
    shorthand_k,
    verbose)

C3_0p_unpolarized = 0.

C3V_0p_unpolarized = 0.

C3A_0p_unpolarized = 0.

curly_C_3_pp_int =  (curly_C_unpolarized_interference_for_pp
            + (C3V_pp_unpolarized * curly_C_V_unpolarized_interference_for_pp / C3_pp_unpolarized)
            + (C3A_pp_unpolarized * curly_C_A_unpolarized_interference_for_pp / C3_pp_unpolarized))

curly_C_3_0p_int =  ((np.sqrt(2. / value_of_Q_squared) * k_tilde / (2. - value_of_x_Bjorken)) * (curly_C_unpolarized_interference_for_0p))

S1_pp_unpolarized = calculate_s_1_plus_plus_unpolarized(
    value_of_lepton_helicity,
    value_of_Q_squared,
    value_of_x_Bjorken,
    epsilon,
    lepton_energy_fraction_y,
    t_prime,
    shorthand_k,
    verbose)

S1V_pp_unpolarized = calculate_s_1_plus_plus_unpolarized_V(
    value_of_lepton_helicity,
    value_of_Q_squared,
    value_of_x_Bjorken,
    value_of_hadron_recoil,
    epsilon,
    lepton_energy_fraction_y,
    shorthand_k,
    verbose)

S1A_pp_unpolarized = calculate_s_1_plus_plus_unpolarized_A(
    value_of_lepton_helicity,
    value_of_Q_squared,
    value_of_x_Bjorken,
    value_of_hadron_recoil,
    epsilon,
    lepton_energy_fraction_y,
    t_prime,
    shorthand_k,
    verbose)

S1_0p_unpolarized = calculate_s_1_zero_plus_unpolarized(
    value_of_lepton_helicity,
    value_of_Q_squared,
    epsilon,
    lepton_energy_fraction_y,
    k_tilde,
    verbose)

S1V_0p_unpolarized = calculate_s_1_zero_plus_unpolarized_V(
    value_of_lepton_helicity,
    value_of_Q_squared,
    value_of_x_Bjorken,
    value_of_hadron_recoil,
    epsilon,
    lepton_energy_fraction_y,
    verbose)

S1A_0p_unpolarized = calculate_s_1_zero_plus_unpolarized_A(
    value_of_lepton_helicity,
    value_of_Q_squared,
    value_of_x_Bjorken,
    value_of_hadron_recoil,
    epsilon,
    lepton_energy_fraction_y,
    shorthand_k,
    verbose)

curly_S_1_pp_int =  (curly_C_unpolarized_interference_for_pp
            + (S1V_pp_unpolarized * curly_C_V_unpolarized_interference_for_pp / S1_pp_unpolarized)
            + (S1A_pp_unpolarized * curly_C_A_unpolarized_interference_for_pp / S1_pp_unpolarized))

print(f"> Curly C: {curly_C_unpolarized_interference_for_pp}")
print(f"> Curly C,V: {curly_C_V_unpolarized_interference_for_pp}")
print(f"> Curly C,A: {curly_C_A_unpolarized_interference_for_pp}")
print(f"> Curly S(n=2)pp: {S1_pp_unpolarized}")
print(f"> Curly S,V(n=2)pp: {S1V_pp_unpolarized}")
print(f"> Curly S,A(n=2)pp: {S1A_pp_unpolarized}")

print(f"> Curly S++: {curly_S_1_pp_int}")

curly_S_1_0p_int =  ((np.sqrt(2. / value_of_Q_squared) * k_tilde / (2. - value_of_x_Bjorken)) * (curly_C_unpolarized_interference_for_0p
            + (S1V_0p_unpolarized * curly_C_V_unpolarized_interference_for_0p / S1_0p_unpolarized)
            + (S1A_0p_unpolarized * curly_C_A_unpolarized_interference_for_0p / S1_0p_unpolarized)))

print(f"> Curly S0+: {curly_S_1_0p_int}")

S2_pp_unpolarized = calculate_s_2_plus_plus_unpolarized(
    value_of_lepton_helicity,
    value_of_Q_squared,
    value_of_x_Bjorken,
    epsilon,
    lepton_energy_fraction_y,
    t_prime,
    verbose)

S2V_pp_unpolarized = calculate_s_2_plus_plus_unpolarized_V(
    value_of_lepton_helicity,
    value_of_Q_squared,
    value_of_x_Bjorken,
    value_of_hadron_recoil,
    epsilon,
    lepton_energy_fraction_y,
    verbose)

S2A_pp_unpolarized = calculate_s_2_plus_plus_unpolarized_A(
    value_of_lepton_helicity,
    value_of_Q_squared,
    value_of_x_Bjorken,
    value_of_hadron_recoil,
    epsilon,
    lepton_energy_fraction_y,
    t_prime,
    verbose)

S2_0p_unpolarized = calculate_s_2_zero_plus_unpolarized(
    value_of_lepton_helicity,
    value_of_Q_squared,
    value_of_x_Bjorken,
    value_of_hadron_recoil,
    epsilon,
    lepton_energy_fraction_y,
    shorthand_k,
    verbose)

S2V_0p_unpolarized = calculate_s_2_zero_plus_unpolarized_V(
    value_of_lepton_helicity,
    value_of_Q_squared,
    value_of_x_Bjorken,
    value_of_hadron_recoil,
    epsilon,
    lepton_energy_fraction_y,
    shorthand_k,
    verbose)

S2A_0p_unpolarized = calculate_s_2_zero_plus_unpolarized_A(
    value_of_lepton_helicity,
    value_of_Q_squared,
    value_of_x_Bjorken,
    value_of_hadron_recoil,
    epsilon,
    lepton_energy_fraction_y,
    shorthand_k,
    verbose)

curly_S_2_pp_int =  (curly_C_unpolarized_interference_for_pp
            + (S2V_pp_unpolarized * curly_C_V_unpolarized_interference_for_pp / S2_pp_unpolarized)
            + (S2A_pp_unpolarized * curly_C_A_unpolarized_interference_for_pp / S2_pp_unpolarized))

print(f"> Curly C: {curly_C_unpolarized_interference_for_pp}")
print(f"> Curly C,V: {curly_C_V_unpolarized_interference_for_pp}")
print(f"> Curly C,A: {curly_C_A_unpolarized_interference_for_pp}")
print(f"> Curly S(n=2)pp: {S2_pp_unpolarized}")
print(f"> Curly S,V(n=2)pp: {S2V_pp_unpolarized}")
print(f"> Curly S,A(n=2)pp: {S2A_pp_unpolarized}")

print(f"> Curly S++: {curly_S_2_pp_int}")

curly_S_2_0p_int =  ((np.sqrt(2. / value_of_Q_squared) * k_tilde / (2. - value_of_x_Bjorken)) * (curly_C_unpolarized_interference_for_0p
            + (S2V_0p_unpolarized * curly_C_V_unpolarized_interference_for_0p / S2_0p_unpolarized)
            + (S2A_0p_unpolarized * curly_C_A_unpolarized_interference_for_0p / S2_0p_unpolarized)))

print(f"> Curly S0+: {curly_S_2_0p_int}")

c_0_interference_coefficient = C0_pp_unpolarized * curly_C_0_pp_int.real + C0_0p_unpolarized * curly_C_0_0p_int.real
c_1_interference_coefficient = C1_pp_unpolarized * curly_C_1_pp_int.real + C1_0p_unpolarized * curly_C_1_0p_int.real
c_2_interference_coefficient = C2_pp_unpolarized * curly_C_2_pp_int.real + C2_0p_unpolarized * curly_C_2_0p_int.real
c_3_interference_coefficient = C3_pp_unpolarized * curly_C_3_pp_int.real + C3_0p_unpolarized * curly_C_3_0p_int.real
s_1_interference_coefficient = S1_pp_unpolarized * curly_S_1_pp_int.imag + S1_0p_unpolarized * curly_S_1_0p_int.imag
s_2_interference_coefficient = S2_pp_unpolarized * curly_S_2_pp_int.imag + S2_0p_unpolarized * curly_S_2_0p_int.imag

print(f"> c_0 = {c_0_interference_coefficient}")
print(f"> c_1 = {c_1_interference_coefficient}")
print(f"> c_2 = {c_2_interference_coefficient}")
print(f"> c_3 = {c_3_interference_coefficient}")
print(f"> s_1 = {s_1_interference_coefficient}")
print(f"> s_2 = {s_2_interference_coefficient}")