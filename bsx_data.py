
normal_ranges = {
    "u": 1,
    "P_a": 100,
    "P_an": 100,
    "Pa_CO2": 40,
    "SaO2sup": 0.96,
    "CBFn": 0.0125,
    "CMRO2": 0.033953,
    "CMRO2_n": 0.033953,
    "Xtot": 91,
    "a": 0.06567,
    "a_n": 0.06567,
    "R_autc": 2.2,
    "bred": 0.001408,
    "bred_n": 0.001408,
    "C_im": 0.00675,
    "cytox_tot_tis": 0.0055,
    "Dpsi": 145,
    "Dpsi_n": 145,
    "H":  0.00003981,
    "H_n":  0.00003981,
    "k_aut": 1,
    "kCV": 0.02047339,
    "O2": 0.024,
    "O2_n": 0.024,
    "phi": 0.036,
    "R_auto": 1.5,
    "R_autp": 4,
    "R_autu": 0.5,
    "r_n": 0.0187,
    "cytox_tot_tis": 0.0055,
    "P_v": 4,
    "P_vn": 4

}
folder_param_dict = {
    f"{param_name}": [param_name] for param_name in normal_ranges.keys() if param_name!="u"
}

study_range_percentage = 0.25
study_ranges =  {
    **{k: ((1-study_range_percentage)*normal_ranges[k], (1+study_range_percentage)*normal_ranges[k]) \
        for k in normal_ranges.keys()}} 