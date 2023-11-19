import argparse
import io
import sys, os
import src
import tempfile
import configparser
import numpy as np
import tempfile
from collections import defaultdict

p_mass = {"[M+H]+": 1.007276, "[M+]": 0.000548579909, "[M-H]-": 1.007276, "[M-]": 0.000548579909}


theoretical_isotope = {"+electron": '', "-electron": '', "+Cl":  "Cl",  "+Br": "Br", "+CH3COO": "CH3COO", "+HCOO": "HCOO", "+H": "H", "-H": '', "+Na": 'Na', "+K": 'K', "+NH4": 'NH4'}

mass_values = {"+electron": 0.000548579909, "-electron": 0.000548579909, "+Cl":  34.969402,  "+Br": 78.918885, "+CH3COO": 59.013851, "+HCOO": 44.998201, "+H": 1.007276, "-H": 1.007276, "+Na": 22.989218, "+K": 38.963158, "+NH4": 18.033823}

function_values = {"+electron": np.subtract, "-electron": np.add, "+Cl": np.add,  "+Br": np.add, "+CH3COO": np.add, "+HCOO": np.add, "+H": np.add, "-H": np.subtract, "+Na": np.add, "+K": np.add, "+NH4": np.add}


if __name__=="__main__":
    import pathlib
    try:
        _temp = tempfile.TemporaryDirectory(prefix = "drift_time_")
        drift = pathlib.Path(_temp.name).as_posix()
        temp = tempfile.TemporaryDirectory(prefix = "rt_isotopic_")
        spectrum = pathlib.Path(temp.name).as_posix()
        message = defaultdict(list)
        formula_list = ['C36H60N6O7S', 'C37H62N6O7S', 'C36H60N6O7S', 'C37H64N6O8']
        primary_data = {"primary_ion": "[M+H]+", "drift_gas": "Nitrogen", "mzml": "C:\\Users\kumar\\Downloads\\CCSfind_Datafiles\\20220808_Mutanobactin_LC-IM-TOFMS_Cal6.d.DeMP.mzML", "beta": 0.13581225201424416, "tfix": 0.428379630615936,  "buffer_text": formula_list}
        secondary_data = {'checked_ions': ['+Na', '+NH4', '+K'], 'mono_combobox': 8, 'c13_combobox': 8, 'abundance_combobox': 15}
        optional_data = defaultdict(lambda: None)

        mapper = src.Final(primary_data, secondary_data, optional_data, mass_values, function_values, drift, spectrum, message)
        mapper.run()

    except Exception as e:
        print(e)
