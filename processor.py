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





def parse_arguments():
    parser = argparse.ArgumentParser(description="""The ccs_find project adds a GUI (graphical user interface) which uses input of mzML input files and molecular formulae to provide a unified set of results within a single data processing step which includes filtering using isotopic confirmation after peak picking.
""",
                        formatter_class=argparse.RawDescriptionHelpFormatter,
                        prog= "python processor.py")

    parser.add_argument("-g", dest = "config_file",
                        metavar = "config .ini file with parameters",
                        default = None,
                        required = True,
                        help = """config .ini file with all the parameters (required)""")
    args = parser.parse_args()
    #
    # try:
    if args.config_file:
        dictall = {}
        config = configparser.ConfigParser()
        # config.read('config.ini')
        config.read(args.config_file)
        cc = config.sections()
        optional_data = defaultdict(lambda: None)
        for c in cc:
            dictall.update(config[c])
        print(dictall)
        molecular_formula = dictall['molecular_formula']
        molecular_formula = [m for m in molecular_formula.split()]
        molecular_formula = [m.strip() for m in molecular_formula]
        xml_file = dictall['mzmlfile_path']
        primary_ion = dictall['primary_ion']
        drift_gas = dictall['gas']
        beta = dictall['beta']
        tfix = dictall.get("tfix")
        path_text = dictall['output_path']
        database_path = dictall['database_path']
        mono_combobox = dictall['mono_combobox']
        c13_combobox = dictall['c13_combobox']
        abundance_combobox = dictall['abundance_combobox']
        string_list = dictall["checked_ions"]
        string_list = [x for x in string_list.split()]
        _temp = tempfile.TemporaryDirectory(prefix = "drift_time_")
        drift = pathlib.Path(_temp.name).as_posix()
        temp = tempfile.TemporaryDirectory(prefix = "rt_isotopic_")
        spectrum = pathlib.Path(temp.name).as_posix()
        optional_data = defaultdict(lambda: None)
        data = defaultdict(lambda: None)
        message = defaultdict(list)
        import src
        secondary_data = defaultdict(list)
        parameters = [abundance_combobox, c13_combobox, mono_combobox]
        print(parameters)
        secondary_data["checked_ions"] = string_list
        secondary_data["ppm_values"] = parameters
        primary_data = {"primary_ion": primary_ion, "drift_gas": drift_gas, "mzml": xml_file, "beta": beta, "tfix": tfix, "buffer_text": molecular_formula, "chargestate": '', "ion_intensity": '', 'use_data': '', "peakwidth": 0.5}
        launcher = src.Final(primary_data, secondary_data, optional_data, mass_values, function_values, drift, spectrum, message, data)
        df22 = launcher.run_commandline()
        ccs_table = df22
        ccs_table["drift_time"] = ccs_table["drift_time"].map(lambda x: float(x)).map(lambda x: round(x, 3))
        ccs_table["retention_time"] = ccs_table["retention_time"].map(lambda x: float(x)).map(lambda x: round(x, 3))
        ccs_table.to_csv(os.path.join(path_text, 'Results.tab'), sep ='\t')
        import shutil
        output = pathlib.Path(drift).as_posix()
        output1 = pathlib.Path(spectrum).as_posix()
        shutil.copytree(output, os.path.join(path_text, "drift"))
        shutil.copytree(output1, os.path.join(path_text, "rt_isotope"))
        
    # except:
    #     print("Please add config file path for all the required information")




if __name__=="__main__":
    import pathlib
    if len(sys.argv) == 1:
        sys.argv.append("-h")
    args = parse_arguments()
