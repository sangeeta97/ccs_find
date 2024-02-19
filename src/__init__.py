import numpy as np
import pandas as pd
from collections import defaultdict, OrderedDict
from .utils import timeit, timeit1




parameters = {"abundance_combobox": 17, "c13_combobox": 20, "mono_combobox": 20}


class Final:

    def __init__(self, primary_data, secondary_data, optional_data, mass_values, function_values, temp_drift, temp_spectrum, message, mzmldata):

        self.primary_data =  primary_data
        self.secondary_data = secondary_data
        self.optional_data = optional_data
        self.ion_species = False
        self.mass_values = mass_values
        self.function_values = function_values
        self.temp_drift = temp_drift
        self.temp_spectrum = temp_spectrum
        self.message = message
        self.ppm_values = list(parameters.values())
        self.isotope_match = defaultdict(lambda: defaultdict(lambda: None))
        self.all_rt = defaultdict(lambda: None)
        self.all_dt = defaultdict(lambda: None)
        self.mzmldata = mzmldata

    @timeit1
    def run(self):
        self.unpack_primary()
        if not self.message['warning']:
            self.unpack_secondary()
            from .formula import Molecular_Formula
            self.formula_df = Molecular_Formula(self)
            self.ions, self.isotope_ratio = self.formula_df.run()
            from .drift_time import Parse_MF
            self.drift_time_df = Parse_MF(self)
            df9 = self.drift_time_df.dataframe_all()
            if "Alert" in self.message.keys():
                return None
            from .ccs import CCS
            self.ccs_df = ccs.CCS(self, df9)
            df22 = self.ccs_df.finish_analysis()
            return df22


    @timeit1
    def run_commandline(self):
        self.ppm_values = self.secondary_data["ppm_values"]
        self.ppm_values = [float(x) for x in self.ppm_values]
        if bool(self.secondary_data.get("checked_ions", 0)):
            self.ion_species = True
            keys_to_search = set(self.secondary_data["checked_ions"])
            mass_set = set(list(self.mass_values.keys()))
            function_set = set(list(self.function_values.keys()))
            mass_set = mass_set.intersection(keys_to_search)
            self.mass_values= {k: self.mass_values[k] for k in mass_set}
            function_set = function_set.intersection(keys_to_search)
            self.function_values= {k: self.function_values[k] for k in function_set}
        from .formula import Molecular_Formula
        self.formula_df = Molecular_Formula(self)
        self.ions, self.isotope_ratio = self.formula_df.run()
        from .drift_time import Parse_MF
        self.drift_time_df = Parse_MF(self)
        df9 = self.drift_time_df.dataframe_all()
        if "Alert" in self.message.keys():
            return None
        from .ccs import CCS
        self.ccs_df = ccs.CCS(self, df9)
        df22 = self.ccs_df.finish_analysis()
        return df22






    def unpack_primary(self):
        xx = list(self.primary_data.values())
        tt = [bool(x) for x in xx]
        bool_value = np.all(tt)
        if not bool_value:
            bb = self.optional_data['formula']
            x = isinstance(bb, list)
            if x:
                self.primary_data['buffer_text'] = self.optional_data['formula']
            bb = self.optional_data['calibration']
            x = isinstance(bb, dict)
            if x:
                self.primary_data['tfix'] = bb['TFix']
                self.primary_data['beta'] = bb['Beta']
            keys_list = ['tfix', 'beta', 'buffer_text', 'drift_gas',  'mzml', 'primary_ion']
            xx = [self.primary_data[x] for x in keys_list]
            tt = [bool(x) for x in xx]
            bool_value = np.all(tt)
            if not bool_value:
                self.message['warning'].append("The compulsory fields should be filled, the mzml file,  primary_ion, tfix, beta value, (or .xml cal file should be uploaded) MF (or formula file should be upload), drift gas must be provided")






    def unpack_secondary(self):

        bool_value = bool(list(self.secondary_data.values()))
        if bool_value:
            self.secondary_data = OrderedDict(sorted(self.secondary_data.items()))

            self.ppm_values = [self.secondary_data.get(x, 0.0) if self.secondary_data.get(x, 0.0) != 0.0 else parameters.get(x, 0.0) for x in parameters.keys()]
            self.ppm_values = [float(x) for x in self.ppm_values]
            if bool(self.secondary_data.get("checked_ions", 0)):
                self.ion_species = True
                keys_to_search = set(self.secondary_data["checked_ions"])
                mass_set = set(list(self.mass_values.keys()))
                function_set = set(list(self.function_values.keys()))
                mass_set = mass_set.intersection(keys_to_search)
                self.mass_values= {k: self.mass_values[k] for k in mass_set}
                function_set = function_set.intersection(keys_to_search)
                self.function_values= {k: self.function_values[k] for k in function_set}
