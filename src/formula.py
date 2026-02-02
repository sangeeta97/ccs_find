import os
import json
import pandas as pd
import statistics
from scipy import interpolate
import numpy as np
from numpy.linalg import norm
import copy
from collections import defaultdict
import operator
from molmass import Formula
from collections import namedtuple
import io



def isotopic_table(molecular_formula):
    f = Formula(molecular_formula)
    dd = f.spectrum()
    tt = str(dd)
    gg = tt.split("\n")
    Ratio = namedtuple('Spectrum', ['mz', 'fraction', 'intensity'])
    yy = []
    for xx in gg[1: -1]:
      ll = [x for x in xx.split("  ") if bool(x)]
      S = Ratio._make(ll[1:])
      yy.append(S)
    df = pd.DataFrame(yy)
    df = df.astype(float)
    return df



p_mass = {"[M+H]+": 1.007276, "[M]+": 0.000548579909, "[M-H]-": 1.007276, "[M]-": 0.000548579909}



p_func = {"[M+H]+": np.add, "[M]+": np.subtract, "[M-H]-": np.subtract, "[M]-": np.add}


charge_two = {"[M+H]+": "[M+2H]2+", "[M-H]-": "[M-2H]2-"}



class Molecular_Formula:
    """
    Methods to compute molecular mass with isotopic distribution from a DataContainer
    Methods
    -------
    monomass: Computes the monoisotpoic mass of each molecular formula and its combination.
    Iso-ratio: computes the isotopic abundance theoretical ratio for each compound.
    df1 class attribute stores the all the combination requested by the user as primary ion and adduct including doubly charged ions with their corresponding monoisiotopic mass.
    """

    def __init__(self, data):
        self.data = data
        self.primary_ion = self.data.primary_data['primary_ion']



    def run(self):
        self.formula_define()
        x = self.table()
        y = self.isotopic_ratio()
        return x, y


    def formula_define(self):
        '''
        Create df1 class attribute which is dataframe and store the information about molecular formula and its mz value
         either obtained from text information or from a file
        '''
        if isinstance(self.data.primary_data["buffer_text"], list):
            data = {"molecular_formula": self.data.primary_data["buffer_text"]}
            self.df1 = pd.DataFrame(data = data)
            self.df1.columns = ["molecular_formula"]
            self.df1['mono'] = [Formula(x) for x in self.df1["molecular_formula"].values]
            self.df1['mono_mass'] = self.df1['mono'].map(lambda x: x.isotope.mass)

        elif isinstance(self.data.optional_data['formula'], list):
            data = {"molecular_formula": self.data.optional_data['formula']}
            self.df1 = pd.DataFrame(data = data)
            self.df1['mono'] = [Formula(x) for x in self.df1["molecular_formula"].values]
            self.df1['mono_mass'] = self.df1['mono'].map(lambda x: x.isotope.mass)

        else:
            self.data.message["warning"].append("formula input must be given either text or file")




    def table(self):
        '''
         adds new columns in the dataframe of formulas as index (named as df1) with adduct type as coulmn name and mass values for all formulas as column values
         and it returns the dictionary with the molecular formula as key and dataframe storing the ion types as index and mz value as column values.
        '''
        p_ion = self.data.primary_data["primary_ion"]
        self.df1[p_ion] = self.df1['mono_mass'].map(lambda x: p_func[p_ion](x, p_mass[p_ion]))

        if self.data.primary_data["chargestate"]:
            if self.data.primary_data["primary_ion"] in ["[M+H]+", "[M-H]-"]:
                doubly = charge_two[p_ion]
                self.df1[doubly] = self.df1[p_ion].map(lambda x: x + 1.007276).map(lambda x: x/2)

        if self.data.ion_species:
            for c in self.data.mass_values.keys():
                self.df1[c] = self.df1['mono_mass'].map(lambda x: self.data.function_values[c](x, self.data.mass_values[c]))

        if bool(self.data.secondary_data.get('positive_formula', 0)):
            for y in self.data.secondary_data['positive_formula']:
                y1 = y.strip("+")
                tt = Formula(y1)
                tt = tt.isotope.mass
                self.df1[y] = self.df1['mono_mass'].map(lambda x: x + tt)


        if bool(self.data.secondary_data.get('negative_formula', 0)):
            for z in self.data.secondary_data['negative_formula']:
                z1 = z.strip("-")
                cc = Formula(z1)
                cc = cc.isotope.mass
                self.df1[z] = self.df1['mono_mass'].map(lambda x: x - cc)
        df_list = []
        keys = []
        for gg, tt in enumerate(self.df1["molecular_formula"].values):
            header = self.df1['mono'].values
            all_gg = self.df1.iloc[gg, ].T
            all_gg = all_gg[3:]
            all_gg.columns = header[gg]
            df_list.append(all_gg)
            keys.append(tt)
        self.all_data = {g.strip("\n"):f for g, f in zip(keys, df_list)}
        return self.all_data



    def isotopic_ratio(self):
        list_isotope = {}
        dict_all = self.all_data
        for kk in dict_all.keys():
            df = isotopic_table(kk)
            list_isotope.update([(kk, df)])
        return list_isotope
