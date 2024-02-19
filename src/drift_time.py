import os
import json
import statistics
from scipy import interpolate
import numpy as np
from numpy.linalg import norm
import copy
from collections import defaultdict
import operator
import math
from collections import namedtuple
import plotly.express as px
from .peak import *
import tempfile
from functools import wraps
import time
import matplotlib.pyplot as plt
from scipy.signal import peak_widths
from .plots import Xic_Eic
import numpy as np
import plotly.express as px
from .utils import *
import pandas as pd
from . import calculation
from biosaur2.utils import MS1OnlyMzML
import polars as pl
from molmass import Formula
from .utils import timeit
import pyarrow





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


def condition_min(target, ppm):
    allowed = np.multiply(ppm, target)
    allowed = np.divide(allowed, 1000000)
    target = target - allowed
    return target




def condition_max(target, ppm):
    allowed = np.multiply(ppm, target)
    allowed = np.divide(allowed, 1000000)
    target = target + allowed
    return target




def test_spec(spec):
    if isinstance(spec, pymzml.spec.MS_Spectrum):
        if spec.ms_level == 1:
            return True
        else:
            return False
    else:
        return False




class Parse_MF:
    '''

    This class runs for one MF at a time, and that MF becomes the class attribute called molecular_formula. this attributes changes by a method in the class called mix_all.


    The main class which calls other classes for plotting and peak-picking and summarizes the results for view.
    This class initilizes dictionary to be used in other classes to store the results and outputs from the other classes.

    Run_all method parse the mzml file into polars dataframe.

    This tool uses the primary ion as pivotol point for scanning the other ions.

    The primary information for each studied molecular formula is dataframe containg the mz value with ppm error along with corresponding rt, dt, intensity values.

    dframe attribute stores the molecular formula as first key with ion as second key and corresponding dataframe.

    rt_found_ion, is dict stores MF+rt as string and as key whose value is list of found ion type.

    All_dt, all_rt is for database entry.

    experimental_isotope is for plotting the eic of isotopes

    temp_spectrum, temp_isotope stores the rsult files and viewed upon calling the tool

    '''
    def __init__(self, data):
        self.data = data
        self.monoPPM = self.data.ppm_values
        self.PPM1 = self.monoPPM[2]
        self.PPM2 = self.monoPPM[1]
        self.mass_df = self.data.ions
        self.primary_ion = self.data.primary_data['primary_ion']
        self.peakwidth = self.data.primary_data['peakwidth']
        self.full_result = defaultdict(lambda: defaultdict(lambda: None))
        self.dframe = defaultdict(lambda: defaultdict(lambda: None))
        self.selected_ions = defaultdict(lambda: None)
        self.ion_dict = {x:y.values[0] for x, y in self.mass_df.items()}
        self.rt_found_ion = defaultdict(list)
        self.ions_data = defaultdict(lambda: defaultdict(lambda: None))
        self.found_molecular_formula = defaultdict(lambda: None)
        self.dt_peaks = defaultdict(lambda: defaultdict(lambda: None))
        self.isotope = self.data.isotope_ratio
        self.isotope_match = self.data.isotope_match
        self.temp_spectrum = self.data.temp_spectrum
        self.temp_drift = self.data.temp_drift
        self.changing = defaultdict(lambda: None)
        self.changing['status'] = False
        self.changing['peaks'] = False
        self.molecular_formula = None
        self.message = self.data.message
        self.noise = self.data.primary_data["ion_intensity"]
        self.all_rt = self.data.all_rt
        self.all_dt = self.data.all_dt
        self.mzmldata = self.data.mzmldata
        self.experimental_isotope = defaultdict(lambda: defaultdict(list))
        self.run_all()


    @timeit
    def formula_define(self, ions):
        ion_not = ["[M+H]+", "[M-H]-", "[M]-", "[M]+", "+electron", "-electron", "[M+2H]2+", "[M-2H]2-", "[M]2-", "[M]2+"]
        for y in ions:
            if not y in ion_not:
                if "+" in y:
                    x = y.strip("+")
                    ions = self.molecular_formula + x
                    ions = ions.strip()
                    df = isotopic_table(ions)
                    self.selected_ions[y] = df
                if "-" in y:
                    x = y.strip("-")
                    df2 = isotopic_table(x)
                    ion_substract = df2.mz.values[0]
                    ions = self.molecular_formula
                    df = isotopic_table(ions)
                    df["mz"] = df["mz"].map(lambda x: x-ion_substract)
                    self.selected_ions[y] = df
            elif y in ["[M+H]+"]:
                ions = self.molecular_formula + "H"
                ions = ions.strip()
                self.selected_ions[y] = isotopic_table(ions)
            elif y in ["[M]-", "[M]+", "+electron", "-electron"]:
                ions = self.molecular_formula
                self.selected_ions[y] = isotopic_table(ions)
            elif y in ["[M-H]-"]:
                ions = self.molecular_formula
                df = isotopic_table(ions)
                df["mz"] = df["mz"].map(lambda x: x - 1.0078)
                self.selected_ions[y] = df
            elif y in ["[M-2H]2-"]:
                ions = self.molecular_formula
                df = isotopic_table(ions)
                df["mz"] = df["mz"].map(lambda x: x - 1.0078 * 2).map(lambda x: x/2)
                self.selected_ions[y] = df
            elif y in ["[M+2H]2+"]:
                ions = self.molecular_formula
                df = isotopic_table(ions)
                df["mz"] = df["mz"].map(lambda x: x + 1.0078 * 2).map(lambda x: x/2)
                self.selected_ions[y] = df
            elif y in ["[M]2+", "[M]2-"]:
                ions = self.molecular_formula
                df = isotopic_table(ions)
                df["mz"] = df["mz"].map(lambda x: x/2)
                self.selected_ions[y] = df



    @timeit
    def run_all(self):
        if self.data.primary_data['use_data'] != "yes":
            self.mzmldata["mzml"] = None
            import gc
            gc.enable()
            gc.collect()
            spectra_generator = MS1OnlyMzML(self.data.primary_data["mzml"])
            ps = list(next(spectra_generator).keys())
            ds = (s for s in spectra_generator)
            df = pd.DataFrame.from_records(ds)
            tt = df.apply(lambda x: extract_decorator(x, ps), axis = 1)
            dtable = pl.concat(tt.values, how = "vertical_relaxed")
            self.mzmldata["mzml"] = dtable
        else:
            dtable = self.mzmldata["mzml"]

        for x, y in self.mass_df.items():
            for i , j in zip(y.values, list(y.index)):
                mass = i
                mass_plus = condition_max(mass, self.PPM1)
                mass_minus = condition_min(mass, self.PPM1)
                sorted_df = dtable.filter((pl.col("mz") < mass_plus) & (pl.col("mz") > mass_minus))
                index_size = sorted_df.select(pl.count()).collect()[0,0]
                theoretical_mass = [i] * index_size
                sorted_df = sorted_df.collect().with_columns(theoretical_mass = pl.lit(theoretical_mass))
                sorted_df = sorted_df.lazy()
                self.dframe[x][j] = sorted_df



    def mix_all(self, x, y):
        '''each molecular formula is studied for the presence of primary ion peak leads to status of changing attribute to true,
        the corresponding rt list is screened for isotope match and for peak picking for each ion queried by the user.
        followed by plotting the drift time and retention time spectrum.
        '''

        self.molecular_formula = x
        self.formula_define(list(y.index))
        calc_object = calculation.Calculation(self)
        calc_object.rt_window_formula()
        if self.changing['status']:
            rt_list = self.found_molecular_formula[self.molecular_formula]
            if rt_list.size > 0:
                calc_object.ion_dataframe()
                for rt in rt_list:
                    label = self.molecular_formula + str(rt)
                    ion_list = self.rt_found_ion[label]
                    for ion in ion_list:
                        plot_object = Xic_Eic(self)
                        plot_object.ticplot(ion, rt)
                        plot_object.plotrt_EIC()
                        plot_object.plotdt_IC()




    def dataframe_all(self):
        xx = [x for x in self.mass_df.keys()]
        yy = [y for y in self.mass_df.values()]
        df = pd.DataFrame({"MF": xx, "ion_list": yy})
        df.apply(lambda x: self.mix_all(x[0], x[1]), axis = 1)
        kk = []
        k = []
        if len(self.full_result.keys()) > 0:
            for x, y in self.full_result.items():
                dd = []
                td = []
                for c, t in y.items():
                    dd.append(t)
                    td.append(c)
                df1 = pd.concat(dd, keys = td)
                df1 = df1.reset_index()
                df1 = df1.drop(['level_1'], axis = 1)
                df1 = df1.rename({"level_0": "ion_type"}, axis = 1)
                kk.append(df1)
                k.append(x)
            df_final = pd.concat(kk, keys = k)
            df_final = df_final.reset_index()
            df_final = df_final.drop(['level_1'], axis = 1)
            df_final = df_final.rename({"level_0": "molecular_formula"}, axis = 1)
            return df_final
        else:
            self.message["Alert"].append(f'no molecular formula found')
            return None
