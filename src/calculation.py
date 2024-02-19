import os
import json
import statistics
from scipy import interpolate
import numpy as np
from numpy.linalg import norm
import copy
import operator
import math
from collections import defaultdict
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
from . import sum_frames
from molmass import Formula
import polars as pl
from scipy.signal import savgol_filter


def filter(xx):
    xx = np.array(xx)
    return savgol_filter(xx, 5, 2, mode='nearest')





'''
ion_dict : dictionary of ion and the corresponding monoisotopic mass value for the ion
selected_ions : dictionary of ion and the corresponding MF with adduct to isotopic confirmation
found_molecular_formula : if the molecular formula is found based on the accurate mass of primary ions, keys are the molecular formula and values are the retention time for those MF.
dt_present : dataframe of combination of rt and ion_type with selected peak in drift time space.
rt_found_ion: key is rt and value is list of ion types found at that rt
ions_data[ion][i]: first key is ion_type and second key is rt and value is a raw dataframe.
dt_peaks[ion][i]: first key is ion_type and second key is rt and value is picked peak.
'''


class Calculation:
    '''
    molecular formula parsing and giving table of drift time and rt information

    '''
    def __init__(self, data):
        self.data = data
        self.monoPPM = self.data.monoPPM
        self.PPM1 = self.data.PPM1
        self.PPM2 = self.data.PPM2
        self.mass_df = self.data.mass_df
        self.primary_ion = self.data.primary_ion
        self.full_result = self.data.full_result
        self.dframe = self.data.dframe
        self.ion_dict = self.data.ion_dict
        self.selected_ions = self.data.selected_ions
        self.rt_found_ion = self.data.rt_found_ion
        self.ions_data = self.data.ions_data
        self.found_molecular_formula = self.data.found_molecular_formula
        self.dt_peaks = self.data.dt_peaks
        self.isotope = self.data.isotope
        self.isotope_match = self.data.isotope_match
        self.temp_spectrum = self.data.temp_spectrum
        self.temp_drift = self.data.temp_drift
        self.peaks = False
        self.molecular_formula = self.data.molecular_formula
        self.message = self.data.message
        self.noise = self.data.noise
        self.all_rt = self.data.all_rt
        self.all_dt = self.data.all_dt
        self.mzmldata = self.data.mzmldata
        self.peakwidth = self.data.peakwidth
        # self.experimental_isotope = self.data.experimental_isotope



    @timeit
    def rt_window_formula(self):
        '''The primary ion of the current molecular formula is being tested for
        rt match and isotopic pattern, it confirms about the rt of the molecular formula.
        that rt value further used to explore the possibility of presence of adduct ions
        which can be confirmed by the presence of isotopic pattern in the next class method ion_dataframe

        dframe leads to LayzFrame
        '''
        self.data.changing['status'] = False
        mm = self.dframe[self.molecular_formula][self.primary_ion]
        index_size = mm.select(pl.count()).collect()[0,0]
        if index_size > 5:
            pp, peak = detect_rt(mm, self.noise)
            if peak:
                target = self.ion_dict[self.molecular_formula]
                test_case = self.selected_ions[self.data.primary_ion]
                bool_rt, _ = sum_frames.spectrum_confirm(self.molecular_formula, pp, target, self.mzmldata["mzml"], test_case, self.PPM1, self.PPM2, self.noise, None, self.peakwidth, self.temp_spectrum)
                final_rt = [x for x in bool_rt.keys() if bool_rt[x] != "no"]
                final_isotope = [bool_rt[x] for x in final_rt]
                if len(final_rt) > 0:
                    # for i, y in zip(final_rt, final_isotope):
                    #     self.experimental_isotope[self.molecular_formula][i].append(y)
                    final_rt = np.array(final_rt)
                    if final_rt.size > 0:
                        self.data.changing['status'] = True
                        pp = final_rt
                        self.found_molecular_formula[self.molecular_formula] = pp
                    else:
                        self.message["warning"].append(f'{self.molecular_formula} not found')
                        self.found_molecular_formula[self.molecular_formula] = np.array([])
                        return np.array([])
                else:
                    self.message["warning"].append(f'{self.molecular_formula} not found')
                    return np.array([])
            else:
                self.message["warning"].append(f'{self.molecular_formula} not found')
                return np.array([])
        else:
            self.message["warning"].append(f'{self.molecular_formula} not found')
            return np.array([])



    def ion_dataframe(self):
        '''The rt values were used to scan for presence of the adduct ions with isotopic pattern match
        including for primary ion, the dataframe with drift time value was used to calculate resolving resolving_power
        of peak, and full information regarding rt, intensity, mz value for each dt peak is summarized for each ion.
        '''
        mass_all = self.dframe[self.molecular_formula]
        ion_dict = self.mass_df[self.molecular_formula]
        ion_dict = {i:j for i, j in zip(list(ion_dict.index), ion_dict.values)}
        pp = self.found_molecular_formula[self.molecular_formula]
        print("printing pp")
        print(pp)
        self.rt = pp
        rt_range_min = pp - float(self.peakwidth)

        rt_range_max = pp + float(self.peakwidth)
        for ion in mass_all.keys():
            mm = mass_all[ion]
            lc = []
            for i in range(len(pp)):
                ll= mm.filter((pl.col('rt') >= rt_range_min[i]) & (pl.col('rt') <= rt_range_max[i]))
                ll = ll.collect().to_pandas()
                if not ll.index.size > 5:
                    continue
                lt = ll[['index' ,'drift_time', 'intensity', 'mz', 'rt', 'theoretical_mass']].sort_values(by = ['intensity'], ascending = False).drop_duplicates(subset = ['drift_time'], keep = "first")
                spec_index = lt['index'].values[0]
                target = ion_dict[ion]
                test_case = self.selected_ions[ion]
                bool_rt, mz_tuple = sum_frames.spectrum_confirm(self.molecular_formula, pp[i], target, self.mzmldata["mzml"], test_case, self.PPM1, self.PPM2, self.noise, ion, self.peakwidth, self.temp_spectrum)
                final_rt = [x for x in bool_rt if bool_rt[x] != "no"]
                final_rt = np.array(final_rt)
                if final_rt.size > 0:

                    ''' Stored rt value of primary ion used as reference (fixed value) the dataframe for the plotting for all the adduction ions as well, as adduct ions might have other rt values as well'''
                    rtE = pp[i]
                    mz = mz_tuple[0]
                    rt = mz_tuple[1]
                    intensity = mz_tuple[2]
                    noise = estimate_noise(intensity)
                    baseline = estimate_baseline(intensity, noise)
                    start, peaks, end = detect_peaks(intensity, noise, baseline, self.noise)
                    peak_intensities  = np.array([intensity[x] for x in peaks])
                    max_intensity = max(peak_intensities)
                    peaks_index = peak_intensities > max_intensity/3
                    peaks = peaks[peaks_index]
                    start = start[peaks_index]
                    end = end[peaks_index]
                    rt = rt[peaks]
                    label = self.molecular_formula + str(rtE)
                    self.ions_data[label][ion] = lt
                    df = self.peak_dt(lt, ion, rt, mz, rtE)
                    if self.peaks:
                        label = self.molecular_formula + str(rtE)
                        self.rt_found_ion[label].append(ion)
                        lc.append(df)

            if len(lc) > 0:
                dd = pd.DataFrame()
                lc.append(dd)
                dt_peaks = pd.concat(lc)
                self.full_result[self.molecular_formula][ion] = dt_peaks




    def peak_dt(self, df, z, i, mz1, rtE):
        '''peak picking using dt dimesion and summarizing information about each ion in the table
        '''
        ll = df.sort_values(by=['drift_time'])
        xnew = ll['drift_time'].values
        ynew = ll['intensity'].values
        ynew = filter(ynew)
        mz = ll['mz'].values
        spec = ll['index'].values
        theoretical_mass = ll['theoretical_mass'].values
        noise = estimate_noise(ynew)
        baseline = estimate_baseline(ynew, noise)
        start, peaks, end = detect_peaks(ynew, noise, baseline, self.noise)
        self.peaks = False
        if peaks.size > 0:
            self.peaks = True
            peak_intensities  = np.array([ynew[x] for x in peaks])
            max_intensity = max(peak_intensities)
            peaks_index = peak_intensities > max_intensity/5
            peaks = peaks[peaks_index]
            start = start[peaks_index]
            end = end[peaks_index]
            max_index = np.argmax(ynew)
            if mz1.size > 0:
                mz_mid = np.median(mz1)
            else:
                mz_mid = np.median(mz)
            end = end - 1
            peak_mid = xnew[peaks]
            peak_mid = peak_mid.round(3)
            spec_number = spec[peaks]
            rt_mid = np.resize(i, (len(peaks), 1))
            rt_mid = rt_mid.round(3)
            rt_mid = rt_mid.flatten()
            intensity_mid = ynew[peaks]
            intensity_mid = intensity_mid.round(3)
            half_max = intensity_mid * 0.5
            theoretical_mass = theoretical_mass[peaks]
            peak_start = xnew[start]
            peak_end = xnew[end]
            fwhm = find_fwhm(half_max, start, end, xnew, ynew)
            dict_all = {"mz_top": mz_mid, "peak_start": peak_start, "peak_end": peak_end, "dt_mid": peak_mid, "fwhm": fwhm, "spec_number": spec_number, "rt_mid": rt_mid, "theoretical_mass": theoretical_mass}
            df = pd.DataFrame(data = dict_all)
            self.dt_peaks[z][i[0]] = df
            df['number of conformers'] = df['dt_mid'].size
            df['Error(PPM)'] = np.subtract(df['mz_top'].values, df['theoretical_mass'].values)
            df['Error(PPM)'] = np.divide(df['Error(PPM)'].values, df['theoretical_mass'].values)
            df['Error(PPM)'] = df['Error(PPM)'] * 1000000
            fwhm = list(fwhm)
            fwhm = [x if x > 0 else 1 for x in fwhm]
            fwhm = np.array(fwhm)
            df['resolving_power'] = np.divide(peak_mid, fwhm)
            df["rt"] = rtE
            df['resolving_power'] = df['resolving_power'].map(lambda x: round(x, 2))
            final_result = df[["mz_top", "Error(PPM)", "number of conformers", "dt_mid", "fwhm", "rt_mid", "resolving_power", "rt"]]
            final_result.columns = ["mz_measured", "Error(PPM)", "#conformer", "drift_time", "fwhm", "retention_time", "resolving_power", "rt"]
            return final_result
        else:
            self.peaks = False
            self.data.changing['peaks'] = False
            return None
