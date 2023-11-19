import os
import pandas as pd
import math
import numpy as np

mg = 28.006148
#
# tfix = 0.428379630615936
#
# beta = 0.13581225201424416


# DriftGas mass="28.006148">N2</DriftGas>
# <TFix>0.428379630615936</TFix>
# <Beta>0.13581225201424416<


# Mass of the drift gas (mg) required. Default is N2 = 28.01 g/mol
# [(tA – tfix)*(z/beta*gamma) = DTCCSN2
# tA = arrival time of the ion of interest (“drift time”)
# mi = mass of the ion (note, this is not equal to m/z for z ≥ 2)
# gamma = sqrt[mi/(mg+mi)]
# z = integer charge state


# [(tA – tfix)*(z/beta*gamma) = DTCCSN2


#primary_ion
#drift_gas
#formula
#mzml
#calibration
#beta
#tfix
#buffer_text

def calculate(row,  tfix, beta, mg = mg):
    mi = row['mz_measured']
    tA = row['drift_time']
    m = mg + mi
    inter_mi = np.divide(mi, m)
    gamma = math.sqrt(inter_mi)
    diff = tA - tfix
    b = beta * gamma
    b = np.divide(1, b)
    ccs = diff * b
    return ccs


class CCS:
    def __init__(self, data, drift_result):
        self.data = data
        self.drift_result = drift_result
        self.tfix = float(self.data.primary_data["tfix"])
        self.beta = float(self.data.primary_data["beta"])
        self.mg = 28.006148



    def tfix_beta(self):
        if float(self.data.primary_data["tfix"]):
            self.tfix = float(self.data.primary_data["tfix"])
        if float(self.data.primary_data["beta"]):
            self.beta = float(self.data.primary_data["beta"])
        jj = self.data.optional_data['calibration']
        if isinstance(jj, dict):
            self.tfix = float(self.data.optional_data["TFix"])
            self.beta = float(self.data.optional_data["Beta"])



    def calculate(self, row):
        mi = row['mz_measured']
        tA = row['drift_time']
        m = self.mg + mi
        inter_mi = np.divide(mi, m)
        gamma = math.sqrt(inter_mi)
        diff = tA - self.tfix
        b = self.beta * gamma
        b = np.divide(1, b)
        ccs = diff * b
        return ccs



    def finish_analysis(self):
        df1 = self.drift_result[['mz_measured', 'drift_time']]
        self.drift_result['ccs'] = df1.apply(self.calculate, axis = 1)
        self.drift_result = self.drift_result.drop(['fwhm'], axis = 1)
        self.drift_result['mz_measured'] = self.drift_result['mz_measured'].map(lambda x: round(x, 4))
        self.drift_result['Error(PPM)'] = self.drift_result['Error(PPM)'].map(lambda x: round(x, 1))
        self.drift_result['ccs'] = self.drift_result['ccs'].map(lambda x: round(x, 2))
        return self.drift_result
