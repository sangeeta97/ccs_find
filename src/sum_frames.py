
import os
import pandas as pd
import statistics
from scipy import interpolate
import numpy as np
from collections import defaultdict
import operator
from collections import namedtuple
import plotly.express as px
from .peak import *
import tempfile
from functools import wraps
import time
import matplotlib.pyplot as plt
from scipy.signal import peak_widths
import numpy as np
import plotly.express as px
from . import iso_pat
import polars as pl



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



def valid_rt(target, query_df, test_case, PPM2, PPM1, ion):
    target1 = test_case['mz'].values[0]
    target2 = test_case['mz'].values[1]
    target3 = test_case['mz'].values[2]
    mz1 = []
    mz2 = []
    intensity1 = []
    intensity2 = []
    min_2 = condition_min(target2,  PPM2)
    max_2 = condition_max(target2, PPM2)
    test_case = test_case['intensity'].values[1]
    test_case = float(test_case)
    query_target2 = query_df.filter((pl.col("mz") >= min_2) & (pl.col("mz") <= max_2))
    index_range =  np.ravel(query_target2.select("index").collect().to_pandas().values)
    if index_range.size > 0:
        query_df = query_df.filter(pl.col("index").is_in(index_range))
        mz1_min = condition_min(target1, PPM1)
        mz1_max = condition_max(target1, PPM1)
        query_target1 =  query_df.filter((pl.col('mz') >= mz1_min) & (pl.col("mz") <= mz1_max))
        query_target1 = query_target1.collect().to_pandas()
        if query_target1.index.size > 0:
            intensity_thresold = np.percentile(query_target1['intensity'].values, 25)
            query_target1 = query_target1.loc[query_target1['intensity'] >= intensity_thresold, :]
            query_mz1 = query_target1["mz"].values
            mz1_mean = np.median(query_mz1)
            query_intensity1 = query_target1["intensity"].values
            sum_intensity1 = np.median(query_intensity1)
            mz2_min = condition_min(target2, PPM2)
            mz2_max = condition_max(target2, PPM2)
            query_target2 =  query_df.filter((pl.col('mz') >= mz2_min) & (pl.col("mz") <= mz2_max))
            query_target2 = query_target2.collect().to_pandas()
            # mz2_mask =  (query_df['mz'] >= mz2_min) & (query_df['mz'] <= mz2_max)
            # query_target2 = query_df.loc[mz2_mask, :]
            intensity_thresold = np.percentile(query_target2['intensity'].values, 25)
            query_target2 = query_target2.loc[query_target2['intensity'] >= intensity_thresold, :]
            query_mz2 = query_target2["mz"].values
            mz2_mean = np.median(query_mz2)
            query_intensity2 = query_target2["intensity"].values
            sum_intensity2 = np.median(query_intensity2)
            intensity_ratio = np.divide(sum_intensity2, sum_intensity1)
            intensity_ratio = np.multiply(intensity_ratio, 100)
            kd = np.absolute(intensity_ratio - test_case)
            kd = np.ravel(kd)
            dk = kd < test_case
            dk = np.ravel(dk)
            ratio_bool = dk[0]
            return ratio_bool, mz1_mean, sum_intensity1, mz2_mean, sum_intensity2, query_df

        else:
            return False, 0, 0, 0, 0, 0

    else:
        return False, 0, 0, 0, 0, 0


def spectrum_confirm(molecular_formula, pp, target, peak_dataframe, test_case, PPM1, PPM2, noise, ion, peakwidth, temp):
    bool_rt = {}
    peakwidth = float(peakwidth)
    pp = np.array(pp)
    pp = np.ravel(pp)
    for p in pp:
        p_min = float(p) - peakwidth
        p_max = float(p) + peakwidth
        query_df = peak_dataframe.filter((pl.col("rt") >= p_min) & (pl.col("rt") <= p_max))
        ratio_bool, mz1, intensity1, mz2, intensity2, query_df = valid_rt(target, query_df, test_case, PPM2, PPM1, ion)
        if ratio_bool:
            gg = iso_pat.Isotope_Plot(query_df, test_case, ion, p, molecular_formula, temp, target)
            gg.process_all()
            mz = gg.mzlist
            bool_rt[p] = "yes"
        else:
            bool_rt[p] = "no"
            mz = np.array([])
    return bool_rt, mz
