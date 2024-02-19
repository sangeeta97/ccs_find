import os
import json
import pandas as pd
import statistics
from scipy import interpolate
import numpy as np
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
import numpy as np
import plotly.express as px
import polars as pl
from scipy.signal import savgol_filter


def filter(xx):
    xx = np.array(xx)
    return savgol_filter(xx, 5, 2, mode='nearest')




def timeit(func):
    @wraps(func)
    def timeit_wrapper(*args, **kwargs):
        start_time = time.perf_counter()
        result = func(*args, **kwargs)
        end_time = time.perf_counter()
        total_time = end_time - start_time
        print(f'Function {func.__name__}{args} {kwargs} Took {total_time:.4f} seconds')
        return result
    return timeit_wrapper





def timeit1(func):
    @wraps(func)
    def timeit_wrapper(*args, **kwargs):
        start_time = time.perf_counter()
        result = func(*args, **kwargs)
        end_time = time.perf_counter()
        total_time = end_time - start_time
        path = os.getcwd()
        path = os.path.join(path, "log.txt")
        f = open(path, "w")
        f.write(f'the total analysis time was {total_time:.4f} seconds')
        f.close()
        print(f'Function {func.__name__}{args} {kwargs} Took {total_time:.4f} seconds')
        return result
    return timeit_wrapper



def find_fwhm(half_max, start, end, xnew, ynew):
    try:

        fwhm_all = []
        for i, j, k in zip(start, end, half_max):
            peak_range = np.arange(i, j, 1)
            peak_range_y = ynew[peak_range]
            peak_range_x = xnew[peak_range]
            peak_range_y = np.ravel(peak_range_y)
            bool_array = peak_range_y > k
            peak_range_x = peak_range_x[bool_array]
            first = peak_range_x.min()
            last = peak_range_x.max()
            fwhm = last - first
            fwhm_all.append(fwhm)
        fwhm_all = np.array(fwhm_all)
        fwhm_all = np.ravel(fwhm_all)
        return fwhm_all

    except:
        return 1


@timeit
def detect_rt(mm, thresold):
    '''mm is polars LazyFrame'''
    confirm_p = np.ravel(mm.select("index").collect().to_pandas().values)
    pp = mm.select(['rt', 'intensity']).collect().to_pandas()
    pp= pp[['rt', 'intensity']].groupby('rt')['intensity'].apply(max)
    rt= np.array(pp.index.values)
    intensity= pp.values
    intensity = filter(intensity)
    noise = estimate_noise(intensity)
    baseline = estimate_baseline(intensity, noise)
    start, peaks, end = detect_peaks(intensity, noise, baseline, thresold)
    peak_intensities  = np.array([intensity[x] for x in peaks])
    if peaks.size > 0:
        peak = True
        max_intensity = max(peak_intensities)
        peaks_index = peak_intensities > max_intensity/3
        peaks = peaks[peaks_index]
        start = start[peaks_index]
        end = end[peaks_index]
        rt_peaks = rt[peaks]
        return rt[peaks], peak
    else:
        peak = False
        return None, peak


def extract_decorator22(spectrum, noise):
    thresold= 2.0 if noise == "yes" else 20.0
    mz = spectrum.mz[spectrum.i > thresold]
    intensity = spectrum.i[spectrum.i > thresold]
    whole = namedtuple('Whole', ['mz', 'intensity', 'dt', 'rt'])
    gg = round(spectrum.get("MS:1002476"), 2)
    if mz.size > 0:
        li = [mz, intensity, gg, round(spectrum.scan_time_in_minutes(), 3)]
        mm = whole._make(li)
        return mm



def extract_decorator(spec, noise):
    pp = list(spec.keys())
    mz = spec[pp[-2]]
    intensity = spec[pp[-1]]
    index = spec['index']
    thresold = 2 if noise == "yes" else 20
    mz = mz[intensity > thresold]
    intensity = intensity[intensity > thresold]
    rt = spec['scanList']['scan'][0]['scan start time']
    rt = round(rt, 3)
    dd = spec['scanList']['scan'][0]['ion mobility drift time']
    dd = round(dd, 3)
    size = mz.size
    if mz.size > 1:
        all = {"mz": mz, "rt": [rt] * size, "drift_time": [dd] * size, "index": [index] * size, "intensity": intensity}
        df = pl.LazyFrame(all, schema = {"mz": pl.Float32, "rt": pl.Float32, "drift_time": pl.Float32, "index": pl.Float32, "intensity": pl.Float32})
        return df
    return pl.LazyFrame({"mz": None, "rt": None, "drift_time": None, "index": None, "intensity": None})
