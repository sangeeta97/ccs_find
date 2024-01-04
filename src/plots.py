import os
import json
import statistics
from scipy import interpolate
import numpy as np
from numpy.linalg import norm
import copy
from collections import defaultdict
import operator
from collections import namedtuple
import plotly.express as px
import tempfile
from functools import wraps
import time
import pandas as pd
from bokeh.io import output_notebook
from bokeh.plotting import figure, show, output_file
from bokeh.models import ColumnDataSource, HoverTool, NumeralTickFormatter, Label
from bokeh.palettes import Category10
from molmass import Formula
from bokeh.plotting import figure, output_file, show
from bokeh.resources import CDN
from bokeh.embed import file_html
from scipy.signal import savgol_filter
import polars as pl



def filter(xx):
    xx = np.array(xx)
    return savgol_filter(xx, 5, 2, mode='nearest')



def add_axis_labels(p):
    p.xaxis.axis_label = 'Fragment m/z'
    p.xaxis.axis_label_text_font_size = '10pt'
    p.xaxis.major_label_text_font_size = '9pt'

    p.yaxis.axis_label = 'Intensity'
    p.yaxis.axis_label_text_font_size = '10pt'
    p.yaxis.major_label_text_font_size = '9pt'
    p.yaxis.formatter = NumeralTickFormatter(format='0.')


def create_p(width=800, height=300,
            main_title='plot'):
    tooltips = [
        ('m/z','@mz{0.0000}'),
        ('Int','@intensity')
        ]
    p = figure(
        plot_width=width, plot_height=height,
        title = main_title,
        tools = 'xwheel_zoom,xpan,box_zoom,undo,reset',
        tooltips=tooltips
        )
    return p


class Xic_Eic:

    def __init__(self, data):
        self.data = data
        self.experimental_isotope = self.data.experimental_isotope
        self.experimental_isotope = self.data.experimental_isotope
        dtable = self.data.mzmldata["mzml"]
        max_rt = dtable.select("rt").collect().max().to_pandas().values[0][0]
        max_drift_time = dtable.select("drift_time").collect().max().to_pandas().values[0][0]
        rt_range = np.arange(0, max_rt, 0.01)
        dt_range = np.arange(0, max_drift_time, 0.05)
        self.plot_rt = pd.DataFrame({"rt": rt_range})
        self.plot_dt = pd.DataFrame({"drift_time": dt_range})


    def plotrt_EIC(self):
        try:

            df_all = []
            rt_list = self.data.found_molecular_formula[self.data.molecular_formula]
            for rt in rt_list:
                label = self.data.molecular_formula + str(rt)
                if label in list(self.data.rt_found_ion.keys()):
                    ion_list = self.data.rt_found_ion[label]
                    for ion in ion_list:
                        mm = self.data.dframe[self.data.molecular_formula][ion]
                        mm = mm.select(['rt', 'intensity']).collect().to_pandas()
                        mm = mm.sort_values(by=['rt'])
                        pp= mm[['rt', 'intensity']].groupby('rt')['intensity'].apply(max)
                        rt= np.array(pp.index.values)
                        intensity= pp.values
                        f = interpolate.interp1d(rt, intensity, kind = "linear", fill_value = "extrapolate")
                        xnew = np.arange(rt.min(), rt.max(), 0.01)
                        ynew = f(xnew)
                        ynew = [max(num, 0) for num in ynew]
                        data = {"rt": xnew, "intensity": ynew}
                        df1 = pd.DataFrame(data = data)
                        df1['ion_type'] = ion
                        min_val = df1["rt"].min()
                        max_val = df1["rt"].max()
                        test_df = self.plot_rt.loc[((self.plot_rt.rt < min_val) | (self.plot_rt.rt > max_val)), :]

                        df22 = df1.merge(test_df, on = "rt", how = "outer")
                        df22["ion_type"] = df22["ion_type"].fillna(ion)
                        df22["intensity"] = df22["intensity"].fillna(0.0)
                        df22 = df22.sort_values(by = ["intensity"], ascending = False)
                        df22 = df22.drop_duplicates(subset = ["rt"])

                        df22["intensity"]= filter(df22["intensity"].values)
                        df_all.append(df22)
            df = pd.concat(df_all)
            df = df.sort_values(by = ["rt"])
            x = df['rt'].values
            y = df['intensity'].values
            f = interpolate.interp1d(x, y, kind = "linear", fill_value = "extrapolate")
            ynew = f(x)
            ynew = [max(num, 0) for num in ynew]
            df['intensity'] = ynew
            fig = px.line(df, x="rt", y="intensity", color='ion_type', title=f'rt EIC plot of {self.data.molecular_formula}')
            fig.write_html(os.path.join(self.data.temp_spectrum, f"{self.data.molecular_formula}_rt_overlay.html"))
            self.data.all_rt[self.data.molecular_formula] = df

        except Exception as e:
            return e




    def plot3d_EIC(self):
        try:

            rt_list = self.data.found_molecular_formula[self.data.molecular_formula]
            for rt in rt_list:
                label = self.data.molecular_formula + str(rt)
                if label in list(self.data.rt_found_ion.keys()):
                    ion_list = self.data.rt_found_ion[label]
                    for ion in ion_list:
                        mm = self.data.dframe[self.data.molecular_formula][ion]
                        mm = mm.collect().to_pandas()
                        fig = px.scatter_3d(mm, x='rt', y='drift_time', z='mz', color='intensity')
                        fig.write_html(os.path.join(self.data.temp_spectrum, f"{self.data.molecular_formula}_{ion}_plot3d.html"))

        except Exception as e:
            return e





    def plotdt_IC(self):
        '''molecular_formula+rt label as key value stores the ion list in the dict rt_found_ion
        ions_data dict stores the drift_time dataframe for a label + ion.
        '''
        try:
            df_all = []
            rt_list = self.data.found_molecular_formula[self.data.molecular_formula]
            for rt in rt_list:
                label = self.data.molecular_formula + str(rt)
                if label in list(self.data.rt_found_ion.keys()):
                    ion_list = self.data.rt_found_ion[label]
                    for ion in ion_list:
                        label = self.data.molecular_formula + str(rt)
                        df = self.data.ions_data[label][ion]
                        min_val = df["drift_time"].min()
                        max_val = df["drift_time"].max()
                        test_df = self.plot_dt.loc[((self.plot_dt.drift_time < min_val) | (self.plot_dt.drift_time > max_val)), :]
                        df22 = df.merge(test_df, on = "drift_time", how = "outer")
                        df22["intensity"] = df22["intensity"].fillna(0.0)
                        x= df22['drift_time'].values
                        y= df22['intensity'].values
                        y = savgol_filter(y, 5, 2, mode='mirror')
                        f = interpolate.interp1d(x, y, kind = "next",  fill_value="extrapolate", bounds_error = False)
                        xnew = np.arange(df22['drift_time'].min(), df22['drift_time'].max(), 0.05)
                        ynew = f(xnew)
                        ynew = np.where(ynew < 0, 0, ynew)
                        ynew = list(ynew)
                        ynew = [max(num, 0) for num in ynew]
                        ynew = savgol_filter(ynew, 5, 2, mode='mirror')
                        ynew = [max(num, 0) for num in ynew]
                        data = {"drift_time": xnew, "intensity": ynew}
                        df = pd.DataFrame(data)
                        df["molecular_formula"] = self.data.molecular_formula
                        df["mf_rt"] = self.data.molecular_formula + "_" + str(rt)
                        df["ion_type"] = ion
                        df["rt"] = rt
                        df_all.append(df)

            df22 = pd.concat(df_all)
            df22 = df22.sort_values(by = ["drift_time"])
            x= df22['drift_time'].values
            y= df22['intensity'].values
            b = str(rt)
            fig = px.line(df22, x="drift_time", y="intensity", color='ion_type', title=f'overlay_IM plot of {self.data.molecular_formula}_{b}')
            fig.write_html(os.path.join(self.data.temp_drift, f"{self.data.molecular_formula}_{b}_IMoverlay.html"))
            self.data.all_dt[self.data.molecular_formula] = df22
        except Exception as e:
            return e



    def ticplot(self, ion, rt):
        try:
            label = self.data.molecular_formula + str(rt)
            ion_list = self.data.rt_found_ion[label]
            if ion in ion_list:
                df = self.data.ions_data[label][ion]
                min_val = df.drift_time.min()
                max_val = df.drift_time.max()
                test_df = self.plot_dt.loc[((self.plot_dt.drift_time < min_val) | (self.plot_dt.drift_time > max_val)), :]
                df = df[["drift_time", "intensity"]]
                df = df.merge(test_df, on = "drift_time", how = "outer")
                df["intensity"] = df["intensity"].fillna(0.0)
                df = df.drop_duplicates(subset = ["drift_time"])
                df = df.sort_values(by = ["drift_time"])
                x= df['drift_time'].values
                y= df['intensity'].values
                y = savgol_filter(y, 5, 2, mode='mirror')
                f = interpolate.interp1d(x, y, kind = "next", fill_value = "extrapolate", bounds_error = False)
                xnew = np.arange(df['drift_time'].min(), df['drift_time'].max(), 0.05)
                ynew = f(xnew)
                ynew = np.where(ynew < 0, 0, ynew)
                ynew = [max(num, 0) for num in ynew]
                ynew = savgol_filter(ynew, 5, 2, mode='mirror')
                ynew = list(ynew)
                ynew = [max(num, 0) for num in ynew]
                rt = float(rt)
                rt = round(rt, 3)
                rt = str(rt)
                fig = px.line(x=xnew, y=ynew, labels={'x':'drift_time', 'y':'intensity'}, title=f'IM plot of {self.data.molecular_formula}_{rt}_{ion}')
                fig.write_html(os.path.join(self.data.temp_drift, f"{self.data.molecular_formula}_{ion}_{rt}_IM.html"))

        except Exception as e:
            print(e)
