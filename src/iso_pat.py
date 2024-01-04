import os
import numpy as np
from collections import defaultdict
import operator
from collections import namedtuple
import pandas as pd
from scipy import interpolate
import plotly.express as px
from bokeh.layouts import column
from bokeh.plotting import figure, output_file, show
from bokeh.resources import CDN
from bokeh.embed import file_html
from bokeh.io import output_notebook
from bokeh.plotting import figure, show, output_file
from bokeh.models import ColumnDataSource, HoverTool, NumeralTickFormatter, Label
from bokeh.palettes import Category10
from bokeh.layouts import column
import polars as pl





def create_p(width=800, height=300,
            main_title='isotopic pattern XIC'):
    tooltips = [
        ('m/z','@mz{0.0000}'),
        ('Int','@intensity'),
        ("RT", '@rt')

        ]
    p = figure(
        width=width, height=height,
        title = main_title,
        tools = 'xwheel_zoom,xpan,box_zoom,undo,reset',
        tooltips=tooltips
        )
    return p

#Format axis labels
def add_axis_labels(p):
    p.xaxis.axis_label = 'Retention time (min)'
    p.xaxis.axis_label_text_font_size = '10pt'
    p.xaxis.major_label_text_font_size = '9pt'

    p.yaxis.axis_label = 'Intensity'
    p.yaxis.axis_label_text_font_size = '10pt'
    p.yaxis.major_label_text_font_size = '9pt'
    p.yaxis.formatter = NumeralTickFormatter(format='0.')

'''query_df is lazyframe'''


class Isotope_Plot(object):
    def __init__(self, query_df, test_case, ion, rt, mf, temp, target):
        self.query_df = query_df.collect()
        self.test_case = test_case
        self.ion = ion
        self.rt = rt
        self.mf = mf
        self.temp = temp
        self.target = target
        self.mzlist = None
        self.rtlist = None
        self.intensitylist = None


    def process_all(self):
        try:
            df = self.find_df()
            ll = self.group_df(df)
            self.plotting_bokeh(ll)
        except Exception as e:
            print(e)


    def find_df(self):
        '''raw data within a rt range were grouped based on the spectrum, each spectrum were looked for all three isotopic peaks
        '''
        self.theoretical_ratio = self.test_case['intensity'].astype(float).values[0:3]
        self.theoretical_masses = self.test_case['mz'].astype(float).values[0:3]
        sorted_df = self.query_df
        grouped = sorted_df.groupby(by = "index")
        list_pat1 = []
        list_pat2 = []
        list_pat3 = []
        whole1 = namedtuple('Whole1', ['mz', 'intensity', 'rt'])
        whole2 = namedtuple('Whole2', ['mz', 'intensity', 'rt'])
        whole3 = namedtuple('Whole3', ['mz', 'intensity', 'rt'])
        for i, g in grouped:
            g = g.to_pandas()
            g = g.sort_values(by = ["mz"])
            mz_exp = g.mz.values
            min = mz_exp.min()
            max = mz_exp.max()
            int_exp = g.intensity.values
            rt_exp = g.rt.values
            mz_index= np.searchsorted(mz_exp, self.theoretical_masses)
            mz_index[mz_index == mz_exp.size] = 0
            mz_index[mz_index == -1] = 0
            test = np.subtract(mz_index, 1)
            num1 = np.absolute(np.subtract(mz_exp[mz_index], self.theoretical_masses))
            num2 = np.absolute(np.subtract(mz_exp[test], self.theoretical_masses))
            mz_index = mz_index if num1[0] < num2[0] else test
            mz_index[mz_index == -1] = 0
            mz_exp = mz_exp[mz_index]
            int_exp = int_exp[mz_index]
            rt_exp = rt_exp[mz_index]
            int_exp = np.array([float(x) for x in int_exp])
            mz_exp = np.array([float(x) for x in mz_exp])
            rt_exp = np.array([float(x) for x in rt_exp])
            test_ratio1 = self.theoretical_ratio[0]/self.theoretical_ratio[1]
            exp_ratio1 = int_exp[0]/int_exp[1]
            diff1 = abs(test_ratio1 - exp_ratio1)
            test_ratio2 = self.theoretical_ratio[0]/self.theoretical_ratio[2]
            exp_ratio2 = int_exp[0]/int_exp[2]
            diff2 = abs(test_ratio2 - exp_ratio2)
            mz_diff1 = abs(mz_exp[1] - self.theoretical_masses[1])
            mz_diff2 = abs(mz_exp[2] - self.theoretical_masses[2])
            if diff1 < test_ratio1 and diff2 < test_ratio2:
                li1 = [mz_exp[0], int_exp[0], rt_exp[0]]
                list_pat1.append(whole1._make(li1))
                li2 = [mz_exp[1], int_exp[1], rt_exp[1]]
                list_pat2.append(whole2._make(li2))
                li3 = [mz_exp[2], int_exp[2], rt_exp[2]]
                list_pat3.append(whole3._make(li3))

        df1 = pd.DataFrame(list_pat1)
        df1["mask"] = df1["mz"].map(lambda x: x - self.theoretical_masses[0]).map(lambda x: abs(x))
        df1 = df1[df1["mask"] < 0.01]
        df2 = pd.DataFrame(list_pat2)
        df2["mask"] = df2["mz"].map(lambda x: x - self.theoretical_masses[1]).map(lambda x: abs(x))
        df2 = df2[df2["mask"] < 0.01]
        df3 = pd.DataFrame(list_pat3)
        if df3.index.size > 0:
            df3["mask"] = df3["mz"].map(lambda x: x - self.theoretical_masses[2]).map(lambda x: abs(x))
            df3 = df3[df3["mask"] < 0.01]
            dfA2 = pd.concat([df1, df2, df3], keys = ["mono", "di", "tri"])
            dfA2 = dfA2.reset_index()
            return dfA2
            # ll = self.group_df(dfA2)
            # plotting_final(ll, self.theoretical_masses, self.theoretical_ratio)

    def transform(self, df):
        df = df.sort_values(by = ["rt"])
        mz = df["mz"].values
        value = df["intensity"].max()
        min_rt = df["rt"].min()
        min_rt_range = np.arange(min_rt - 3.5, min_rt, 0.1)
        max_rt = df["rt"].max()
        max_rt_range = np.arange(max_rt, max_rt + 3.5, 0.1)
        full_range = np.concatenate([min_rt_range, max_rt_range])
        df11 = pd.DataFrame({"rt": full_range, "intensity1": 0.0})
        df1 = df.merge(df11, on = "rt", how = "outer")
        df1 = df1.fillna(0.0)
        df1["intensity"] = np.add(df1["intensity"].values, df1["intensity1"].values)
        df1 = df1.sort_values(by = ["intensity"], ascending = False)
        df1 = df1.sort_values(by = ["rt"])
        df1 = df1.fillna(0.0)
        x= df1['rt'].values
        y= df1['intensity'].values
        f = interpolate.interp1d(x, y, kind = "linear",  fill_value="extrapolate", bounds_error = False)
        xnew = np.arange(df1['rt'].min(), df1['rt'].max(), 0.1)
        ynew = f(xnew)
        f = interpolate.interp1d(xnew, ynew, kind = "linear",  fill_value="extrapolate", bounds_error = False)
        ynew = f(xnew)
        df = pd.DataFrame({"rt": xnew, "intensity": ynew})
        diff = value - df.intensity.max()
        df["intensity"] = df["intensity"].map(lambda x: np.add(x, diff))
        mz = np.resize(mz, ynew.size)
        df["mz"] = mz
        return df


    def group_df(self, df):
        from bokeh.plotting import figure, output_file, show
        from bokeh.resources import CDN
        from bokeh.embed import file_html
        max_intensity = df["intensity"].max()
        df["intensity"] = df["intensity"].map(lambda x: np.divide(x, max_intensity)).map(lambda x: np.multiply(x, 100))
        df = df.sort_values(by = ["intensity"], ascending = False)
        df = df[["level_0", "mz", "intensity", "rt", "mask"]]
        df.columns = ["label", "mz", "intensity", "rt", "mask"]
        df = df.sort_values(by = ["mz"])
        grouped = df.groupby("label")

        ll = []
        for i, g in grouped:
            ll.append(self.transform(g))
        return ll


    def plotting_bokeh(self, ll):
        if len(ll) == 3:
            a, b, c = ll
            self.mzlist = b.mz.values
            self.rtlist = b.rt.values
            self.intensitylist = b.intensity.values
            t1, t2, t3 = list(self.theoretical_ratio)
            m1, m2, m3 = list(self.theoretical_masses)
            df = ll[0]
            c1, c2, c3 = ColumnDataSource(pd.DataFrame({"rt": np.median(df.rt.values), "intensity": t1, "mz": m1}, index = [0])), ColumnDataSource(pd.DataFrame({"rt": np.median(df.rt.values), "intensity": t2, "mz": m2}, index = [0])), ColumnDataSource(pd.DataFrame({"rt": np.median(df.rt.values), "intensity": t3, "mz": m3}, index = [0]))
            p1 = create_p(main_title=f'XIC of monoisotopic theoretical {m1}')
            p2 = create_p(main_title= f'XIC of A+1 theoretical {m2}')
            p3 = create_p(main_title=f'XIC of A+2 theoretical {m3}')
            bottom = b.intensity.min()

            p1.line(x = "rt", y = "intensity", color="firebrick", alpha=0.3, line_width=1, source = ColumnDataSource(b))
            p1.vbar(x = "rt", width = 0.001, bottom = bottom, top = "intensity", color = '#324ea8', fill_alpha = 0, source = c1)
            bottom = a.intensity.min()

            p2.line(x = "rt", y = "intensity", color="firebrick", alpha=0.3, line_width=1, source = ColumnDataSource(a))
            p2.vbar(x = "rt", width = 0.001, bottom = bottom, top = "intensity", color = '#324ea8', fill_alpha = 0, source = c2)
            bottom = c.intensity.min()

            p3.line(x = "rt", y = "intensity", color="firebrick", alpha=0.3, line_width=1, source = ColumnDataSource(c))
            p3.vbar(x = "rt", width = 0.001, bottom = bottom, top = "intensity", color = '#324ea8', fill_alpha = 0, source = c3)

            add_axis_labels(p1)
            add_axis_labels(p2)
            add_axis_labels(p3)

            html = file_html(column(p1, p2, p3), CDN, "my plot")
            self.rt = float(self.rt)
            print(self.rt)
            self.rt = round(self.rt, 3)
            self.rt = str(self.rt)
            f = open(os.path.join(self.temp, f"{self.mf}_{self.ion}_{self.rt}.html"), "w")
            f.write(html)
            f.close()
        elif len(ll) == 2:
            a, b = ll
            self.mzlist = b.mz.values
            t1, t2, t3 = list(self.theoretical_ratio)
            m1, m2, m3 = list(self.theoretical_masses)
            df = ll[0]
            c1, c2 = ColumnDataSource(pd.DataFrame({"rt": np.median(df.rt.values), "intensity": t1, "mz": m1}, index = [0])), ColumnDataSource(pd.DataFrame({"rt": np.median(df.rt.values), "intensity": t2, "mz": m2}, index = [0]))
            p1 = create_p(main_title=f'XIC of monoisotopic theoretical {m1}')
            p2 = create_p(main_title= f'XIC of A+1 theoretical {m2}')
            bottom = b.intensity.min()
            p1.line(x = "rt", y = "intensity", color="firebrick", alpha=0.3, line_width=1, source = ColumnDataSource(b))
            p1.vbar(x = "rt", width = 0.001, bottom = bottom, top = "intensity", color = '#324ea8', fill_alpha = 0, source = c1)
            bottom = a.intensity.min()
            p2.line(x = "rt", y = "intensity", color="firebrick", alpha=0.3, line_width=1, source = ColumnDataSource(a))
            p2.vbar(x = "rt", width = 0.001, bottom = bottom, top = "intensity", color = '#324ea8', fill_alpha = 0, source = c2)
            add_axis_labels(p1)
            add_axis_labels(p2)
            html = file_html(column(p1, p2), CDN, "my plot")
            self.rt = float(self.rt)
            print(self.rt)
            self.rt = round(self.rt, 3)
            self.rt = str(self.rt)
            f = open(os.path.join(self.temp, f"{self.mf}_{self.ion}_{self.rt}.html"), "w")
            f.write(html)
            f.close()
