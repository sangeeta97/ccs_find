import os
import numpy as np
from functools import wraps
import time
from biosaur2.utils import MS1OnlyMzML
import polars as pl
import os
# import modin.config
# modin.config.NPartitions.put(4)
#
# os.environ["MODIN_ENGINE"] = "dask"
import pandas as pd
# import swifter





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


def extract_decorator(row, ps):
    index = row['index']
    rt = row['scanList']['scan'][0]['scan start time']
    rt = round(rt, 3)
    dd = row['scanList']['scan'][0]['ion mobility drift time']
    dd = round(dd, 3)
    # ss = MS1OnlyMzML("20230310_pos_4bit_01.d.DeMP.mzML")
    # ps = list(next(ss).keys())
    mz = row[ps[-2]]
    intensity = row[ps[-1]]

    thresold = 20
    mz = mz[intensity > thresold]
    intensity = intensity[intensity > thresold]

    size = mz.size
    if mz.size > 1:
        all = {"mz": mz, "rt": [rt] * size, "drift_time": [dd] * size, "index": [index] * size, "intensity": intensity}
        df = pl.LazyFrame(all, schema = {"mz": pl.Float32, "rt": pl.Float32, "drift_time": pl.Float32, "index": pl.Float32, "intensity": pl.Float32})
        return df
    return pl.LazyFrame({"mz": None, "rt": None, "drift_time": None, "index": None, "intensity": None})



@timeit
def main_all():
    ss = MS1OnlyMzML("20230310_pos_4bit_01.d.DeMP.mzML")
    ps = list(next(ss).keys())
    print(ps)
    ds = [s for s in ss]
    # dtable = dt.Frame(mz = [], rt = [], drift_time = [], index = [], intensity = [])
    df = pd.DataFrame.from_records(ds)
    # from pandarallel import pandarallel
    # pandarallel.initialize(progress_bar=True)
    tt = df.apply(lambda x: extract_decorator(x, ps), axis = 1)
    dfg = pl.concat(tt.values, how = "vertical")
    # df = pl.from_dicts(ds)
    # ks = df.swifter.apply()
    # print(type(ks))
    # dm = pl.concat(ks, how = "vertical_relaxed")

    # pp = list(next(ss).keys())
    # xx = ([s] for s in ss)
    # df = pd.DataFrame({"spectrum": xx})

    # ss = (extract_decorator(s, pp) for s in ss)
    # import pandas as pd
    # df = pd.DataFrame({"spec": ss})
    # print(df.head())
    # from pandarallel import pandarallel
    # pandarallel.initialize()
    # ps = df["spectrum"].swifter.apply(extract_decorator)
    # ps = (extract_decorator(s, pp) for s in ss)
    # df3 = pl.concat(ps, how = "vertical_relaxed")
    # print(df3.head())





if __name__ == '__main__':
    main_all()
