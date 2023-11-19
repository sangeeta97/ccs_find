import sqlite3
import os
from os import listdir
from os.path import isfile, join


if __name__ == '__main__':
    import sqlite3

    conn = sqlite3.connect('ccs_database.db')
    c = conn.cursor()

    c.execute('''CREATE TABLE ccs_table([id]  INTEGER PRIMARY KEY AUTOINCREMENT,
                                   [pubchemID] text,
                                   [compound_name] text,
                                   [molecular_formula] text,
                                   [ion_type] text,
                                   [mz_measured] real,
                                   [Error_PPM] real,
                                   [conformer] real,
                                   [drift_time] real,
                                   [ccs] real,
                                   [resolving_power] real,
                                   [retention_time] real,
                                   [mf_rt] text,
                                   [time] text
                                   )'''
                                   )



    c.execute('''CREATE TABLE EIC_spectra (
                                          id  INTEGER PRIMARY KEY AUTOINCREMENT,
                                          ion_type text,
                                          molecular_formula text,
                                          rt real,
                                          intensity real,
                                          time text

                                          )'''
              )

    c.execute('''CREATE TABLE DT_spectra (
                                          id  INTEGER PRIMARY KEY AUTOINCREMENT,
                                          ion_type text,
                                          rt real,
                                          molecular_formula text,
                                          intensity real,
                                          drift_time real,
                                          mf_rt text,
                                          time text

                                          )'''
              )

    c.execute('''CREATE TABLE Analytical_Method (
                                          id  INTEGER PRIMARY KEY AUTOINCREMENT,
                                          Analytical_Method text,
                                          Comments text,
                                          time text

                                          )'''
              )

    conn.commit()
