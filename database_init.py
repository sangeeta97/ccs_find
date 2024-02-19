import sqlite3
import os
from os import listdir
from os.path import isfile, join


if __name__ == '__main__':
    import sqlite3

    conn = sqlite3.connect('ccs_database.db')
    c = conn.cursor()

    c.execute('''CREATE TABLE CCS_table([id]  INTEGER PRIMARY KEY AUTOINCREMENT,
                                   [PubChem_CID] text,
                                   [Compound_name] text,
                                   [Molecular_formula] text,
                                   [Ion_type] text,
                                   [mz_measured] real,
                                   [Error_ppm] real,
                                   [Conformer] real,
                                   [Arrival_time_ms] real,
                                   [CCS_sqA] real,
                                   [Resolving_power] real,
                                   [RT_min] real,
                                   [Formula_RT] text,
                                   [Time] text
                                   )'''
                                   )



    c.execute('''CREATE TABLE EICs (
                                          id  INTEGER PRIMARY KEY AUTOINCREMENT,
                                          Ion_type text,
                                          Molecular_formula text,
                                          RT_min real,
                                          Intensity real,
                                          Time text

                                          )'''
              )


    c.execute('''CREATE TABLE IM_spectra (
                                          id  INTEGER PRIMARY KEY AUTOINCREMENT,
                                          Ion_type text,
                                          RT_min real,
                                          Molecular_formula text,
                                          Intensity real,
                                          Arrival_time_ms real,
                                          Formula_RT text,
                                          Time text

                                          )'''
              )

    c.execute('''CREATE TABLE Analytical_Method (
                                          id  INTEGER PRIMARY KEY AUTOINCREMENT,
                                          Drift_gas text,
                                          beta_slope real,
                                          tfix_intercept real,
                                          Field_strength_V_cm text,
                                          Drift_tube_pressure_Pa text,
                                          Drift_tube_temperature_Kelvin text,
                                          LC_method_number text,
                                          processing_applied text,
                                          Time text

                                          )'''
              )

    conn.commit()
