from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
import sys
from PyQt5.uic import loadUiType
import datetime
from xlrd import *
from xlsxwriter import *
import sqlite3
from collections import defaultdict
import numpy as np
import pandas as pd
import io
import tempfile
import pathlib
import asyncio
from PyQt5.QtCore import QThread
import contextlib
import time
import matplotlib.pyplot as plt
from file_parser import *
import plotly.express as px
from scipy import interpolate
import popup



def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """
    base_path = getattr(sys, '_MEIPASS', os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(base_path, relative_path)



ui,_ = loadUiType(resource_path('ccs.ui'))


p_mass = {"[M+H]+": 1.007276, "[M]+": 0.000548579909, "[M-H]-": 1.007276, "[M]-": 0.000548579909}


theoretical_isotope = {"+electron": '', "-electron": '', "+Cl":  "Cl",  "+Br": "Br", "+CH3COO": "CH3COO", "+HCOO": "HCOO", "+H": "H", "-H": '', "+Na": 'Na', "+K": 'K', "+NH4": 'NH4'}

mass_values = {"+electron": 0.000548579909, "-electron": 0.000548579909, "+Cl":  34.969402,  "+Br": 78.918885, "+CH3COO": 59.013851, "+HCOO": 44.998201, "+H": 1.007276, "-H": 1.007276, "+Na": 22.989218, "+K": 38.963158, "+NH4": 18.033823}

function_values = {"+electron": np.subtract, "-electron": np.add, "+Cl": np.add,  "+Br": np.add, "+CH3COO": np.add, "+HCOO": np.add, "+H": np.add, "-H": np.subtract, "+Na": np.add, "+K": np.add, "+NH4": np.add}




def make_clickable_both(val):
    name, url = val.split('#')
    return f'<a href="{url}">{name}</a>'




class MyAbstract(QThread):
    """Base export thread"""
    done = pyqtSignal(object)
    fail = pyqtSignal(object)
    loop = asyncio.get_event_loop()

    def __init__(self, func, parent=None):
        super().__init__(parent)
        self.func= func

    def run(self):
        try:
            result= self.loop.run_until_complete(self.func())

        except Exception as exception:
            self.fail.emit(exception)
        else:
            self.done.emit(result)


class MainApp(QMainWindow , ui):
    def __init__(self):
        QMainWindow.__init__(self)
        self.setupUi(self)
        self.primary_ion.currentTextChanged.connect(self.primary_ion_text)
        self.drift_gas.currentTextChanged.connect(self.drift_gas_text)
        self.mono_combobox.currentTextChanged.connect(self.mono_combobox_text)
        self.c13_combobox.currentTextChanged.connect(self.c13_combobox_text)
        self.abundance_combobox.currentTextChanged.connect(self.abundance_combobox_text)
        self.upload_mzml.clicked.connect(self.open_file_mzml)
        self.upload_formula.clicked.connect(self.open_file_formula)
        self.upload_calibration.clicked.connect(self.open_file_calibration)
        self.findccs.clicked.connect(self.download3)
        # self.findccs.clicked.connect(self.find_ccs)
        self.tfix.textEdited.connect(self.onTextEdited)
        self.beta.textEdited.connect(self.onTextEdited)
        self.formula.textChanged.connect(self.plaintextdata)
        self.positive_listwidget.itemClicked.connect(self.remove1)
        self.negative_listwidget.itemClicked.connect(self.remove2)

        self.progressbar.setStyle(QStyleFactory.create("windows"))
        self.progressbar.setRange(0, 1)
        self.progress = MyAbstract(self.find_ccs)
        self.delete_row.clicked.connect(self.table_deleterow)
        self.overlay_im.clicked.connect(self.table_xic)
        self.overlay_eic.clicked.connect(self.table_EIC)
        self.all_dt.clicked.connect(self.plot_all_dt)
        self.all_rt.clicked.connect(self.plot_all_rt)
        self.im.clicked.connect(self.table_im)
        self.distribution.clicked.connect(self.table_spectrum)
        self.download.clicked.connect(self.download_result)
        self.plot3d.clicked.connect(self.table_plot3d)
        self.view.clicked.connect(self.view_database)
        self.database.clicked.connect(self.commit_database)
        self.lab.clicked.connect(self.options)

        self.positive_electron.stateChanged.connect(lambda:self.btnstate(self.positive_electron))
        self.positive_hydrogen.stateChanged.connect(lambda:self.btnstate(self.positive_hydrogen))
        self.positive_sodium.stateChanged.connect(lambda:self.btnstate(self.positive_sodium))
        self.positive_potassium.stateChanged.connect(lambda:self.btnstate(self.positive_potassium))
        self.positive_ammonium.stateChanged.connect(lambda:self.btnstate(self.positive_ammonium))

        self.negative_electron.stateChanged.connect(lambda:self.btnstate(self.negative_electron))
        self.negative_hydrogen.stateChanged.connect(lambda:self.btnstate(self.negative_hydrogen))
        self.negative_chlorine.stateChanged.connect(lambda:self.btnstate(self.negative_chlorine))
        self.negative_bromine.stateChanged.connect(lambda:self.btnstate(self.negative_bromine))
        self.negative_formate.stateChanged.connect(lambda:self.btnstate(self.negative_formate))
        self.negative_acetate.stateChanged.connect(lambda:self.btnstate(self.negative_acetate))

        self.positive_pushbutton.clicked.connect(lambda:self.mass_added("positive_pushbutton"))
        self.negative_pushbutton.clicked.connect(lambda:self.mass_added("negative_pushbutton"))
        #
        reg_ex = QRegExp("[A-Z][a-z]?\d*|\(.*?\)\d+")
        input_validator = QRegExpValidator(reg_ex, self.formula)
        self.primary_data = {"primary_ion": "", "drift_gas": "", "mzml": "", "beta": "", "tfix": "", "buffer_text": "", "chargestate": '', "ion_intensity": '', 'use_data': ''}
        self.optional_data = defaultdict(lambda: None)
        self.secondary_data = defaultdict(list)
        self.message = defaultdict(list)
        self.data = defaultdict(lambda: None)
        self.lq = popup.Advanced()


    def options(self):
        self.lq.show()
        self.lq.run()



    @staticmethod
    def write_df_to_qtable(df, table):
        headers = list(df)
        table.setRowCount(df.shape[0])
        table.setColumnCount(df.shape[1])
        table.setHorizontalHeaderLabels(headers)

        # getting data from df is computationally costly so convert it to array first
        df_array = df.values
        for row in range(df.shape[0]):
            for col in range(df.shape[1]):
                table.setItem(row, col, QTableWidgetItem(str(df_array[row,col])))



    def btnstate(self, b):
        if b.isChecked() == True:
            self.secondary_data["checked_ions"].append(b.text())
        elif b.isChecked() == False:
            index = self.secondary_data["checked_ions"].index(b.text())
            self.secondary_data["checked_ions"].pop(index)



    @pyqtSlot()
    def download3(self):

        def fail(exception):
            if not bool(self.n):
                self.n += 1
                QMessageBox.warning(self, "Warning", f"The output not obtained, because of {exception} {str(self.message['warning'])}")



        def done(result):
            if not bool(self.n):
                self.n += 1
                path = os.getcwd()
                path = os.path.join(path, "log.txt")
                f = open(path, "r")
                QMessageBox.information(self, "Information", f"The output result obtained successfully and {f.read()}")
                f.close()


        def started():
            self.progressbar.setRange(0, 0)


        def finished():
            self.progressbar.setRange(0, 1)
            pass


        self.progress.started.connect(started)
        self.progress.finished.connect(finished)
        self.progress.done.connect(done)
        self.progress.fail.connect(fail)
        self.progress.start()




    def onTextEdited(self):
        tfix = self.tfix.text()
        beta = self.beta.text()
        self.primary_data["tfix"] = tfix
        self.primary_data["beta"] = beta



    def mass_added(self, b):
        if b == "positive_pushbutton":
            value = self.positive_lineedit.text()
            value = "+" + value
            self.positive_listwidget.addItem(value)
            self.secondary_data["positive_formula"].append(value)

        if b == "negative_pushbutton":
            value = self.negative_lineedit.text()
            value = "-" + value
            self.negative_listwidget.addItem(value)
            self.secondary_data["negative_formula"].append(value)




    def primary_ion_text(self, s):

        self.primary_data["primary_ion"] = s


    def drift_gas_text(self, s):

        self.primary_data["drift_gas"] = s


    def mono_combobox_text(self, s):
        if s == "PPM Tolerance":
            s = 0.0
        self.secondary_data["mono_combobox"] = s


    def c13_combobox_text(self, s):
        if s == "PPM Tolerance":
            s = 0.0
        self.secondary_data["c13_combobox"] = s

    def abundance_combobox_text(self, s):
        if s == "Percentage tolerance":
            s = 0.0
        self.secondary_data["abundance_combobox"] = s




    def input_text(self):
        if self.chargestate.currentIndex() == 2:
            self.primary_data['chargestate'] = "yes"
        else:
            self.primary_data['chargestate'] = ""

        if self.ion_intensity.currentIndex() == 2:
            self.primary_data['ion_intensity'] = "yes"
        else:
            self.primary_data['ion_intensity'] = ""

        if self.use_data.currentIndex() == 1:
            self.primary_data['use_data'] = "yes"
        else:
            self.primary_data['use_data'] = ""

        if self.peakwidth.currentIndex() == 0:
            self.primary_data['peakwidth'] = 0.5
        else:
            self.primary_data['peakwidth'] = float(self.peakwidth.currentText())




    def file_dialog(self, msg, path):
        return QFileDialog.getOpenFileName(self, msg, path)[0]


    def remove1(self, item):
         value = self.positive_listwidget.currentRow()
         session_name = self.positive_listwidget.currentItem().text()
         self.positive_listwidget.takeItem(value)
         self.secondary_data["positive_formula"].remove(session_name)


    def remove2(self, item):
         value = self.negative_listwidget.currentRow()
         session_name = self.negative_listwidget.currentItem().text()
         self.negative_listwidget.takeItem(value)
         self.secondary_data["negative_formula"].remove(session_name)




    def open_file_formula(self):
        try:

            msg = '1) Select the formula tab (txt) or excel file \n'
            QMessageBox.information(self, 'Add input file', msg)
            sel = 'Select formula file'
            formula = self.file_dialog(sel, ".")
            obj1 = MF_Parser(formula)
            formula_list = obj1.run()
            self.optional_data['formula'] = formula_list
        except:
            QMessageBox.warning(self, "Warning", "Please upload a correct file")
            pass


    def open_file_mzml(self):
        try:

            msg = '1) Select the mzml file \n'
            QMessageBox.information(self, 'Add input file', msg)
            sel = 'Select mzml file'
            mzml = self.file_dialog(sel, ".")
            self.primary_data['mzml'] = mzml
            self.mzlabel.setWordWrap(True)
            _ , mzml_text = os.path.split(mzml)
            self.mzlabel.setText("The uploaded : " + str(mzml_text))


        except:

            QMessageBox.warning(self, "Warning", "Please upload a correct file")
            pass




    def open_file_calibration(self):
        try:

            msg = '1) Select the calibration xml file \n'
            QMessageBox.information(self, 'Add input file', msg)
            sel = 'Select calibration xml file'
            calibration = self.file_dialog(sel, ".")
            obj1 = XML_parser(calibration)
            obj1.run()
            self.optional_data['calibration'] = obj1.result
            self.tfix.setText(obj1.result["TFix"])
            self.beta.setText(obj1.result["Beta"])
            self.primary_data["tfix"] = float(obj1.result["TFix"])
            self.primary_data["beta"] = float(obj1.result["Beta"])

        except:

            QMessageBox.warning(self, "Warning", "Please upload a correct file")
            pass




    def plaintextdata(self):
        text = self.formula.toPlainText()
        formulas = text.split('\n')
        formula_list = [x for x in formulas if bool(x)]
        formula_list = [x.strip() for x in formula_list]
        self.primary_data['buffer_text'] = formula_list



    async def find_ccs(self):
            import gc
            gc.enable()
            gc.collect()
            self.n = 0
            self.result.clear()
            self.input_text()
            self._temp = tempfile.TemporaryDirectory(prefix = "drift_time_")
            self.drift = pathlib.Path(self._temp.name).as_posix()
            self.temp = tempfile.TemporaryDirectory(prefix = "rt_isotopic_")
            self.spectrum = pathlib.Path(self.temp.name).as_posix()
            self.message = defaultdict(list)
            if "mzml" in self.data.keys():
                if not self.primary_data["use_data"]:
                    self.launcher = None
            import src
            self.launcher = src.Final(self.primary_data, self.secondary_data, self.optional_data, mass_values, function_values, self.drift, self.spectrum, self.message, self.data)
            self.ccs_table = self.launcher.run()
            self.ccs_table["compound_name"] = "Unknown"
            self.ccs_table["pubchemID"] = " "
            self.ccs_table["drift_time"] = self.ccs_table["drift_time"].map(lambda x: float(x)).map(lambda x: round(x, 3))
            self.ccs_table["retention_time"] = self.ccs_table["retention_time"].map(lambda x: float(x)).map(lambda x: round(x, 3))
            if "Alert" in self.message.keys():
                self.message["warning"].append(f"The output not obtained, because of {str(self.message['Alert'])}")
                return None
            self.optional_data = defaultdict(lambda: None)
            show_table = self.ccs_table[['pubchemID', 'compound_name', 'molecular_formula',  'ion_type',  'mz_measured', "Error(PPM)", "#conformer", "drift_time", "ccs",  "resolving_power", "retention_time"]]
            self.write_df_to_qtable(show_table, self.result)
            print(self.ccs_table)



    def table_deleterow(self):
        try:

            current_row = self.result.currentRow()
            if current_row < 0:
                return QMessageBox.warning(self, 'Warning','Please select a record to delete')

            button = QMessageBox.question(
                self,
                'Confirmation',
                'Are you sure that you want to delete the selected row?',
                QMessageBox.StandardButton.Yes |
                QMessageBox.StandardButton.No
            )
            if button == QMessageBox.StandardButton.Yes:
                self.result.removeRow(current_row)
                self.ccs_table = self.ccs_table.reset_index(drop = True)
                self.ccs_table = self.ccs_table.drop(current_row)
                self.ccs_table = self.ccs_table.reset_index(drop = True)
        except Exception as e:
            print(e)

            pass






    def table_im(self):
        try:

            row = self.result.currentRow()
            ff = self.ccs_table.loc[row, 'molecular_formula']
            ss = self.ccs_table.loc[row, 'rt']
            ss = str(ss)
            tt = self.ccs_table.loc[row, 'ion_type']
            tt = str(tt)
            import webbrowser
            new = 2
            xx= os.path.join(self.drift, f'{ff}_{tt}_{ss}_IM.html')
            url = f"file://{xx}"
            webbrowser.open(url,new=new)

        except Exception as e:
            print(e)
            pass



    def table_xic(self):
        try:

            row = self.result.currentRow()
            ff = self.ccs_table.loc[row, 'molecular_formula']
            ss = self.ccs_table.loc[row, 'rt']
            ss = str(ss)
            import webbrowser
            new = 2
            xx= os.path.join(self.drift, f'{ff}_{ss}_IMoverlay.html')
            url = f"file://{xx}"
            webbrowser.open(url,new=new)
        except Exception as e:
            print(e)
            pass


    def table_plot3d(self):
        try:

            row = self.result.currentRow()
            ff = self.ccs_table.loc[row, 'molecular_formula']
            tt = self.ccs_table.loc[row, 'ion_type']
            tt = str(tt)
            import webbrowser
            new = 2
            xx= os.path.join(self.spectrum, f'{ff}_{tt}_plot3d.html')
            url = f"file://{xx}"
            webbrowser.open(url,new=new)

        except Exception as e:
            print(e)
            pass



    def table_EIC(self):
        try:

            row = self.result.currentRow()
            ff = self.ccs_table.loc[row, 'molecular_formula']
            import webbrowser
            new = 2
            xx= os.path.join(self.spectrum, f'{ff}_rt_overlay.html')
            url = f"file://{xx}"
            webbrowser.open(url,new=new)
        except Exception as e:
            print(e)
            pass


    def plot_all_rt(self):
        try:
            ll = []
            xx = self.launcher.all_rt
            if len(xx) > 1:
                for x in xx:
                    mm = xx[x]
                    mm = mm[mm['ion_type'] == self.primary_data["primary_ion"]]
                    mm = mm[["rt", "intensity"]]
                    mm['molecular_formula'] = x
                    ll.append(mm)
                df = pd.concat(ll)
                x = df['rt'].values
                y = df['intensity'].values
                f = interpolate.interp1d(x, y, kind = "linear", fill_value = "extrapolate")
                ynew = f(x)
                df['intensity'] = ynew
                fig = px.line(df, x="rt", y="intensity", color='molecular_formula', title=f'overlay_all_rt plot')
                fig.write_html(os.path.join(self.spectrum, f"all_rt_overlay.html"))
                import webbrowser
                new = 2
                xx= os.path.join(self.spectrum, f'all_rt_overlay.html')
                url = f"file://{xx}"
                webbrowser.open(url,new=new)
        except Exception as e:
            print(e)
            pass


    def plot_all_dt(self):
        try:
            ll = []
            xx = self.launcher.all_dt
            if len(xx) > 1:
                for x in xx:
                    mm = xx[x]
                    mm = mm[mm['ion_type'] == self.primary_data["primary_ion"]]
                    mm = mm[["drift_time", "intensity", "mf_rt"]]
                    ll.append(mm)
                df = pd.concat(ll)
                print(df)
                fig = px.line(df, x="drift_time", y="intensity", color='mf_rt', title=f'overlay_all_dt plot')
                fig.write_html(os.path.join(self.spectrum, f"all_dt_overlay.html"))
                import webbrowser
                new = 2
                xx= os.path.join(self.spectrum, 'all_dt_overlay.html')
                url = f"file://{xx}"
                webbrowser.open(url,new=new)
        except Exception as e:
            print(e)
            pass





    def table_spectrum(self):
        try:
            row = self.result.currentRow()
            ff = self.ccs_table.loc[row, 'molecular_formula']
            ss = self.ccs_table.loc[row, 'rt']
            ss = str(ss)
            tt = self.ccs_table.loc[row, 'ion_type']
            tt = str(tt)
            import webbrowser
            new = 2
            xx= os.path.join(self.spectrum, f'{ff}_{tt}_{ss}.html')
            url = f"file://{xx}"
            webbrowser.open(url,new=new)
        except Exception as e:
            print(e)
            pass




    def view_database(self):
        try:

            import sqlite3
            con = sqlite3.connect(resource_path("ccs_database.db"))
            df = pd.read_sql_query("SELECT * from CCS_TABLE", con)
            outtext2 = os.path.join(self.spectrum, 'cal.html')
            fish = open(outtext2, 'w+')
            result = df.to_html()
            fish.write(result)
            fish.close()
            import webbrowser
            new = 2
            xx= os.path.join(self.spectrum, 'cal.html')
            url = f"file://{xx}"
            webbrowser.open(url,new=new)
        except Exception as e:
            print(e)
            pass


    def download_result(self):
        try:

            from os.path import expanduser
            path_text = str(QFileDialog.getExistingDirectory(self, "Select Directory", expanduser("~"), QFileDialog.ShowDirsOnly))

            import string
            import random
            num = random.choice(string.ascii_letters)
            path_text = os.path.join(path_text, f"result_{num}")

            os.mkdir(path_text)
            rowCount = self.result.rowCount()
            columnCount = self.result.columnCount()
            make_all = []
            whole = namedtuple('Whole', ['pubchemID', 'compound_name', 'molecular_formula',  'ion_type',  'mz_measured', "Error_PPM", "conformer", "drift_time", "ccs",  "resolving_power", "retention_time"])
            for row in range(rowCount):
                rowData = []
                for column in range(columnCount):
                    widgetItem = self.result.item(row, column)
                    if (widgetItem and widgetItem.text):
                        rowData.append(widgetItem.text())
                    else:
                        rowData.append('NULL')

                make_all.append(whole._make(rowData))
            df = pd.DataFrame(make_all, columns=whole._fields)
            show_table = self.ccs_table[['pubchemID', 'compound_name', 'molecular_formula',  'ion_type',  'mz_measured', "Error(PPM)", "#conformer", "drift_time", "ccs",  "resolving_power", "retention_time"]]
            df.to_csv(os.path.join(path_text, 'Results.tab'), sep ='\t')
            import shutil
            output = pathlib.Path(self.drift).as_posix()
            output1 = pathlib.Path(self.spectrum).as_posix()
            database_path = resource_path("ccs_database.db")
            shutil.copytree(output, os.path.join(path_text, "drift"))
            shutil.copytree(output1, os.path.join(path_text, "rt_isotope"))
            shutil.copy2(database_path, path_text)
        except Exception as e:
            print(e)
            pass

# https://stackoverflow.com/questions/50209206/clickable-link-in-pandas-dataframe

    def commit_database(self):
        try:

            button = QMessageBox.question(
                    self,
                    'Confirmation',
                    'Are you sure you typed the pubchemID (only numeric) in the table and want to add current data to the database?',
                    QMessageBox.StandardButton.Yes |
                    QMessageBox.StandardButton.No
                )
            if button == QMessageBox.StandardButton.Yes:

                rowCount = self.result.rowCount()
                columnCount = self.result.columnCount()
                make_all = []
                whole = namedtuple('Whole', ['pubchemID', 'compound_name', 'molecular_formula',  'ion_type',  'mz_measured', "Error_PPM", "conformer", "drift_time", "ccs",  "resolving_power", "retention_time"])
                for row in range(rowCount):
                    rowData = []
                    for column in range(columnCount):
                        widgetItem = self.result.item(row, column)
                        if (widgetItem and widgetItem.text):
                            rowData.append(widgetItem.text())
                        else:
                            rowData.append('NULL')

                    mm = whole._make(rowData)
                    make_all.append(mm)
                df = pd.DataFrame(make_all, columns=whole._fields)
                df['mz_measured'] = df['mz_measured'].map(lambda x: float(x))
                df['Error_PPM'] = df['Error_PPM'].map(lambda x: float(x))
                df['conformer'] = df['conformer'].map(lambda x: float(x))
                df['drift_time'] = df['drift_time'].map(lambda x: float(x))
                df['ccs'] = df['ccs'].map(lambda x: float(x))
                df['resolving_power'] = df['resolving_power'].map(lambda x: float(x))
                df["mf_rt"] = df["molecular_formula"] + df["retention_time"].astype(str)
                df['retention_time'] = df['retention_time'].map(lambda x: float(x))
                df["pubchemID"] = df["pubchemID"].map(lambda x: f'https://pubchem.ncbi.nlm.nih.gov/compound/{x}')
                import datetime
                now = datetime.datetime.now()
                df["time"] = now
                conn = sqlite3.connect(resource_path('ccs_database.db'))
                df.to_sql(name='ccs_table', con=conn, if_exists='append', index=False)
                ll = []
                xx = self.launcher.all_rt
                for x in xx:
                    mm = xx[x]
                    mm['molecular_formula'] = x
                    ll.append(mm)
                df2 = pd.concat(ll)
                df2["rt"] = df2["rt"].map(lambda x: float(x)).map(lambda x: round(x, 3))
                df2["intensity"] = df2["intensity"].map(lambda x: float(x)).map(lambda x: round(x, 3))
                ll = []
                xx = self.launcher.all_dt
                for x in xx.keys():
                    mm = xx[x]
                    ll.append(mm)
                df3 = pd.concat(ll)
                df3["rt"] = df3["rt"].map(lambda x: float(x)).map(lambda x: round(x, 3))
                df3["intensity"] = df3["intensity"].map(lambda x: float(x)).map(lambda x: round(x, 3))
                df3["drift_time"] = df3["drift_time"].map(lambda x: float(x)).map(lambda x: round(x, 3))
                df3["mf_rt"] = df3["molecular_formula"].astype(str) + df3["rt"].astype(str)
                df2 = df2[df2["molecular_formula"].isin(df["molecular_formula"])]
                import datetime
                now = datetime.datetime.now()
                df2 = df2[["molecular_formula", "rt", "ion_type", "intensity"]]
                df3 = df3[["molecular_formula", "rt", "ion_type", "intensity", "drift_time", "mf_rt"]]
                df2["time"] = now
                df3["time"] = now
                print(self.lq.information)
                xx = self.lq.information["method"]
                yy = self.lq.information["comments"]
                x1 = self.lq.information["field_strength"]
                x2 = self.lq.information["drift_tube_pressure"]
                x3 = self.lq.information["drift_tube_temperature"]
                x4 = self.lq.information["LC_method_number"]
                x5 = self.lq.information["processing_applied"]
                zz = self.primary_data["drift_gas"]
                kk = float(self.primary_data["beta"])
                gg = float(self.primary_data["tfix"])

                df4 = pd.DataFrame({"Analytical_Method": xx, "Comments": yy, "drift_gas": zz, "beta_slope": kk, "tfix_intercept": gg, "Field_strength_V_cm": x1, "Drift_tube_pressure_Pa": x2, "Drift_tube_temperature_Kelvin": x3, "LC_method_number": x4, "processing_applied_hrDM": x5, "time": now}, index = [0])
                df2.to_sql(name='EIC_spectra', con=conn, if_exists='append', index=False)
                df3.to_sql(name='DT_spectra', con=conn, if_exists='append', index=False)
                df4.to_sql(name='Analytical_Method', con=conn, if_exists='append', index=False)
        except Exception as e:
            print(e)
            pass


def main():
    app = QApplication(sys.argv)
    window = MainApp()
    window.show()
    app.exec_()


if __name__ == '__main__':
    import multiprocessing
    main()
