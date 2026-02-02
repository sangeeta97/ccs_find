from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
import sys, os, re
import pandas as pd
import numpy as np
from collections import defaultdict



def searchButtonStyle():
    return """
    QPushButton{
    background-color: #fcc324;
    border-style:outset;
    border-width:2px;
    border-radius:10px;
    border-color:beige;
    font:12px;
    padding:6px;
    min-width:6em;
    }
     """


def listButtonStyle():
    return """
        QPushButton{
        background-color: #9bc9ff;
        border-style:outset;
        border-width:2px;
        border-radius:10px;
        border-color:beige;
        font:12px;
        padding:6px;
        min-width:6em;
        }

       """

def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """
    base_path = getattr(sys, '_MEIPASS', os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(base_path, relative_path)



class Advanced(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Enter method details")
        self.setWindowIcon(QIcon(resource_path(os.path.join('icon', 'lab.png'))))
        self.setGeometry(100,100,800,600)
        self.information = defaultdict(lambda: None)
        self.UI()
        # self.show()

    def UI(self):
      self.widgets()
      self.c1.textEdited.connect(lambda x: self.onTextEdited(self.c1, "field_strength"))
      self.c2.textEdited.connect(lambda x: self.onTextEdited(self.c2, "drift_tube_pressure"))
      self.c3.textEdited.connect(lambda x: self.onTextEdited(self.c3, "drift_tube_temperature"))
      self.c4.textEdited.connect(lambda x: self.onTextEdited(self.c4, "LC_method_number"))
      self.c5.textEdited.connect(lambda x: self.onTextEdited(self.c5, "processing_applied"))



    def widgets(self):
        self.titleText=QLabel("Please add below parameters for lab information for database")
        self.titleText.setAlignment(Qt.AlignCenter)
        self.submit=QPushButton("Proceed")
        self.cancel=QPushButton("Cancel")
        self.submit.setStyleSheet(searchButtonStyle())
        self.cancel.setStyleSheet(listButtonStyle())
        self.mainLayout=QVBoxLayout()
        self.bottomLayout=QFormLayout()
        self.bottomLayout.addWidget(self.titleText)
        self.c1=QLineEdit()
        self.c1.setText("")
        self.c1.resize(QSize(300, 30))
        self.c2=QLineEdit()
        self.c2.setText("")
        self.c2.resize(QSize(300, 30))
        self.c3=QLineEdit()
        self.c3.setText("")
        self.c3.resize(QSize(300, 30))
        self.c4=QLineEdit()
        self.c4.setText("")
        self.c4.resize(QSize(300, 30))
        self.c5=QLineEdit()
        self.c5.setText("")
        self.c5.resize(QSize(300, 30))
        self.bottomLayout.addRow(QLabel("Field strength (V/cm): "),self.c1)
        self.bottomLayout.addRow(QLabel("Drift tube pressure (Pa) : "),self.c2)
        self.bottomLayout.addRow(QLabel("Drift tube temperature (Kelvin): "),self.c3)
        self.bottomLayout.addRow(QLabel("LC method number: "),self.c4)
        self.bottomLayout.addRow(QLabel("Level of processing applied (e.g., HRdm) : "),self.c5)
        self.bottomLayout.addRow(QLabel(""),self.submit)
        self.bottomLayout.addRow(QLabel(""),self.cancel)
        self.mainLayout.addLayout(self.bottomLayout)
        self.setLayout(self.mainLayout)


    def run(self):
        self.n = 0
        self.submit.clicked.connect(self.Handle)
        self.cancel.clicked.connect(self.nothing)




    def Handle(self):

        if not bool(self.n):
            QMessageBox.information(self, "Information", "Response is submitted")
            self.n += 1




    def nothing(self):
        self.c1.setText("")
        self.c2.setText("")
        self.c3.setText("")
        self.c4.setText("")
        self.c5.setText("")



    def onTextEdited(self, x, y):
        try:
            if x:
                yy = x.text()
                self.information[y] = yy

        except:
            pass
