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
        self.setWindowTitle("lab information data Details")
        self.setWindowIcon(QIcon(resource_path(os.path.join('icon', 'lab.png'))))
        self.setGeometry(100,100,800,600)
        self.information = defaultdict(lambda: None)
        self.UI()
        # self.show()

    def UI(self):
      self.widgets()
      self.method.textEdited.connect(lambda x: self.onTextEdited(self.method, "method"))
      self.comments.textEdited.connect(lambda x: self.onTextEdited(self.comments, "comments"))



    def widgets(self):
        self.titleText=QLabel("Please add below parameters for lab information for database")
        self.titleText.setAlignment(Qt.AlignCenter)
        self.method=QLineEdit()
        self.method.setText("")
        self.comments=QLineEdit()
        self.comments.setText("")
        self.comments.resize(QSize(300, 30))
        self.submit=QPushButton("proceed")
        self.cancel=QPushButton("cancel")
        self.submit.setStyleSheet(searchButtonStyle())
        self.cancel.setStyleSheet(listButtonStyle())
        self.mainLayout=QVBoxLayout()
        self.bottomLayout=QFormLayout()
        self.bottomLayout.addWidget(self.titleText)
        self.bottomLayout.addRow(QLabel("Analytical Method: "),self.method)
        self.bottomLayout.addRow(QLabel("Comments: "),self.comments)
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
        self.method.setText("")
        self.comments.setText("")



    def onTextEdited(self, x, y):
        try:
            if x:
                if y == "method":
                    yy = x.text()
                    self.information[y] = yy
                elif y == "comments":
                    yy = x.text()
                    self.information[y] = yy
        except:
            pass
