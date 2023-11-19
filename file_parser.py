import pandas as pd
import numpy as np
import re
import xml.etree.ElementTree as ET
from collections import namedtuple, defaultdict
from lxml import etree
import os
import sys, re
import csv
import re






class MF_Parser:
    def __init__(self, file):
        self.file = file


    def run(self):
        if self.file.endswith((".txt", "tab", "csv")):
            xx = self.text_all()
            return xx
        if self.file.endswith((".xls", ".xlsx")):
            xx = self.excel_all()
            return xx


    def text_all(self):
        with open(self.file, mode= r) as csvfile:
            self.delimit = csv.Sniffer().sniff(csvfile.read(3))
        df = pd.read.table(self.file, sep = self.delimit)
        formulas = list(df.values)
        regex = re.compile("[A-Z][a-z]?\d*|\(.*?\)\d+")
        formula_list = [regex.search(x) for x in formulas if bool(x)]
        formula_list = [x.group() for x in formula_list if x]
        return formula_list


    def excel_all(self):
        df = pd.read_excel(self.file)
        formulas = df.values
        regex = re.compile(r"[A-Z][a-z]?\d*|\((?:[^()]*(?:\(.*\))?[^()]*)+\)\d+")
        formulas = [x[0] for x in formulas if bool(x[0])]
        formula_list = [x for x in formulas if re.match(regex, str(x))]
        return formula_list


#
# DriftGas
# TFix
# Beta

class XML_parser:

    def __init__(self, file):
        self.file = file
        self.tree = ET.parse(self.file)
        self.root = self.tree.getroot()
        self.result = {}



    def run(self):
        for item in self.root.findall('./SingleFieldCcsCalibration/'):
            print(item)
            self.result[item.tag] = item.text
