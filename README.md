# CCSfind
The most recent and up-to date version (February 2024).
We introduce **CCSfind**, a new tool for building of comprehensive databases from experimental LC-IM-MS measurements of small molecules. 
CCSfind allows predicted ion species to be chosen for input chemical formulae.
The tool uses a mzML input files to provide a unified set of results within a single data processing step which includes filtering using isotopologue matching confirmation after peak picking.
CCSfind can handle both chromatographically separated isomers and IM separation of isomeric ions (e.g., “protomers” or conformers of the same ion species) with simple user control over the output for new database entries in SQL format. 

# Compulsory input

* The program needs 5 compulsary inputs.
* A profile data mzML input file from LC-DTIM-MS.
* A list of molecular formulae to be searched for in the mzML file.
* The single-field CCS calibration coefficients (Tfix and Beta) for the data file considered.
* The identity of the drift gas (e.g., nitrogen).
* The primary ion species to be searched for.



# Additional features

* Supports screening of additional adducts ions for both positive and negative modes.
* Option to adjust mass extraction window for isotopologue ion confirmation (A, A+1, A+2).
* Option to adjust the expected LC peak width.
* Option to adjust the expected abundance of the peak (noisy or well-defined peak).
* Option to re-run the analysis faster for an already loaded (parsed) mzML file with different settings or additional molecular formulae.
* Option to search for multiply charged ions.



## Command line tool and create environment using conda (for pip use instructions below), python 3.9 is a must.

```
conda create --name <env> python=3.9
pip install -r requirements.txt
python processor.py

    usage:  python processor.py [-help] 
    [-g config.ini, this file has all parameter to run the analysis, can be changed as per analysis criteria]

     optional arguments:
      -h, --help            show this help message and exit
      --cmd                 activates the command-line interface
      -g                   config.ini file with all the options (refer to config.ini file) (required)

```
## Generating an executable
Using [PyInstaller](http://www.pyinstaller.org) is recommended. You should first clone the repository and install all dependencies.


```
conda config --prepend channels conda-forge
conda create --name ccs_find  python=3.9
conda activate ccs_find
git clone https://github.com/sangeeta97/ccs_find.git
cd ccs_find
pip install -r requirements.txt
python gui.py

```

or download the codes as zip and unpack them and cd into the ccs_find folder and install packages individually name-wise or using requirements.txt

```
pip install polars XlsxWriter xlrd numpy pandas pyarrow PyQt5 molmass matplotlib plotly lxml bokeh biosaur2 scipy
python gui.py
pyinstaller ccs.spec

```

After the following instruction, the directory `dist` will be created (among others) and the executable will be inside it:
```
ccs_find.exe
```


# How to use it


* To use it as GUI tool; install on your terminal and follow the instructions.
```
python gui.py
```


