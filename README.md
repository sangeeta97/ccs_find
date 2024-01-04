# ccs_find
The most recent and up-to date (Jan 2024).
We introduce CCSfind, a new tool for building of comprehensive databases from experimental IM-MS measurements of small molecules. 
CCSfind allows predicted ion species to be chosen for input chemical formulae.
The tool uses a mzML input files to provide a unified set of results within a single data processing step which includes filtering using isotopic confirmation after peak picking.
CCSfind can handle both chromatographically separated isomers and IM separation of isomeric ions (e.g., “protomers” or conformers of the same ion species) with 
simple user control over the output for new database entries in SQL format. 

# Compulsory Input

* Then the program needs compulsory 5 inputs to start the analysis.
* A profile mzml input file
* A list of molecular formulas to be screened in the mzml file
* The Tfix and Beta values of calibration of the experimental analysis
* The nature of gas
* The primary ion for screening (either positive or negative primary ion)



# Optional or additional inputs

* support for screening of additional adducts ions including both positive and negative modes
* Option to select mass screening window for isotopic ions confirmation (A, A+1, A+2)
* Option to select the expected peakwidth of the peak 
* Option to selct abundance of the peak (noisy or well-defined peak)
* Option to rerun the analysis faster for already loaded mzml file with different setting or new molecular formulas
* Option to search for multiple charged ions



## Commandline Tool and create environment using conda (for pip use instructions below)

```
conda create --name <env> python=3.9
pip install -r requirements.txt
python processor.py

    usage:  python processor.py [-help] 
    [-g config.ini, this file has all parameter to run the analysis, can be changed as per analysis criteria]

     optional arguments:
      -h, --help            show this help message and exit
      --cmd                 activates the command-line interface
      -g                   config.ini file with all the options (refer to config.ini file)

```
## Generating an executable
Using [PyInstaller](http://www.pyinstaller.org) is recommended. You should first clone the repository and install all dependencies.


```
conda config --prepend channels conda-forge
conda create --name ccs_find  python=3.9
pip install polars XlsxWriter xlrd numpy pandas pyarrow PyQt5 molmass matplotlib plotly lxml bokeh biosaur2 scipy
conda activate ccs_find
git clone https://github.com/sangeeta97/ccs_find.git
cd ccs_find
pyinstaller ccs.spec

```

After the following instruction, the directory `dist` will be created (among others) and the executable will be inside it:
```
ccs_find.exe
```


# How to use it


* To use it as GUI tool; Please on your terminal and follow the instructions.
```
python gui.py
```


