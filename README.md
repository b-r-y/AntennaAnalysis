# Prepare the environment and prerequisites

This is best run in a python virtual environment.

The external packages needed (beyond the python module dependencies) are:

* a TeX compiler
* python virtual environment

On Linux this can be done with the `setup_evnironment.sh` script.

# Install

To install, navigate to the cloned folder and run:

``` bash
python -m pip install .
```

# Uninstall 

To uninstall run:

```bash
python -m pip uninstall antenna_analysis
```

# Running examples

The `example` folder contains a number of examples for typical usage.
To use them the module must be installed.
For trying things out without isntalling the module check the `test` folder, where the module is loaded on the fly.

# Getting help

From the cloned folder run:

```bash
python ./buildall.py
```
Then either:

* open in browser ./docs/build/html/index.html 
    
OR     
   
* open in PDF viewer ./docs/build/latex/AntennaAnalysis.PDF

NOTE: 
To be albe to build PDF version of the documentation a TEX compiler is required.
The recommended one is [TeX Live](https://tug.org/texlive/).
