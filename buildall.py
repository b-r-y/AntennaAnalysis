#!/usr/bin/env python3
'''
Simple common build script. Builds the documentation in HTLM, latex and PDF formats
'''
import sys
import os

WORKING_DIR = os.getcwd()

os.system("python -m pip install -r requirements.txt")

#os.environ["SPHINX_APIDOC_OPTIONS"]="members,show-inheritance,special-members"

os.chdir("docs")

# Create Autodoc API files
os.system("sphinx-apidoc -f -e -a -T -o ./source/modules ../includes")

# Create the covertrage reports
os.system("sphinx-build -M coverage ./source ./source")

# Build the HTLM Documentation
os.system("sphinx-build -M html ./source ./build")

# Build the Latex Documentation
os.system("sphinx-build -M latex ./source ./build")

os.chdir("build/latex")
os.system("pdflatex AntennaAnalysis.tex")

os.chdir(WORKING_DIR)
