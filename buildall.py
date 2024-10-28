#!/usr/bin/env python3
"""Build the documentation in HTLM, latex and PDF formats."""
import sys
import os

if sys.platform == "win32":
    env_keyword = "set"
    os.system("python -m pip install -r requirements.txt")
elif sys.platform == "linux":
    env_keyword = "export"
    os.system("python3 -m pip install -r requirements.txt")

WORKING_DIR = os.getcwd()

# os.environ["SPHINX_APIDOC_OPTIONS"]="members,show-inheritance,special-members"

os.chdir("docs")

# Create Autodoc API files
os.system(f"{env_keyword} PYTHONPATH={WORKING_DIR} && sphinx-apidoc -f -e -a -T -o ./source/modules ../antenna_analysis")

# Create the covertrage reports
os.system(f"{env_keyword} PYTHONPATH={WORKING_DIR} && sphinx-build -M coverage ./source ./source")

# Build the HTLM Documentation
os.system(f"{env_keyword} PYTHONPATH={WORKING_DIR} && sphinx-build -M html ./source ./build")

# Build the Latex Documentation
os.system(f"{env_keyword} PYTHONPATH={WORKING_DIR} && sphinx-build -M latex ./source ./build")

os.chdir("build/latex")
os.system("pdflatex AntennaAnalysis.tex")

os.chdir(WORKING_DIR)
