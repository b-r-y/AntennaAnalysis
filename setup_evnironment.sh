#!/bin/bash

# Install system dependencies
sudo apt-get install -y texlive texlive-latex-extra python3-venv

python3 -m venv .venv
source .venv/bin/activate

python3 buildall.py