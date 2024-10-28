"""Antenna Analysis module set up."""

from pathlib import Path
from setuptools import setup, find_packages

_root = Path(__file__).parent

setup(
    name="antenna_analysis",
    version="0.1.0",
    url="https://github.com/b-r-y/AntennaAnalysis",
    author="b-r-y",
    description="Python library for antenna pattern analysis and plotting",
    packages=find_packages(include=["antenna_analysis", "antenna_analysis.*"]),
    install_requires=[
        "openpyxl >= 3.1.5",
        "numpy >= 1.22.1",
        "scipy >= 1.8.0",
        "matplotlib >= 3.4.1",
        "sphinx_rtd_theme >= 3.0.1",
    ],
)
