# PSPS
## _Persistent Scatterer selection based on Phase Similarity_
The Phase Similarity Persistent Scatterer (PSPS) is an open-source package for Interferometric Synthetic Aperture Radar (InSAR) data analysis. It reads a stack of interferograms (coregistered and wrapped), and produces a set of persistent scatterer (PS) pixels.

This is research code provided to you "as is" with NO WARRANTIES OF CORRECTNESS. Use at your own risk.

## Installation
PSPS requires Python3.6+ to run.
Clone to build/install
```
pip install -e .
```
To obtain the pybind11 headers (included as a submodule)

`git submodule update --init --recursive`

https://pybind11.readthedocs.io/en/stable/installing.html#include-as-a-submodule

These will be located in the `extern` folder

## Tutorial in Jupyter Notebook
[A tutorial](tests/tutorial.ipynb) is provided to demonstrate how to select PS pixels using the PSPS package.
