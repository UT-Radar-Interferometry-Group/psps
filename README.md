# PSPS
## _Persistent Scatterer selection based on Phase Similarity_
The Phase Similarity Persistent Scatterer (PSPS) is an open-source package for Interferometric Synthetic Aperture Radar (InSAR) data analysis. It reads a stack of interferograms (coregistered and wrapped), and produces a set of persistent scatterer (PS) pixels.

This is research code provided to you "as is" with NO WARRANTIES OF CORRECTNESS. Use at your own risk.

## Installation
A C++ complier is required to install this package.
psps requires Python3.6+ to run.

```
git clone https://github.com/UT-Radar-Interferometry-Group/psps.git
cd psps
pip install -e .
```

## Run a test example
```
cd tests
python test.py
```

## Tutorial in Jupyter Notebook
[How to select PS pixels using the PSPS package](tests/tutorial.ipynb)

[Description on the supported data format](tests/data_format.ipynb)
