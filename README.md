# GaussianBeam2DFit

## Overview
GaussianBeam2DFit is an interactive application for analyzing laser beam profiles. It allows users to load camera images, crop regions of interest, and perform two-dimensional Gaussian fits to extract key beam parameters such as beam waists, orientation, and intensity. The tool provides both numerical results and visual plots to make beam diagnostics straightforward.

The application is a Python reimplementation of a MATLAB program, extended with a modern PyQt5 interface and Matplotlib integration.

## Features
- Load and crop beam profile images directly from file.  
- Perform 2D Gaussian fits with optional offset handling.  
- Extract beam parameters: waist sizes, orientation, peak intensity, chi-square fit quality.  
- Display images and fits side by side with adjustable pixel size.  
- View intensity cross-sections along beam axes.  
- Generate beam propagation profiles (w(z)) and extract parameters such as Rayleigh length and M² values.  
- Save results to text files for later analysis.  
- Interactive GUI with navigation toolbars (zoom, pan, save).  

## Requirements
- Python 3.8+  
- PyQt5  
- Matplotlib  
- NumPy  
- SciPy  

Install dependencies with:
```bash
pip install -r requirements.txt
```

## Usage
Run the main application:
```bash
python Gaussianbeam2Dfit.py
```

### Workflow
1. Load a beam image.  
2. Crop to the desired region of interest.  
3. Perform Gaussian fitting to extract parameters.  
4. Add more images at different propagation distances.  
5. Generate intensity profiles or beam propagation plots.  
6. Save the results to text files for documentation or further processing.  

## Motivation
Precise knowledge of beam parameters such as waist size, divergence, and M² is essential for laser optics experiments. This tool provides an open-source, user-friendly alternative to proprietary beam profiler software.


## Acknowledgments
- Based on a MATLAB implementation for Gaussian beam analysis.  
- Uses PyQt5 and Matplotlib for the graphical interface.  
