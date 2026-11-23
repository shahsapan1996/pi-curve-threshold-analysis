# PI Curve Threshold Analysis (MATLAB)

Batch analysis of laser power–current (P–I) measurements. Determines threshold current and slope efficiency for several files at once.

## Methods
- **Linear fit** above the threshold knee (threshold from the x-intercept; uncertainty from varying the fit window)
- **Second derivative (d²P/dI²)** peak in the knee region

## Features
- Multi-file selection, optional splitting of sweeps by a column
- Dark-offset correction and Savitzky-Golay smoothing
- Figures saved as PNG, summary saved as `PI_threshold_summary.csv`

## Requirements
MATLAB R2020a or newer

## Usage
Edit the settings at the top (column numbers, units, fit window), run the script and select your data files.

## Input format
Columns for current and power (default: columns 1 and 2). Units are set with `unitScale`.
