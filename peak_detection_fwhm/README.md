Peak Detection and FWHM Analysis in Time-Series Data

Language: C
Domain: Signal Processing / Numerical Methods
Type: Offline Data Analysis Tool
Authors: Shreya Ashok, Dhruvadeepa Reddy, Sumaiya Fatima (As part of numerical methods course)

Overview
This project implements an automated peak detection and characterisation system for time-series sensor data.
Given an input CSV file containing (x, y) data points, the program:

-Detects peak regions using threshold-based segmentation
-Fits a quadratic curve to each detected peak
-Estimates peak position and height analytically
-Identifies rising and falling slopes using linear regression
-Computes Full Width at Half Maximum (FWHM)
-Outputs peak position, FWHM, and inter-peak separation

The goal was to design a structured algorithm for extracting meaningful physical parameters from noisy experimental data.

Input:
CSV file expected
Skips the first line.
Takes 2 columns

Algorithm Overview
1. Data Loading
Dynamically allocates memory for large datasets
Resizes arrays using realloc() when needed
Parses input file into (x, y) data points

2. Peak Region Detection
Peak regions are identified using:
Threshold condition (y > 1.5)
Minimum number of data points per region
Separation constraint between peaks
Region stagnation detection to determine peak boundary
This isolates candidate peak segments in the signal.

3. Quadratic Least-Squares Peak Fitting
A quadratic model is fitted to each selected region using least squares regression with summation.
Peak location and height are computed from the quadratic formulae

4. Rising and Falling Slope Detection
Linear regression is applied to estimate local slopes.

Slopes are classified as:
Rising (positive slope)
Falling (negative slope)

Only slopes exceeding a minimum threshold are retained.

5. FWHM Calculation
For each valid peak:
Compute half maximum: y_half = y_peak / 2
Intersect rising and falling linear fits with y_half
Compute: FWHM = x₂ − x₁
Constraints are applied to reject unrealistic matches.

Output
x_peak      FWHM        separation
where x_peak is the location of the peak, fwhm is the full width at half maximum and separation is the distance from the previous valid peak 




Key Concepts Demonstrated
-Time-series signal segmentation
-Threshold-based peak detection
-Quadratic least-squares regression
-Linear regression for slope estimation
-Analytical extraction of peak characteristics
-Dynamic memory allocation in C
-Numerical stability handling

Design Considerations
-Peaks are filtered by minimum inter-peak distance
-Slopes are matched using proximity constraints
-Singular matrix checks prevent invalid regression results
-Memory is dynamically resized for scalability

Limitations
This code was written for a pre-defined experimental dataset, hence assumptions with respect to threshold and quadratic peaks were made for simplicity, which reduces reusability
-Fixed threshold values (not adaptive)
-Offline processing only
-Assumes peaks are approximately quadratic
-No smoothing or filtering stage prior to detection
-Computational complexity increases with dataset size

Possible Extensions
-Adaptive thresholding
-Noise filtering (moving average / low-pass filtering)
-Gaussian peak fitting instead of quadratic
-Vectorised or optimized implementation
-Visualization module

Learning Outcomes
This project explores:
-Numerical regression techniques
-Peak characterisation methods
-Signal feature extraction
-Structured algorithm design in C
-Handling real-world noisy datasets
