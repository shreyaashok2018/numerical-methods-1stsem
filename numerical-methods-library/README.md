EE1103 Numerical Methods Library (Static C Library)

Overview
This project is a modular static library implementing core numerical methods and computational physics algorithms in C.

It was designed as a reusable scientific computing toolkit with clean header separation and Makefile-based build automation.

Features
Numerical Methods
-Euler, Heun, Midpoint, RK4 ODE solvers
-Logistic map simulation (chaotic dynamics)
-Lagrange & Newton interpolation
-Natural cubic spline interpolation
-Gaussian elimination & LU decomposition

Computational Physics
-2D Ising model (Metropolis-Hastings Monte Carlo)
-Magnetization computation
-Periodic boundary conditions

Signal & Data Analysis
-Peak detection
-Full Width at Half Maximum (FWHM)
-Droplet analysis
-Moving average smoothing

Statistics
-Box-Muller normal distribution generator
-Mean, variance, standard deviation
-Student’s t-test

Bit & Communication Utilities
-Hamming distance calculation
-Bit flipping with probability
-Codeword detection in bitstreams

Build Instructions
make

Generates:
ee25b126_ee1103.a

Optional:
make lib

Creates:
libee25b126_ee1103.a

Clean:
make clean

Design Philosophy
-Modular architecture
-Reusable header interface
-Explicit memory management
-Numerical stability checks
-Minimal external dependencies

Concepts Covered
-Numerical integration
-Monte Carlo simulation
-Chaos theory
-Linear algebra solvers
-Signal peak detection
-Statistical analysis
-Bit-level error detection
-Static library compilation
