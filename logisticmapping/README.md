# Logistic Map Simulation and Bifurcation Plot (C + Gnuplot)

## Overview

This program simulates the **logistic map**, a classic nonlinear discrete dynamical system, and generates data for a **bifurcation diagram** by varying the control parameter ( r ).

The logistic map is defined as:

[
x_{n+1} = r,x_n(1 - x_n)
]

For each value of ( r ) in a specified range, the system is iterated multiple times, initial transients are discarded, and the remaining values are plotted using **gnuplot**.

---

## Features

* Implements the logistic map in **C**
* Accepts simulation parameters via **command-line arguments**
* Skips the first 100 iterations to remove transient behavior
* Stores data points in a text file for plotting
* Automatically generates a **PNG bifurcation plot** using gnuplot

---

## Files Generated

* `myfile.txt`
  Contains data points in the format:

  ```
  r  x
  ```
* `myfile.png`
  Bifurcation diagram generated using gnuplot

---

## Compilation

Make sure **gcc** and **gnuplot** are installed.

```bash
gcc logistic_map.c -o logistic_map
```

---

## Usage

```bash
./logistic_map rmin rmax rstep Niter
```

### Arguments

| Argument | Description                      |
| -------- | -------------------------------- |
| `rmin`   | Minimum value of parameter ( r ) |
| `rmax`   | Maximum value of parameter ( r ) |
| `rstep`  | Step size for ( r )              |
| `Niter`  | Number of iterations per ( r )   |

### Example

```bash
./logistic_map 0.3 4.0 0.05 300
```

This generates a bifurcation diagram for ( r \in [0.3, 4.0] ).

---

## Methodology

1. Initialize the system with ( x_0 = 0.3 )
2. Iterate the logistic map for `Niter` steps
3. Discard the first 100 iterations to remove transient effects
4. Record subsequent values of ( x_n ) for each ( r )
5. Plot the results using gnuplot

---

## Output Plot

* **x-axis**: control parameter ( r )
* **y-axis**: population variable ( x )
* Each point represents a long-term state of the system

The resulting plot visualizes:

* Fixed points
* Period doubling
* Onset of chaos

---

## Notes & Limitations

* The program assumes valid command-line input (no argument checking).
* A variable-length array is used for storing iterations.
* The transient-skipping logic skips printing, not computation (future improvement).
* Intended as an **educational / exploratory** implementation.

---

## Possible Extensions

* Replace stack array with dynamic memory allocation
* Improve transient-handling logic
* Add argument validation
* Support output to CSV
* Add Lyapunov exponent calculation

Just tell me.

