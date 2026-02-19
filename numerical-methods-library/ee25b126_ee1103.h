/**
 * @file ee25b126_ee1103.h
 * @brief Master Header File for EE1103 Library
 * @author EE25B126
 * @date November 2025
 *
 * This header consolidates all functions for the static library,
 * based on the ee1103_master.h template.
 */

#ifndef EE25B126_EE1103_H
#define EE25B126_EE1103_H

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
    STANDARD INCLUDES
   ============================================================================ */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include <stdbool.h> // Added for peak.c

/* ============================================================================
    CONSTANTS AND MACROS
   ============================================================================ */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define MAX_SIZE 256           // Maximum matrix size
#define PBC(i, N) (((i) + (N)) % (N))  // Periodic boundary conditions
#define IDX(i, j, N) ((i) * (N) + (j)) // 2D to 1D array indexing
#define EPSILON 1e-9           // Small value for zero comparisons

/* ============================================================================
    TYPE DEFINITIONS AND STRUCTURES
   ============================================================================ */

/**
 * @brief Function pointer for ODE right-hand side: dy/dx = f(x, y)
 */
typedef double (*ODE_RHS)(double x, double y);

/**
 * @brief 2D Point structure
 */
typedef struct {
    double x;
    double y;
} Point;

/**
 * @brief State vector for ODE solvers (position and velocity)
 */
typedef struct {
    double x, y;     // Position
    double vx, vy;   // Velocity
} State;

/**
 * @brief Peak information for data analysis
 */
typedef struct {
    int peak_index;    // Array index of peak maximum
    double location;   // x-value of peak
    double height;     // y-value (height) of peak
    double width;      // Full Width at Half Maximum (FWHM)
    double separation; // Distance from previous peak
} PeakInfo;

/**
 * @brief Droplet detection result
 */
typedef struct {
    double location;
    double width;
    double separation;
} DropletResult;

/**
 * @brief Hamming distance result
 */
typedef struct {
    int location;    // Bit position
    int distance;    // Hamming distance (0 = perfect match)
} HammingResult;

/* ============================================================================
    1. NORMAL DISTRIBUTION & STATISTICS
   ============================================================================ */

/**
 * @brief Box-Muller transform to generate normally distributed random numbers
 * @param arr1 First output array (will have mean1, stdev)
 * @param arr2 Second output array (will have mean2, stdev)
 * @param N Number of samples
 * @param mean1 Mean for first distribution
 * @param mean2 Mean for second distribution
 * @param stdev Standard deviation for both distributions
 */
void box_muller_transform(double *arr1, double *arr2, int N,
                          double mean1, double mean2, double stdev);

/**
 * @brief Calculate mean (average) of array
 * @param arr Input array
 * @param N Number of elements
 * @return Mean value
 */
double calculate_mean(const double *arr, int N);

/**
 * @brief Calculate variance of array
 * @param arr Input array
 * @param N Number of elements
 * @return Variance
 */
double calculate_variance(const double *arr, int N);

/**
 * @brief Calculate standard deviation of array
 * @param arr Input array
 * @param N Number of elements
 * @return Standard deviation
 */
double calculate_stdev(const double *arr, int N);

/**
 * @brief Student's t-test for two samples
 * @param x First sample array
 * @param y Second sample array
 * @param N Number of elements in each sample
 * @return t-test statistic
 */
double t_test(const double *x, const double *y, int N);

/**
 * @brief Shift array to have specified mean
 * @param arr Array to modify (in-place)
 * @param N Number of elements
 * @param target_mean Desired mean
 */
void shift_to_mean(double *arr, int N, double target_mean);

/**
 * @brief Rescale array to have specified standard deviation
 * @param arr Array to modify (in-place)
 * @param N Number of elements
 * @param target_stdev Desired standard deviation
 */
void rescale_stdev(double *arr, int N, double target_stdev);

/**
 * @brief Initialize random number generator with high-precision seed
 */
void seed_random_generator(void);

/* ============================================================================
    2. TRIGONOMETRIC FUNCTIONS (TAYLOR SERIES)
   ============================================================================ */

/**
 * @brief Calculate factorial
 * @param n Non-negative integer
 * @return n! (factorial of n)
 */
long long factorial(int n);

/**
 * @brief Taylor series approximation of sine
 * @param N Number of terms in series
 * @param theta Angle in radians
 * @return Approximate sin(theta)
 */
double taylor_sine(int N, double theta);

/**
 * @brief Taylor series approximation of cosine
 * @param N Number of terms in series
 * @param theta Angle in radians
 * @return Approximate cos(theta)
 */
double taylor_cosine(int N, double theta);

/**
 * @brief Calculate error between Taylor approximation and true sine
 * @param N Number of terms
 * @param theta Angle in radians
 * @return Absolute error |taylor_sine - sin|
 */
double sine_error(int N, double theta);

/**
 * @brief Calculate error between Taylor approximation and true cosine
 * @param N Number of terms
 * @param theta Angle in radians
 * @return Absolute error |taylor_cosine - cos|
 */
double cosine_error(int N, double theta);

/**
 * @brief Convert degrees to radians
 * @param degrees Angle in degrees
 * @return Angle in radians
 */
double degrees_to_radians(double degrees);

/**
 * @brief Convert radians to degrees
 * @param radians Angle in radians
 * @return Angle in degrees
 */
double radians_to_degrees(double radians);

/* ============================================================================
    3. LOGISTIC MAP
   ============================================================================ */

/**
 * @brief Generate logistic map data
 * @param rmin Minimum growth rate parameter
 * @param rmax Maximum growth rate parameter
 * @param rstep Step size for r
 * @param Niter Number of iterations per r value
 * @param Ntransient Number of initial iterations to discard (burn-in)
 *
 * Outputs to stdout: r, iteration_number, x_value
 * Demonstrates chaos and bifurcation in the logistic equation: x[n+1] = r*x[n]*(1-x[n])
 */
void logistic_map(double rmin, double rmax, double rstep, int Niter, int Ntransient);

/* ============================================================================
    4. ISING MODEL (2D MAGNETIC LATTICE)
   ============================================================================ */

/**
 * @brief Perform one Metropolis-Hastings sweep on Ising lattice
 * @param Mold Current lattice configuration (N×N)
 * @param Mnew Updated lattice configuration (output, N×N)
 * @param N Lattice size (N×N spins)
 * @param J Coupling constant (J<0: ferromagnetic, J>0: antiferromagnetic)
 * @param T Temperature
 */
void metropolis_sweep(int **Mold, int **Mnew, int N, double J, double T);

/**
 * @brief Complete Ising model simulation over temperature range
 * @param N Lattice size (N×N)
 * @param J Coupling constant
 * @param Tstart Starting temperature
 * @param Tstop Ending temperature
 * @param dT Temperature step
 * @param max_iterations Maximum iterations per temperature
 * @param convergence_threshold Magnetization convergence criterion
 *
 * Outputs to file: Temperature, Average Magnetization, Steps to Convergence
 */
void ising_simulation(int N, double J, double Tstart, double Tstop,
                      double dT, int max_iterations, double convergence_threshold);

/**
 * @brief Calculate average magnetization of lattice
 * @param M Lattice configuration (N×N)
 * @param N Lattice size
 * @return Average magnetization (between -1 and +1)
 */
double calculate_magnetization(int **M, int N);

/**
 * @brief Calculate sum of neighboring spins (for energy calculation)
 * @param M Lattice configuration
 * @param i Row index of spin
 * @param j Column index of spin
 * @param N Lattice size
 * @return Sum of 4 nearest neighbors
 */
int neighbor_sum(int **M, int i, int j, int N);

/**
 * @brief Allocate 2D grid for Ising model
 * @param N Lattice size
 * @return Pointer to N×N integer array
 */
int **alloc_ising_grid(int N);

/**
 * @brief Free 2D grid allocated by alloc_ising_grid
 * @param grid Grid to free
 * @param N Lattice size
 */
void free_ising_grid(int **grid, int N);

/**
 * @brief Initialize lattice with random spin configuration
 * @param M Lattice to initialize
 * @param N Lattice size
 */
void initialize_random_spins(int **M, int N);

/* ============================================================================
    5. INTERPOLATION (LAGRANGE & NEWTON)
   ============================================================================ */

/**
 * @brief Lagrange polynomial interpolation
 * @param x Known x-coordinates
 * @param y Known y-coordinates
 * @param M Number of known points
 * @param xp Point at which to evaluate
 * @return Interpolated value at xp
 */
double lagrange_interpolate(const double *x, const double *y, int M, double xp);

/**
 * @brief Build divided difference table for Newton interpolation
 * @param x Known x-coordinates
 * @param y Known y-coordinates
 * @param M Number of points
 * @param table Output divided difference table (M×M)
 *
 * Note: table must be pre-allocated
 */
void newton_build_table(const double *x, const double *y, int M, double **table);

/**
 * @brief Evaluate Newton polynomial using divided difference table
 * @param x Known x-coordinates
 * @param table Divided difference table (from newton_build_table)
 * @param M Number of points
 * @param xp Point at which to evaluate
 * @return Interpolated value at xp
 */
double newton_interpolate(const double *x, double **table, int M, double xp);

/**
 * @brief Evaluate polynomial given coefficients
 * @param x Point at which to evaluate
 * @param coeffs Polynomial coefficients [a0, a1, ..., aN] for a0 + a1*x + ... + aN*x^N
 * @param N Degree of polynomial
 * @return Polynomial value at x
 */
double polynomial_eval(double x, const double *coeffs, int N);

/**
 * @brief Generate Gaussian noise with specified standard deviation
 * @param stdev Standard deviation
 * @return Random value from N(0, stdev²)
 */
double gaussian_noise(double stdev);

/**
 * @brief Calculate Root Mean Square Error between two arrays
 * @param y_true True values
 * @param y_pred Predicted values
 * @param N Number of values
 * @return RMSE = sqrt(mean((y_true - y_pred)²))
 */
double calculate_rmse(const double *y_true, const double *y_pred, int N);

/* ============================================================================
    6. ODE SOLVERS
   ============================================================================ */

/**
 * @brief Euler method ODE solver
 * @param f Right-hand side function dy/dx = f(x,y)
 * @param x0 Initial x value
 * @param y0 Initial y value
 * @param h Step size
 * @param num_steps Number of steps to take
 * @return Final y value after num_steps
 */
double ode_euler(ODE_RHS f, double x0, double y0, double h, int num_steps);

/**
 * @brief Heun's method (improved Euler) ODE solver
 * @param f Right-hand side function dy/dx = f(x,y)
 * @param x0 Initial x value
 * @param y0 Initial y value
 * @param h Step size
 * @param num_steps Number of steps to take
 * @return Final y value after num_steps
 */
double ode_heun(ODE_RHS f, double x0, double y0, double h, int num_steps);

/**
 * @brief Midpoint method (RK2) ODE solver
 * @param f Right-hand side function dy/dx = f(x,y)
 * @param x0 Initial x value
 * @param y0 Initial y value
 * @param h Step size
 * @param num_steps Number of steps to take
 * @return Final y value after num_steps
 */
double ode_midpoint(ODE_RHS f, double x0, double y0, double h, int num_steps);

/**
 * @brief Runge-Kutta 4th order ODE solver
 * @param f Right-hand side function dy/dx = f(x,y)
 * @param x0 Initial x value
 * @param y0 Initial y value
 * @param h Step size
 * @param num_steps Number of steps to take
 * @return Final y value after num_steps
 */
double ode_rk4(ODE_RHS f, double x0, double y0, double h, int num_steps);

/**
 * @brief Single Euler step (updates state in-place)
 * @param theta Current angular position (in-out)
 * @param v Current velocity (in-out)
 * @param a Acceleration
 * @param R Radius
 * @param dt Time step
 */
void euler_step(double *theta, double *v, double a, double R, double dt);

/**
 * @brief Single Heun step (updates state in-place)
 * @param theta Current angular position (in-out)
 * @param v Current velocity (in-out)
 * @param a Acceleration
 * @param R Radius
 * @param dt Time step
 */
void heun_step(double *theta, double *v, double a, double R, double dt);

/**
 * @brief Single midpoint step (updates state in-place)
 * @param theta Current angular position (in-out)
 * @param v Current velocity (in-out)
 * @param a Acceleration
 * @param R Radius
 * @param dt Time step
 */
void midpoint_step(double *theta, double *v, double a, double R, double dt);

/**
 * @brief Single RK4 step (updates state in-place)
 * @param theta Current angular position (in-out)
 * @param v Current velocity (in-out)
 * @param a Acceleration
 * @param R Radius
 * @param dt Time step
 */
void rk4_step(double *theta, double *v, double a, double R, double dt);

/**
 * @brief Normalize angle to [0, 2π)
 * @param angle Angle in radians (any value)
 * @return Equivalent angle in [0, 2π)
 */
double normalize_angle(double angle);

/**
 * @brief Calculate RMS error between two trajectories
 * @param trajectory1 First trajectory
 * @param trajectory2 Second trajectory
 * @param num_points Number of points
 * @return RMS error
 */
double trajectory_rmse(const double *trajectory1, const double *trajectory2, int num_points);

/* ============================================================================
    7. PEAK/DROPLET ANALYSIS
   ============================================================================ */

/**
 * @brief Find peaks in 1D data
 * @param y Data array
 * @param N Number of points
 * @param peak_indices Output array for peak indices
 * @param max_peaks Maximum number of peaks to find
 * @return Number of peaks found
 */
int find_peaks(const double *y, int N, int *peak_indices, int max_peaks);

/**
 * @brief Calculate Full Width at Half Maximum (FWHM)
 * @param x X-coordinates
 * @param y Y-coordinates
 * @param N Number of points
 * @param peak_index Index of peak maximum
 * @return FWHM value
 */
double calculate_fwhm(const double *x, const double *y, int N, int peak_index);

/**
 * @brief Analyze droplets in data file
 * @param filename Input CSV file with x,y data
 * @param results Output array of droplet information
 * @param max_results Maximum number of droplets to detect
 * @return Number of droplets found
 */
int analyze_droplets(const char *filename, DropletResult *results, int max_results);

/**
 * @brief Smooth data using moving average
 * @param input Input data
 * @param output Smoothed output data
 * @param N Number of points
 * @param window_size Smoothing window size (should be odd)
 */
void smooth_data(const double *input, double *output, int N, int window_size);

/**
 * @brief Linear interpolation between two points
 * @param x1 First x coordinate
 * @param y1 First y coordinate
 * @param x2 Second x coordinate
 * @param y2 Second y coordinate
 * @param y_target Target y value
 *Standardize
 * @return Interpolated x value at y_target
 */
double linear_interpolate(double x1, double y1, double x2, double y2, double y_target);

/* ============================================================================
    8. MATRIX OPERATIONS (LINEAR SOLVERS)
   ============================================================================ */

/**
 * @brief Solve linear system using Gaussian Elimination with partial pivoting
 * @param A Augmented matrix [A|b] of size N×(N+1)
 * @param x Solution vector (output)
 * @param N System size
 *
 * Modifies A in-place. Solution returned in x.
 */
void gaussian_elimination(double **A, double *x, int N);

/**
 * @brief Solve linear system using LU Decomposition
 * @param A Augmented matrix [A|b] of size N×(N+1)
 * @param x Solution vector (output)
 * @param N System size
 *
 * Performs LU decomposition with forward/backward substitution.
 */
void lu_decomposition(double **A, double *x, int N);

/**
 * @brief Allocate 2D matrix
 * @param rows Number of rows
 * @param cols Number of columns
 * @return Pointer to rows×cols double array
 */
double **alloc_matrix(int rows, int cols);

/**
 * @brief Free 2D matrix allocated by alloc_matrix
 * @param matrix Matrix to free
 * @param rows Number of rows
 */
void free_matrix(double **matrix, int rows);

/**
 * @brief Read matrix from file
 * @param filename Input file
 * @param N System size
 * @return Augmented matrix [A|b] of size N×(N+1)
 */
double **read_matrix_from_file(const char *filename, int N);

/**
 * @brief Print solution vector
 * @param x Solution vector
 * @param N Size
 */
void print_solution(const double *x, int N);

/**
 * @brief Swap two rows in matrix
 * @param A Matrix
 * @param row1 First row index
 * @param row2 Second row index
 * @param cols Number of columns
 */
void swap_matrix_rows(double **A, int row1, int row2, int cols);

/**
 * @brief Interpolate using natural cubic splines, solving with a matrix function
 * @param x Known x-coordinates (must be sorted ascending)
 * @param y Known y-coordinates
 * @param M Number of known points
 * @param xp Point at which to evaluate
 * @return Interpolated value at xp
 */
double cubic_spline_interpolate(const double *x, const double *y, int M, double xp);

/* ============================================================================
    9. HAMMING DISTANCE & BIT OPERATIONS
   ============================================================================ */

/**
 * @brief Calculate Hamming distance between bit arrays
 * @param bits1 First bit array (packed in unsigned char)
 * @param bits2 Second bit array (packed in unsigned char)
 * @param num_bits Number of bits to compare
 * @return Hamming distance (number of differing bits)
 */
int hamming_distance(const unsigned char *bits1, const unsigned char *bits2, int num_bits);

/**
 * @brief Get bit value at specified position
 * @param arr Bit array (packed in unsigned char)
 * @param bit_index Bit position (0-indexed)
 * @return Bit value (0 or 1)
 */
int get_bit(const unsigned char *arr, int bit_index);

/**
 * @brief Set bit value at specified position
 * @param arr Bit array (packed in unsigned char)
 * @param bit_index Bit position (0-indexed)
 * @param value Bit value to set (0 or 1)
 */
void set_bit(unsigned char *arr, int bit_index, int value);

/**
 * @brief Flip bit at specified position
 * @param arr Bit array (packed in unsigned char)
 * @param bit_index Bit position (0-indexed)
 */
void flip_bit(unsigned char *arr, int bit_index);

/**
 * @brief Generate random bit array
 * @param arr Output bit array (packed in unsigned char)
 * @param num_bits Number of bits to generate
 */
void generate_random_bits(unsigned char *arr, int num_bits);

/**
 * @brief Flip each bit with given probability
 * @param arr Bit array (packed in unsigned char)
 * @param num_bits Number of bits
 * @param probability Probability of flipping each bit [0,1]
 */
void flip_bits_with_probability(unsigned char *arr, int num_bits, double probability);

/**
 * @brief Find codeword in bit stream
 * @param bitstream Bit stream to search (packed in unsigned char)
 * @param stream_length Length of bit stream in bits
 * @param codeword Codeword to find (packed in unsigned char)
 * @param codeword_length Length of codeword in bits
 * @param results Output array of top matches
 * @param max_results Maximum number of results to return
 * @return Number of results found
 *
 * Finds positions with minimum Hamming distance to codeword
 */
int find_codeword(const unsigned char *bitstream, int stream_length,
                  const unsigned char *codeword, int codeword_length,
                  HammingResult *results, int max_results);

/**
 * @brief Convert string of '0' and '1' characters to bit array
 * @param str String containing '0' and '1' characters
 * @param arr Output bit array (packed in unsigned char)
 * @param length Number of bits (length of str)
 */
void string_to_bits(const char *str, unsigned char *arr, int length);

/**
 * @brief Write bit array to binary file
 * @param filename Output filename
 * @param arr Bit array (packed in unsigned char)
 * @param num_bits Number of bits to write
 */
void write_bits_to_file(const char *filename, const unsigned char *arr, int num_bits);

/**
 * @brief Read bit array from binary file
 * @param filename Input filename
 * @param arr Output bit array (packed in unsigned char)
 * @param num_bits Number of bits to read
 */
void read_bits_from_file(const char *filename, unsigned char *arr, int num_bits);

/* ============================================================================
    10. UTILITY FUNCTIONS
   ============================================================================ */

/**
 * @brief Read CSV file with two columns
 * @param filename Input CSV file
 * @param x Output x-coordinates
 * @param y Output y-coordinates
 * @param max_points Maximum number of points to read
 * @return Number of points read
 */
int read_csv_data(const char *filename, double *x, double *y, int max_points);

/**
 * @brief Print array to stdout
 * @param arr Array to print
 * @param N Number of elements
 * @param label Label to print before array
 */
void print_array(const double *arr, int N, const char *label);

// ---------------- EXTRA FUNCTIONS -----------------------------
/**
 * @brief Print bits from a binary file (MSB first)
 * @param filename Binary file path
 * @param num_bits Number of bits to print; -1 = all
 * @return 0 on success, non-zero on error
 */
int print_binary_bits(const char *filename, long num_bits);

/**
 * @brief Plot CSV data (x,y) using GNUplot
 * @param csv_filename Input CSV file (two columns: x,y)
 * @param title Plot title
 * @param xlabel X-axis label
 * @param ylabel Y-axis label
 * @param output_png Output image filename (e.g., "plot.png")
 * @return 0 on success, non-zero on failure
 *
 * Requires GNUplot to be installed.
 */
int plot_csv_data(const char *csv_filename,
                  const char *title,
                  const char *xlabel,
                  const char *ylabel,
                  const char *output_png);

#ifdef __cplusplus
}
#endif

#endif /* EE25B126_EE1103_H */
