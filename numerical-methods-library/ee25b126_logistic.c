/**
 * @file ee25b126_logistic.c
 * @brief Implementation of Logistic Map Bifurcation Diagram Generator
 * @author EE25B126
 * @date November 2025
 */

#include "ee25b126_ee1103.h"
#include <stdio.h>
#include <stdlib.h>

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
void logistic_map(double rmin, double rmax, double rstep, int Niter, int Ntransient) {
    // Print header
    printf("r n x\n");

    // Temporary storage for x (only need current and previous)
    double x = 0.5;  // Initial condition

    for (double r = rmin; r <= rmax + EPSILON; r += rstep) {
        x = 0.5;  // Reset initial condition for each r

        // Burn-in (transient) phase - discard first Ntransient iterations
        for (int i = 0; i < Ntransient; i++) {
            x = r * x * (1.0 - x);
        }

        // Output phase - print Niter values after transient
        for (int n = 1; n <= Niter; n++) {
            x = r * x * (1.0 - x);
            printf("%f %d %f\n", r, n, x);
        }
    }
}
