// ising.c
// 2D Ising Model using Metropolis algorithm with convergence check
// EE25B126, EE25B055, EE25B058

#include "ee25b126_ee1103.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

// Function: neighbor_sum
int neighbor_sum(int **M, int i, int j, int N) {
    int up    = M[(i - 1 + N) % N][j];
    int down  = M[(i + 1) % N][j];
    int left  = M[i][(j - 1 + N) % N];
    int right = M[i][(j + 1) % N];
    return up + down + left + right;
}

/* ----- calculate_magnetization -------------------------------- */
double calculate_magnetization(int **M, int N)
{
    long long sum = 0;
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            sum += M[i][j];
    return (double)sum / (N * N);
}

// Function: alloc_ising_grid
int **alloc_ising_grid(int N) {
    int **grid = (int **)malloc(N * sizeof(int *));
    if (!grid) return NULL;
    for (int i = 0; i < N; i++) {
        grid[i] = (int *)malloc(N * sizeof(int));
        if (!grid[i]) {
            // Clean up previously allocated rows
            for (int k = 0; k < i; k++) free(grid[k]);
            free(grid);
            return NULL;
        }
    }
    return grid;
}

// Function: free_ising_grid
void free_ising_grid(int **grid, int N) {
    if (!grid) return;
    for (int i = 0; i < N; i++) {
        free(grid[i]);
    }
    free(grid);
}

// Function: initialize_random_spins
void initialize_random_spins(int **M, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            M[i][j] = (rand() / (double)RAND_MAX > 0.5) ? 1 : -1;
        }
    }
}


/* ----- metropolis_sweep --------------------------------------- */
void metropolis_sweep(int **Mold, int **Mnew, int N, double J, double T)
{
    /* copy the whole lattice */
    for (int i = 0; i < N; ++i)
        memcpy(Mnew[i], Mold[i], N * sizeof(int));

    /* N*N spin-flip attempts */
    for (int k = 0; k < N * N; ++k) {
        int i = rand() % N;
        int j = rand() % N;

        int s           = Mnew[i][j];
        int sum_neigh   = neighbor_sum(Mnew, i, j, N);
        double delta_E  = 2.0 * J * s * sum_neigh;   /* J can be negative */

        double prob = (delta_E <= 0.0) ? 1.0 : exp(-delta_E / T);
        if ((rand() / (double)RAND_MAX) < prob)
            Mnew[i][j] = -s;
    }
}

// Function: ising_simulation
void ising_simulation(int N, double J, double Tstart, double Tstop,
                      double dT, int max_iterations, double convergence_threshold) {

    // Seed random number generator
    srand(time(NULL));

    // Allocate two lattices
    int **Mold = alloc_ising_grid(N);
    int **Mnew = alloc_ising_grid(N);
    if (!Mold || !Mnew) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(1);
    }

    // Initialize starting lattice with random spins
    initialize_random_spins(Mold, N);

    // Header
    printf("# Temperature   AvgMagnetization   StepsToConvergence\n");

    // Temperature loop
    for (double T = Tstart; T <= Tstop + 1e-9; T += dT) {  // +1e-9 to include Tstop due to floating point
        double M_accum = 0.0;
        int measurements = 0;
        int steps_taken = 0;
        int min_steps = max_iterations / 4;     // Ensure some equilibration
        int equilibration_steps = max_iterations / 2;

        double M_old = 0.0, M_new = 0.0;

        // Copy initial config to Mnew for first iteration
        for (int i = 0; i < N; i++) {
            memcpy(Mnew[i], Mold[i], N * sizeof(int));
        }

        for (int p = 1; p <= max_iterations; p++) {
            steps_taken = p;

            // Perform one Metropolis sweep
            metropolis_sweep(Mold, Mnew, N, J, T);

            // After equilibration, start collecting data
            if (p > equilibration_steps) {
                double M_abs = fabs(calculate_magnetization(Mnew, N));
                M_accum += M_abs;
                measurements++;
            }

            // Convergence check after sufficient measurements
            if (measurements > min_steps) {
                M_old = calculate_magnetization(Mold, N);
                M_new = calculate_magnetization(Mnew, N);
                if (fabs(M_new - M_old) < convergence_threshold) {
                    break;  // Converged
                }
            }

            // Swap lattices for next iteration
            int **temp = Mold;
            Mold = Mnew;
            Mnew = temp;
        }

        // Compute final average
        double final_M_avg = (measurements > 0) ? (M_accum / measurements) : 0.0;

        // Print result only if we measured something
        if (measurements > 0) {
            printf("%f   %f   %d\n", T, final_M_avg, steps_taken);
        } else {
            printf("%f   %f   %d\n", T, 0.0, steps_taken);
        }

        // Carry over final configuration to next temperature
        for (int i = 0; i < N; i++) {
            memcpy(Mold[i], Mnew[i], N * sizeof(int));
        }
    }

    // Cleanup
    free_ising_grid(Mold, N);
    free_ising_grid(Mnew, N);
}
