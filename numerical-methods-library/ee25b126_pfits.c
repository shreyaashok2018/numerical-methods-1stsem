/**
 * @file ee25b126_pfits.c
 * @brief Implementation of Interpolation Methods (Lagrange & Newton)
 * @author EE25B126
 * @date November 2025
 */

#include "ee25b126_ee1103.h"
#include <math.h>
#include <stdlib.h>

/* ============================================================================
   5. INTERPOLATION (LAGRANGE & NEWTON)
   ============================================================================ */

double lagrange_interpolate(const double *x, const double *y, int M, double xp) {
    double result = 0.0;
    for (int i = 0; i < M; i++) {
        double term = y[i];
        for (int j = 0; j < M; j++) {
            if (j == i) continue;
            term *= (xp - x[j]) / (x[i] - x[j]);
        }
        result += term;
    }
    return result;
}

void newton_build_table(const double *x, const double *y, int M, double **table) {
    // First column: f[x_i]
    for (int i = 0; i < M; i++) {
        table[i][0] = y[i];
    }

    // Fill divided differences
    for (int j = 1; j < M; j++) {
        for (int i = 0; i < M - j; i++) {
            table[i][j] = (table[i + 1][j - 1] - table[i][j - 1]) / (x[i + j] - x[i]);
        }
    }
}

double newton_interpolate(const double *x, double **table, int M, double xp) {
    double result = table[0][0];
    double term = 1.0;
    for (int j = 1; j < M; j++) {
        term *= (xp - x[j - 1]);
        result += table[0][j] * term;
    }
    return result;
}

double polynomial_eval(double x, const double *coeffs, int N) {
    double result = 0.0;
    double power = 1.0;
    for (int i = 0; i <= N; i++) {
        result += coeffs[i] * power;
        power *= x;
    }
    return result;
}

double gaussian_noise(double stdev) {
    static int have_spare = 0;
    static double spare;

    if (have_spare) {
        have_spare = 0;
        return stdev * spare;
    }

    double u1, u2;
    do {
        u1 = (double)rand() / RAND_MAX;
        u2 = (double)rand() / RAND_MAX;
    } while (u1 <= 0.0);

    double r = sqrt(-2.0 * log(u1));
    double theta = 2.0 * M_PI * u2;
    spare = r * sin(theta);
    have_spare = 1;

    return stdev * r * cos(theta);
}

double calculate_rmse(const double *y_true, const double *y_pred, int N) {
    if (N <= 0) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < N; i++) {
        double diff = y_true[i] - y_pred[i];
        sum += diff * diff;
    }
    return sqrt(sum / N);
}
