/**
 * @file ee25b126_normal.c
 * @brief Implementation of Normal Distribution & Statistics functions
 * @author EE25B126
 * @date November 2025
 */

#include "ee25b126_ee1103.h"
#include <math.h>
#include <stdlib.h>
#include <time.h>

/* ============================================================================
   1. NORMAL DISTRIBUTION & STATISTICS
   ============================================================================ */

void seed_random_generator(void) {
    srand((unsigned int)time(NULL));
}

void box_muller_transform(double *arr1, double *arr2, int N,
                         double mean1, double mean2, double stdev) {
    int i = 0;
    while (i < N) {
        double u1 = (double)rand() / RAND_MAX;
        double u2 = (double)rand() / RAND_MAX;

        // Avoid log(0)
        if (u1 <= 0.0) u1 = 1e-12;

        double r = sqrt(-2.0 * log(u1));
        double theta = 2.0 * M_PI * u2;

        double z1 = r * cos(theta);
        double z2 = r * sin(theta);

        if (i < N) arr1[i++] = mean1 + stdev * z1;
        if (i < N) arr1[i++] = mean1 + stdev * z2;
    }

    i = 0;
    while (i < N) {
        double u1 = (double)rand() / RAND_MAX;
        double u2 = (double)rand() / RAND_MAX;

        if (u1 <= 0.0) u1 = 1e-12;

        double r = sqrt(-2.0 * log(u1));
        double theta = 2.0 * M_PI * u2;

        double z1 = r * cos(theta);
        double z2 = r * sin(theta);

        if (i < N) arr2[i++] = mean2 + stdev * z1;
        if (i < N) arr2[i++] = mean2 + stdev * z2;
    }
}

double calculate_mean(const double *arr, int N) {
    if (N <= 0) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < N; i++) {
        sum += arr[i];
    }
    return sum / N;
}

double calculate_variance(const double *arr, int N) {
    if (N <= 0) return 0.0;
    double mean = calculate_mean(arr, N);
    double sum_sq_diff = 0.0;
    for (int i = 0; i < N; i++) {
        double diff = arr[i] - mean;
        sum_sq_diff += diff * diff;
    }
    return sum_sq_diff / N;  // Population variance
}

double calculate_stdev(const double *arr, int N) {
    return sqrt(calculate_variance(arr, N));
}

double t_test(const double *x, const double *y, int N) {
    if (N < 2) return 0.0;

    double mean_x = calculate_mean(x, N);
    double mean_y = calculate_mean(y, N);

    double var_x = calculate_variance(x, N);
    double var_y = calculate_variance(y, N);

    // Convert to sample variance: s² = Σ(x-μ)² / (N-1)
    double s2_x = (N > 1) ? var_x * N / (N - 1.0) : 0.0;
    double s2_y = (N > 1) ? var_y * N / (N - 1.0) : 0.0;

    double pooled_var = (s2_x + s2_y) / 2.0;
    if (pooled_var < EPSILON) return 0.0;

    double se = sqrt(pooled_var * 2.0 / N);
    return (mean_x - mean_y) / se;
}

void shift_to_mean(double *arr, int N, double target_mean) {
    if (N <= 0) return;
    double current_mean = calculate_mean(arr, N);
    double shift = target_mean - current_mean;
    for (int i = 0; i < N; i++) {
        arr[i] += shift;
    }
}

void rescale_stdev(double *arr, int N, double target_stdev) {
    if (N <= 0) return;

    double mean = calculate_mean(arr, N);
    double current_stdev = calculate_stdev(arr, N);

    if (current_stdev < EPSILON) return;  // Avoid division by zero

    double scale = target_stdev / current_stdev;
    for (int i = 0; i < N; i++) {
        arr[i] = mean + (arr[i] - mean) * scale;
    }
}
