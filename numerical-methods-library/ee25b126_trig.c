/**
 * @file ee25b126_trig.c
 * @brief Implementation of Taylor Series Trigonometric Functions
 * @author EE25B126
 * @date November 2025
 */

#include "ee25b126_ee1103.h"
#include <math.h>

/* ============================================================================
   2. TRIGONOMETRIC FUNCTIONS (TAYLOR SERIES)
   ============================================================================ */

long long factorial(int n) {
    if (n < 0) return 0;
    if (n == 0 || n == 1) return 1;
    long long fact = 1;
    for (int i = 2; i <= n; i++) {
        fact *= i;
    }
    return fact;
}

double taylor_sine(int N, double theta) {
    double sum = 0.0;
    int sign = 1;
    for (int n = 0; n < N; n++) {
        int term = 2 * n + 1;
        long long fact = factorial(term);
        double power = 1.0;
        for (int i = 0; i < term; i++) {
            power *= theta;
        }
        sum += sign * power / fact;
        sign = -sign;
    }
    return sum;
}

double taylor_cosine(int N, double theta) {
    double sum = 0.0;
    int sign = 1;
    for (int n = 0; n < N; n++) {
        int term = 2 * n;
        long long fact = factorial(term);
        double power = 1.0;
        for (int i = 0; i < term; i++) {
            power *= theta;
        }
        sum += sign * power / fact;
        sign = -sign;
    }
    return sum;
}

double sine_error(int N, double theta) {
    double approx = taylor_sine(N, theta);
    double actual = sin(theta);
    return fabs(approx - actual);
}

double cosine_error(int N, double theta) {
    double approx = taylor_cosine(N, theta);
    double actual = cos(theta);
    return fabs(approx - actual);
}

double degrees_to_radians(double degrees) {
    return degrees * M_PI / 180.0;
}

double radians_to_degrees(double radians) {
    return radians * 180.0 / M_PI;
}
