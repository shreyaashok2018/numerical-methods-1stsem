/**
 * @file ee25b126_ode.c
 * @brief Implementation of Generic ODE Solvers and Pendulum Single-Step Updates
 * @author EE25B126
 * @date November 2025
 */

#include "ee25b126_ee1103.h"
#include <math.h>

/* ============================================================================
   6. ODE SOLVERS - GENERIC (1st-order: dy/dx = f(x,y))
   ============================================================================ */

double ode_euler(ODE_RHS f, double x0, double y0, double h, int num_steps) {
    double x = x0;
    double y = y0;
    for (int i = 0; i < num_steps; i++) {
        y += h * f(x, y);
        x += h;
    }
    return y;
}

double ode_heun(ODE_RHS f, double x0, double y0, double h, int num_steps) {
    double x = x0;
    double y = y0;
    for (int i = 0; i < num_steps; i++) {
        double k1 = f(x, y);
        double k2 = f(x + h, y + h * k1);
        y += h * (k1 + k2) / 2.0;
        x += h;
    }
    return y;
}

double ode_midpoint(ODE_RHS f, double x0, double y0, double h, int num_steps) {
    double x = x0;
    double y = y0;
    for (int i = 0; i < num_steps; i++) {
        double k1 = f(x, y);
        double k2 = f(x + h/2.0, y + h/2.0 * k1);
        y += h * k2;
        x += h;
    }
    return y;
}

double ode_rk4(ODE_RHS f, double x0, double y0, double h, int num_steps) {
    double x = x0;
    double y = y0;
    for (int i = 0; i < num_steps; i++) {
        double k1 = f(x, y);
        double k2 = f(x + h/2.0, y + h/2.0 * k1);
        double k3 = f(x + h/2.0, y + h/2.0 * k2);
        double k4 = f(x + h, y + h * k3);
        y += h * (k1 + 2.0*k2 + 2.0*k3 + k4) / 6.0;
        x += h;
    }
    return y;
}

/* ============================================================================
   PENDULUM SINGLE-STEP (Second-order: d²θ/dt² = a/R)
   ============================================================================ */

void euler_step(double *theta, double *v, double a, double R, double dt) {
    *theta += *v * dt;
    *v += (a / R) * dt;
}

void heun_step(double *theta, double *v, double a, double R, double dt) {
    double k1_theta = *v;
    double k1_v = a / R;

    double v_mid = *v + k1_v * dt;

    double k2_theta = v_mid;
    double k2_v = a / R;

    *theta += (k1_theta + k2_theta) * dt / 2.0;
    *v += (k1_v + k2_v) * dt / 2.0;
}

void midpoint_step(double *theta, double *v, double a, double R, double dt) {
    double k1_theta = *v;
    double k1_v = a / R;

    double theta_mid = *theta + k1_theta * dt / 2.0;
    double v_mid = *v + k1_v * dt / 2.0;

    *theta += v_mid * dt;
    *v += (a / R) * dt;
}

void rk4_step(double *theta, double *v, double a, double R, double dt) {
    double accel = a / R;

    double k1_theta = *v;
    double k1_v = accel;

    double v2 = *v + k1_v * dt / 2.0;
    double k2_theta = v2;
    double k2_v = accel;

    double theta3 = *theta + k2_theta * dt / 2.0;
    double v3 = *v + k2_v * dt / 2.0;
    double k3_theta = v3;
    double k3_v = accel;

    double theta4 = *theta + k3_theta * dt;
    double v4 = *v + k3_v * dt;
    double k4_theta = v4;
    double k4_v = accel;

    *theta += (k1_theta + 2.0*k2_theta + 2.0*k3_theta + k4_theta) * dt / 6.0;
    *v += (k1_v + 2.0*k2_v + 2.0*k3_v + k4_v) * dt / 6.0;
}

/* ============================================================================
   UTILITIES
   ============================================================================ */

double normalize_angle(double angle) {
    double two_pi = 2.0 * M_PI;
    angle = fmod(angle, two_pi);
    if (angle < 0.0) angle += two_pi;
    return angle;
}

double trajectory_rmse(const double *trajectory1, const double *trajectory2, int num_points) {
    if (num_points <= 0) return 0.0;
    double sum_sq = 0.0;
    for (int i = 0; i < num_points; i++) {
        double diff = trajectory1[i] - trajectory2[i];
        sum_sq += diff * diff;
    }
    return sqrt(sum_sq / num_points);
}
