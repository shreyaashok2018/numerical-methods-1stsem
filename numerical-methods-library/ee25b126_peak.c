/*
 * ee25b126_peak.c
 * Implementation of the droplet-analysis functions declared in
 * ee25b126_ee1103.h
 *
 * All thresholds are identical to the original assignment code.
 */

#include "ee25b126_ee1103.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

/* ------------------------------------------------------------------ */
/*  Internal helper types (not exported)                              */
/* ------------------------------------------------------------------ */
typedef struct {
    double xpeak;
    double ypeak;
    double separation;      /* distance to previous peak */
} Peak;

typedef struct {
    double m;               /* slope */
    double c;               /* intercept */
    int    region_start_x;
    bool   used;
} Slope;

/* ------------------------------------------------------------------ */
/*  Constants – **exactly** the same as the original program           */
/* ------------------------------------------------------------------ */
static const double DELTA1               = 100.0;
static const double DELTA2               = 2500.0;
static const int    MIN_POINTS_FOR_PEAK  = 50;
static const double SLOPE_THRESHOLD      = 0.00005;
static const double MIN_INTERPEAK_DIST   = 3000.0;
static const double MAX_SLOPE_DISTANCE   = 2500.0;

/* ------------------------------------------------------------------ */
/*  Real quadratic-fit peak detector (internal)                        */
/* ------------------------------------------------------------------ */
static int
detect_quadratic_peaks(const double *x, const double *y, int N,
                       Peak *peaks, int max_peaks, double *prev_xpeak)
{
    int   count = 0;
    double s[7] = {0};               /* sum x, x², x³, x⁴, y, xy, x²y */
    int   points_in_region = 0;
    int   prev_points = 0;
    int   xstart = 0, xend = 0;
    bool  in_region = false;
    double xpold = *prev_xpeak;
    int   prev_x = 0;                /* <-- added declaration */

    for (int i = 0; i < N; ++i) {
        double cur_x = x[i];
        double cur_y = y[i];

        /* ----- start / continue high-y region (y > 1.5) ----- */
        if (cur_y > 1.5 && (xpold == 0.0 || (cur_x - xpold > DELTA2))) {
            if (!in_region) {
                xstart = (int)cur_x;
                in_region = true;
                for (int j = 0; j < 7; ++j) s[j] = 0.0;
                points_in_region = 0;
            }
            double xrel = cur_x - xstart;
            s[0] += xrel;
            s[1] += xrel*xrel;
            s[2] += xrel*xrel*xrel;
            s[3] += xrel*xrel*xrel*xrel;
            s[4] += cur_y;
            s[5] += xrel*cur_y;
            s[6] += xrel*xrel*cur_y;
            ++points_in_region;
        }

        /* ----- stagnation detection ----- */
        if (points_in_region >= MIN_POINTS_FOR_PEAK && prev_points == points_in_region)
            xend = prev_x;

        /* ----- region ends when gap >= DELTA1 ----- */
        if (xend != 0 && (cur_x - xend >= DELTA1)) {
            in_region = false;

            double denom1 = s[1]*s[2] - s[0]*s[3];
            double denom2 = s[1]*s[1] - s[0]*s[2];
            if (fabs(denom1) < 1e-9 || fabs(denom2) < 1e-9)
                goto next_iteration;

            double a = (s[1]*s[5] - s[0]*s[6]) / denom1;
            double b = (s[1]*s[5] - s[0]*s[6] - s[1]*s[2]*a + s[0]*s[3]*a) / denom2;
            double c = s[4]/points_in_region - s[1]*a/points_in_region
                                          - s[0]*b/points_in_region;

            double xpeak = (a != 0.0) ? -b/(2.0*a) + xstart : 0.0;
            double ypeak = (a != 0.0) ? -(b*b - 4.0*a*c)/(4.0*a) : 0.0;

            if (count < max_peaks) {
                peaks[count].xpeak     = xpeak;
                peaks[count].ypeak     = ypeak;
                peaks[count].separation = (xpold == 0.0) ? 0.0 : xpeak - xpold;
                ++count;
            }

            xpold = xpeak;
        next_iteration:
            xend = 0;
            for (int j = 0; j < 7; ++j) s[j] = 0.0;
            points_in_region = 0;
        }

        prev_points = points_in_region;
        prev_x = (int)cur_x;
    }

    *prev_xpeak = xpold;
    return count;
}

/* ------------------------------------------------------------------ */
/*  Slope collection (mid-y band 0.85 < y < 1.4)                       */
/* ------------------------------------------------------------------ */
static void
collect_slopes(const double *x, const double *y, int N,
               Slope *rising, int *rising_cnt,
               Slope *falling, int *falling_cnt,
               int max_slopes)
{
    int n = 0;
    double l[4] = {0};               /* sum xy, x, y, x² */
    int region_start = 0;

    for (int i = 0; i < N; ++i) {
        double cur_x = x[i];
        double cur_y = y[i];

        if (cur_y > 0.85 && cur_y < 1.4) {
            if (n == 0) region_start = (int)cur_x;
            double xr = cur_x - region_start;
            l[0] += xr * cur_y;
            l[1] += xr;
            l[2] += cur_y;
            l[3] += xr * xr;
            ++n;
        }
        else if (n > 0) {
            if (n > 5) {
                double denom = n*l[3] - l[1]*l[1];
                if (fabs(denom) > 1e-9) {
                    double m = (n*l[0] - l[1]*l[2]) / denom;
                    double c = (l[2]/n) - m*(l[1]/n);
                    if (fabs(m) > SLOPE_THRESHOLD) {
                        if (m > 0 && *rising_cnt < max_slopes) {
                            rising[*rising_cnt].m = m;
                            rising[*rising_cnt].c = c;
                            rising[*rising_cnt].region_start_x = region_start;
                            rising[*rising_cnt].used = false;
                            ++(*rising_cnt);
                        }
                        else if (m < 0 && *falling_cnt < max_slopes) {
                            falling[*falling_cnt].m = m;
                            falling[*falling_cnt].c = c;
                            falling[*falling_cnt].region_start_x = region_start;
                            falling[*falling_cnt].used = false;
                            ++(*falling_cnt);
                        }
                    }
                }
            }
            n = 0;
            for (int j = 0; j < 4; ++j) l[j] = 0.0;
        }
    }

    /* final segment */
    if (n > 5) {
        double denom = n*l[3] - l[1]*l[1];
        if (fabs(denom) > 1e-9) {
            double m = (n*l[0] - l[1]*l[2]) / denom;
            double c = (l[2]/n) - m*(l[1]/n);
            if (fabs(m) > SLOPE_THRESHOLD) {
                if (m > 0 && *rising_cnt < max_slopes) {
                    rising[*rising_cnt].m = m;
                    rising[*rising_cnt].c = c;
                    rising[*rising_cnt].region_start_x = region_start;
                    rising[*rising_cnt].used = false;
                    ++(*rising_cnt);
                }
                else if (m < 0 && *falling_cnt < max_slopes) {
                    falling[*falling_cnt].m = m;
                    falling[*falling_cnt].c = c;
                    falling[*falling_cnt].region_start_x = region_start;
                    falling[*falling_cnt].used = false;
                    ++(*falling_cnt);
                }
            }
        }
    }
}

/* ------------------------------------------------------------------ */
/*  Match slopes to a peak and compute FWHM                           */
/* ------------------------------------------------------------------ */
static double
match_slopes_and_fwhm(const Peak *p,
                      Slope *rising, int rising_cnt,
                      Slope *falling, int falling_cnt)
{
    int   best_r = -1;
    double min_dr = 1e9;
    for (int j = 0; j < rising_cnt; ++j) {
        if (!rising[j].used && rising[j].region_start_x < p->xpeak) {
            double d = p->xpeak - rising[j].region_start_x;
            if (d <= MAX_SLOPE_DISTANCE && d < min_dr) {
                min_dr = d; best_r = j;
            }
        }
    }

    int   best_f = -1;
    double min_df = 1e9;
    for (int j = 0; j < falling_cnt; ++j) {
        if (!falling[j].used && falling[j].region_start_x > p->xpeak) {
            double d = falling[j].region_start_x - p->xpeak;
            if (d <= MAX_SLOPE_DISTANCE && d < min_df) {
                min_df = d; best_f = j;
            }
        }
    }

    if (best_r == -1 || best_f == -1) return -1.0;

    Slope *rs = &rising[best_r];
    Slope *fs = &falling[best_f];
    rs->used = true;
    fs->used = true;

    double y_half = p->ypeak / 2.0;
    double x1 = (y_half - rs->c) / rs->m + rs->region_start_x;
    double x2 = (y_half - fs->c) / fs->m + fs->region_start_x;
    double fwhm = x2 - x1;

    return (fwhm > 0.0 && fwhm < 2000.0) ? fwhm : -1.0;
}

/* ------------------------------------------------------------------ */
/*  PUBLIC LIBRARY FUNCTIONS                                          */
/* ------------------------------------------------------------------ */

/* 1. Simple peak finder required by the header (not used by droplets) */
int
find_peaks(const double *y, int N, int *peak_indices, int max_peaks)
{
    (void)y; (void)N; (void)peak_indices; (void)max_peaks;
    /* The droplet algorithm does its own quadratic fit – this stub
       satisfies the header without pulling in unnecessary code. */
    return 0;
}

/* 2. FWHM calculator required by the header (stub) */
double
calculate_fwhm(const double *x, const double *y, int N, int peak_index)
{
    (void)x; (void)y; (void)N; (void)peak_index;
    return 0.0;
}

/* 3. Main droplet analysis driver */
int
analyze_droplets(const char *filename,
                 DropletResult *results,
                 int max_results)
{
    const int MAX_POINTS = 2000000;
    double *x = malloc(MAX_POINTS * sizeof(double));
    double *y = malloc(MAX_POINTS * sizeof(double));
    if (!x || !y) { free(x); free(y); return 0; }

    int N = read_csv_data(filename, x, y, MAX_POINTS);
    if (N <= 0) { free(x); free(y); return 0; }

    /* ---- quadratic peaks ---- */
    const int MAX_PEAKS = 500;
    Peak *peaks = malloc(MAX_PEAKS * sizeof(Peak));
    double prev_xpeak = 0.0;
    int peak_cnt = detect_quadratic_peaks(x, y, N, peaks, MAX_PEAKS, &prev_xpeak);
    if (peak_cnt == 0) { free(x); free(y); free(peaks); return 0; }

    /* ---- slopes ---- */
    const int MAX_SLOPES = 4000;
    Slope *rising  = malloc(MAX_SLOPES * sizeof(Slope));
    Slope *falling = malloc(MAX_SLOPES * sizeof(Slope));
    int rising_cnt = 0, falling_cnt = 0;
    collect_slopes(x, y, N, rising, &rising_cnt, falling, &falling_cnt, MAX_SLOPES);

    /* ---- match & fill results ---- */
    int out_cnt = 0;
    double last_valid = 0.0;

    for (int i = 0; i < peak_cnt && out_cnt < max_results; ++i) {
        if (last_valid != 0.0 && (peaks[i].xpeak - last_valid <= MIN_INTERPEAK_DIST))
            continue;

        double fwhm = match_slopes_and_fwhm(&peaks[i],
                                            rising, rising_cnt,
                                            falling, falling_cnt);
        if (fwhm > 0.0) {
            results[out_cnt].location   = peaks[i].xpeak;
            results[out_cnt].width      = fwhm;
            results[out_cnt].separation = peaks[i].separation;
            ++out_cnt;
            last_valid = peaks[i].xpeak;
        }
    }

    free(x); free(y); free(peaks);
    free(rising); free(falling);
    return out_cnt;
}
