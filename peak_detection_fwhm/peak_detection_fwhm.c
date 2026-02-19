// EE25B094, EE25B126, EE25B138
// USED GEMINI AND GROK AI TO DEBUG AND MAKE THE SYNTAX MORE READABLE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

// Structure to hold x,y data points from input
typedef struct {
    int x;
    double y;
} DataPoint;

// Structure to hold peak data
typedef struct {
    double xpeak;
    double ypeak;
    double separation;
} Peak;

// Structure to hold a single calculated slope
typedef struct {
    double m;             // Slope
    double c;             // Intercept
    int region_start_x;   // The x-value where this slope started
    bool used;            // Track if slope is used
} Slope;

// Comparator for sorting peaks by xpeak ascending
int compare_peaks(const void *a, const void *b) {
    double xa = ((Peak*)a)->xpeak;
    double xb = ((Peak*)b)->xpeak;
    return (xa > xb) - (xa < xb);
}

int main(int argc, char **argv) {
    // Check command-line arguments
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <input_filename>\n", argv[0]);
        return 1;
    }

    // Open input file
    char *name = argv[1];
    FILE *fp = fopen(name, "r");
    if (!fp) {
        perror("Failed to open input file");
        return 1;
    }

    // Constants for peak detection
    const double delta1 = 100.0;      // Min distance to end a region after stagnation
    const double delta2 = 2500.0;     // Min separation to start a new peak region
    const int min_points_for_peak = 50; // Min data points required to fit a peak

    // Allocate dynamic array for data points
    int data_capacity = 1000000; // Adjust based on expected input size
    int data_count = 0;
    DataPoint *data = malloc(data_capacity * sizeof(DataPoint));
    if (!data) {
        perror("Memory allocation failed for data");
        fclose(fp);
        return 1;
    }

    // Read input file into DataPoint array, skip header
    char line[256];
    if (!fgets(line, sizeof(line), fp)) {
        fprintf(stderr, "Empty file or read error\n");
        free(data);
        fclose(fp);
        return 1;
    }

    // Load data points
    while (fgets(line, sizeof(line), fp)) {
        int x;
        double y;
        if (sscanf(line, "%d,%lf", &x, &y) != 2) {
            fprintf(stderr, "Invalid data format at line %d\n", data_count + 2);
            continue;
        }
        if (data_count >= data_capacity) {
            data_capacity *= 2;
            DataPoint *temp = realloc(data, data_capacity * sizeof(DataPoint));
            if (!temp) {
                perror("Reallocation failed for data");
                free(data);
                fclose(fp);
                return 1;
            }
            data = temp;
        }
        data[data_count].x = x;
        data[data_count].y = y;
        data_count++;
    }
    fclose(fp);

    // --- Peak Detection ---
    // Dynamic array for peaks
    int peak_capacity = 250;
    int peak_count = 0;
    Peak *peaks = malloc(peak_capacity * sizeof(Peak));
    if (!peaks) {
        perror("Memory allocation failed for peaks");
        free(data);
        return 1;
    }

    // Variables for peak fitting
    double s[7] = {0}; // Summations: sum x, x^2, x^3, x^4, y, x*y, x^2*y
    int count1 = 0;    // Points in current peak region
    int prev_count = 0; // Previous count to detect stagnation
    int xstart = 0, xend = 0; // Region boundaries
    bool in_region = false; // Flag for peak region
    double xpold = 0.0; // Previous peak x-position
    int prev_x = 0;     // To track x for delta checks

    // Process data points for peak detection
    for (int i = 0; i < data_count; i++) {
        int x = data[i].x;
        double y = data[i].y;

        // Start or continue a high-y region (y > 1.5) for peak fitting
        if (y > 1.5 && (xpold == 0 || (x - xpold > delta2))) {
            if (!in_region) {
                xstart = x;
                in_region = true;
                for (int j = 0; j < 7; j++) s[j] = 0.0; // Reset summations
                count1 = 0;
            }
            double xrel = x - xstart;
            s[0] += xrel;
            s[1] += xrel * xrel;
            s[2] += xrel * xrel * xrel;
            s[3] += xrel * xrel * xrel * xrel;
            s[4] += y;
            s[5] += xrel * y;
            s[6] += xrel * xrel * y;
            count1++;
        }

        // Detect end of region: stagnation and sufficient gap
        if (count1 >= min_points_for_peak && prev_count == count1) {
            xend = prev_x;
        }
        if (xend != 0 && (x - xend >= delta1)) {
            in_region = false;
            // Compute quadratic coefficients a, b, c for y = a*x^2 + b*x + c
            double denom1 = s[1] * s[2] - s[0] * s[3];
            double denom2 = s[1] * s[1] - s[0] * s[2];
            if (fabs(denom1) < 1e-9 || fabs(denom2) < 1e-9) {
                fprintf(stderr, "Singular matrix in peak calculation at x=%d\n", xstart);
                xend = 0;
                prev_count = count1;
                continue;
            }
            double a = (s[1] * s[5] - s[0] * s[6]) / denom1;
            double b = (s[1] * s[5] - s[0] * s[6] - s[1] * s[2] * a + s[0] * s[3] * a) / denom2;
            double c = s[4] / count1 - s[1] * a / count1 - s[0] * b / count1;

            // Calculate peak location and height
            double xpeak = (a != 0) ? -b / (2.0 * a) + xstart : 0;
            double ypeak = (a != 0) ? -(b * b - 4 * a * c) / (4 * a) : 0;

            // Resize peaks array if needed
            if (peak_count >= peak_capacity) {
                peak_capacity *= 2;
                Peak *temp = realloc(peaks, peak_capacity * sizeof(Peak));
                if (!temp) {
                    perror("Reallocation failed for peaks");
                    free(data);
                    free(peaks);
                    return 1;
                }
                peaks = temp;
            }
            peaks[peak_count].xpeak = xpeak;
            peaks[peak_count].ypeak = ypeak;
            peaks[peak_count].separation = (xpold == 0) ? 0 : (xpeak - xpold);
            peak_count++;

            // Reset for next region
            xend = 0;
            xpold = xpeak;
            for (int j = 0; j < 7; j++) s[j] = 0;
            count1 = 0;
        }
        prev_count = count1;
        prev_x = x;
    }

    // Sort peaks by xpeak
    qsort(peaks, peak_count, sizeof(Peak), compare_peaks);

    // --- Slope Calculation and FWHM ---
    // Allocate arrays for rising and falling slopes
    int slope_capacity = 2000;
    int rising_count = 0, falling_count = 0;
    Slope *rising_slopes = malloc(slope_capacity * sizeof(Slope));
    Slope *falling_slopes = malloc(slope_capacity * sizeof(Slope));
    if (!rising_slopes || !falling_slopes) {
        perror("Memory allocation failed for slopes");
        free(data);
        free(peaks);
        free(rising_slopes);
        free(falling_slopes);
        return 1;
    }

    // Variables for slope detection in mid-y range (0.85 < y < 1.4)
    int n = 0;
    double l[4] = {0}; // sum x*y, sum x, sum y, sum x^2
    int region_start_x = 0;
    const double slope_threshold = 0.00005;

    // Process data points for slope calculation
    for (int i = 0; i < data_count; i++) {
        int current_x = data[i].x;
        double y = data[i].y;

        if (y > 0.85 && y < 1.4) {
            if (n == 0) region_start_x = current_x;
            double x_rel = current_x - region_start_x;
            l[0] += x_rel * y;
            l[1] += x_rel;
            l[2] += y;
            l[3] += x_rel * x_rel;
            n++;
        } else if (n > 0) {
            // Calculate slope when exiting valid y-range
            if (n > 5) {
                double denom = n * l[3] - l[1] * l[1];
                if (fabs(denom) > 1e-9) {
                    double m = (n * l[0] - l[1] * l[2]) / denom;
                    double c = (l[2] / n) - m * (l[1] / n);
                    if (fabs(m) > slope_threshold) {
                        if (m > 0 && rising_count < slope_capacity) {
                            rising_slopes[rising_count].m = m;
                            rising_slopes[rising_count].c = c;
                            rising_slopes[rising_count].region_start_x = region_start_x;
                            rising_slopes[rising_count].used = false;
                            rising_count++;
                        } else if (m < 0 && falling_count < slope_capacity) {
                            falling_slopes[falling_count].m = m;
                            falling_slopes[falling_count].c = c;
                            falling_slopes[falling_count].region_start_x = region_start_x;
                            falling_slopes[falling_count].used = false;
                            falling_count++;
                        }
                    }
                } else {
                    fprintf(stderr, "Singular matrix in slope calculation at x=%d\n", region_start_x);
                }
            }
            // Reset for next region
            n = 0;
            for (int j = 0; j < 4; j++) l[j] = 0;
        }
    }

    // Process any remaining slope at end of data
    if (n > 5) {
        double denom = n * l[3] - l[1] * l[1];
        if (fabs(denom) > 1e-9) {
            double m = (n * l[0] - l[1] * l[2]) / denom;
            double c = (l[2] / n) - m * (l[1] / n);
            if (fabs(m) > slope_threshold) {
                if (m > 0 && rising_count < slope_capacity) {
                    rising_slopes[rising_count].m = m;
                    rising_slopes[rising_count].c = c;
                    rising_slopes[rising_count].region_start_x = region_start_x;
                    rising_slopes[rising_count].used = false;
                    rising_count++;
                } else if (m < 0 && falling_count < slope_capacity) {
                    falling_slopes[falling_count].m = m;
                    falling_slopes[falling_count].c = c;
                    falling_slopes[falling_count].region_start_x = region_start_x;
                    falling_slopes[falling_count].used = false;
                    falling_count++;
                }
            }
        } else {
            fprintf(stderr, "Singular matrix in final slope calculation at x=%d\n", region_start_x);
        }
    }

    // --- Match Slopes to Peaks and Compute FWHM ---
    double last_valid_xpeak = 0.0;
    const double min_interpeak_distance = 3000.0;
    const double max_slope_distance = 2500.0;

    for (int i = 0; i < peak_count; i++) {
        if (last_valid_xpeak != 0.0 && (peaks[i].xpeak - last_valid_xpeak) <= min_interpeak_distance) {
            continue;
        }

        // Find closest unused rising slope before peak
        int best_rising_idx = -1;
        double min_dist_x1 = 1e9;
        for (int j = 0; j < rising_count; j++) {
            if (!rising_slopes[j].used && rising_slopes[j].region_start_x < peaks[i].xpeak) {
                double dist = peaks[i].xpeak - rising_slopes[j].region_start_x;
                if (dist < min_dist_x1 && dist <= max_slope_distance) {
                    min_dist_x1 = dist;
                    best_rising_idx = j;
                }
            }
        }

        // Find closest unused falling slope after peak
        int best_falling_idx = -1;
        double min_dist_x2 = 1e9;
        for (int j = 0; j < falling_count; j++) {
            if (!falling_slopes[j].used && falling_slopes[j].region_start_x > peaks[i].xpeak) {
                double dist = falling_slopes[j].region_start_x - peaks[i].xpeak;
                if (dist < min_dist_x2 && dist <= max_slope_distance) {
                    min_dist_x2 = dist;
                    best_falling_idx = j;
                }
            }
        }

        // Compute FWHM if both slopes are found
        if (best_rising_idx != -1 && best_falling_idx != -1) {
            Slope rising = rising_slopes[best_rising_idx];
            Slope falling = falling_slopes[best_falling_idx];

            double y_half = peaks[i].ypeak / 2.0;
            double x1_rel = (y_half - rising.c) / rising.m;
            double x1 = x1_rel + rising.region_start_x;
            double x2_rel = (y_half - falling.c) / falling.m;
            double x2 = x2_rel + falling.region_start_x;
            double fwhm = x2 - x1;

            if (fwhm > 0 && fwhm < 2000) {
                printf("%-15.2f %-15.2f %-15.2f\n", peaks[i].xpeak, fwhm, peaks[i].separation);
                last_valid_xpeak = peaks[i].xpeak;
                rising_slopes[best_rising_idx].used = true;
                falling_slopes[best_falling_idx].used = true;
            }
        }
    }

    // Cleanup
    free(data);
    free(peaks);
    free(rising_slopes);
    free(falling_slopes);
    return 0;
}
