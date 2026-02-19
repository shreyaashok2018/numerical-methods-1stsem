/**
 * @file ee25b126_utility.c
 * @brief Utility Functions (Part 10) – header-required + binary-bit printer
 * @author EE25B126
 * @date November 2025
 */

#include "ee25b126_ee1103.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ============================================================================
   10. REQUIRED HEADER FUNCTIONS
   ============================================================================ */

int read_csv_data(const char *filename, double *x, double *y, int max_points) {
    if (!filename || !x || !y || max_points <= 0) return 0;

    FILE *fp = fopen(filename, "r");
    if (!fp) return 0;

    int n = 0;
    char line[1024];

    /* skip possible header line */
    if (fgets(line, sizeof(line), fp) == NULL) {
        fclose(fp);
        return 0;
    }

    while (n < max_points && fgets(line, sizeof(line), fp)) {
        if (sscanf(line, "%lf,%lf", &x[n], &y[n]) == 2)
            ++n;
    }
    fclose(fp);
    return n;
}

void print_array(const double *arr, int N, const char *label) {
    if (!arr || N <= 0) return;
    if (label) printf("%s: ", label);
    for (int i = 0; i < N; ++i) {
        printf("%.6f", arr[i]);
        if (i < N - 1) putchar(' ');
    }
    putchar('\n');
}

/* ============================================================================
   10. BONUS – BINARY BIT-PRINTER (modular, reusable) & GNUPLOT PLOTTER
   ============================================================================ */

/**
 * @brief Print the first *num_bits* bits of a binary file (MSB first).
 *
 * @param filename   Path to the binary file.
 * @param num_bits   Number of bits to print; -1 means “all bits”.
 * @return 0 on success, non-zero on error.
 *
 * The function prints the bits without any separator and ends with a newline.
 */
int print_binary_bits(const char *filename, long num_bits) {
    if (!filename) return 1;

    FILE *fp = fopen(filename, "rb");
    if (!fp) {
        perror("print_binary_bits: fopen");
        return 1;
    }

    long printed = 0;
    int byte;

    while ((byte = fgetc(fp)) != EOF) {
        for (int i = 7; i >= 0; --i) {               // MSB to LSB
            if (num_bits != -1 && printed >= num_bits)
                goto done;

            putchar( (byte & (1 << i)) ? '1' : '0' );
            ++printed;
        }
    }

done:
    putchar('\n');
    fclose(fp);
    return 0;
}

/**
 * @brief Generic GNUplot plotter for CSV (x,y) data
 */
int plot_csv_data(const char *csv_filename,
                  const char *title,
                  const char *xlabel,
                  const char *ylabel,
                  const char *output_png) {
    if (!csv_filename || !output_png) return 1;

    /* Default labels if NULL */
    const char *t = title ? title : "Data Plot";
    const char *xl = xlabel ? xlabel : "X";
    const char *yl = ylabel ? ylabel : "Y";

    /* Temporary GNUplot script */
    const char *gp_script = ".tmp_plot.gp";
    FILE *fp = fopen(gp_script, "w");
    if (!fp) {
        perror("Failed to create GNUplot script");
        return 1;
    }

    fprintf(fp,
        "set terminal pngcairo enhanced font 'arial,10' size 800,600\n"
        "set output '%s'\n"
        "set title '%s'\n"
        "set xlabel '%s'\n"
        "set ylabel '%s'\n"
        "set grid\n"
        "set key off\n"
        "plot '%s' using 1:2 with lines linewidth 1.5 linecolor rgb '#0066CC'\n"
        "set terminal pop\n"
        "replot\n",
        output_png, t, xl, yl, csv_filename
    );
    fclose(fp);

    /* Run GNUplot */
    char cmd[1024];
    snprintf(cmd, sizeof(cmd), "gnuplot %s", gp_script);
    int ret = system(cmd);

    /* Cleanup */
    remove(gp_script);

    if (ret != 0) {
        fprintf(stderr, "GNUplot failed. Is it installed?\n");
        fprintf(stderr, "Install: sudo apt install gnuplot  (or brew install gnuplot)\n");
        fprintf(stderr, "Data file: %s\n", csv_filename);
        return 1;
    }

    printf("Plot saved: %s\n", output_png);
    return 0;
}
