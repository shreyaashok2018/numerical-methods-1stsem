/*
 * ee25b126_matrix.c
 * -----------------------------------------------------------------------------
 * Matrix operations and cubic spline interpolation for EE1103 library.
 *
 * Implements the following functions (as declared in ee25b126_ee1103.h):
 *  - alloc_matrix()
 *  - free_matrix()
 *  - swap_matrix_rows()
 *  - read_matrix_from_file()
 *  - print_solution()
 *  - gaussian_elimination()
 *  - lu_decomposition()
 *  - cubic_spline_interpolate()
 *
 * The cubic spline solver constructs the natural spline system and solves it
 * using the LU decomposition routine for numerical stability.
 * -----------------------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ee25b126_ee1103.h"

/* ============================================================================
   1. BASIC MATRIX UTILITIES
   ============================================================================ */

double **alloc_matrix(int rows, int cols)
{
    double **matrix = (double **)malloc(rows * sizeof(double *));
    if (!matrix) { perror("alloc_matrix: malloc rows"); return NULL; }
    if (rows <= 0 || cols <= 0) { fprintf(stderr, "alloc_matrix: invalid dimensions %d x %d\n", rows, cols); free(matrix); return NULL; }
    if (!matrix) {
        fprintf(stderr, "Memory allocation failed for matrix rows.\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < rows; ++i) {
        matrix[i] = (double *)calloc(cols, sizeof(double));
        if (!matrix[i]) {
            fprintf(stderr, "Memory allocation failed for matrix row %d.\n", i);
            for (int j = 0; j < i; ++j)
                free(matrix[j]);
            free(matrix);
            exit(EXIT_FAILURE);
        }
    }
    return matrix;
}

void free_matrix(double **matrix, int rows)
{
    if (!matrix) return;
    for (int i = 0; i < rows; ++i)
        free(matrix[i]);
    free(matrix);
}

void swap_matrix_rows(double **A, int row1, int row2, int cols)
{
    if (row1 == row2) return;
    for (int j = 0; j < cols; ++j) {
        double temp = A[row1][j];
        A[row1][j] = A[row2][j];
        A[row2][j] = temp;
    }
}

double **read_matrix_from_file(const char *filename, int N)
{
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Error: could not open file %s\n", filename);
        exit(EXIT_FAILURE);
    }

    double **A = alloc_matrix(N, N + 1);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N + 1; ++j) {
            if (fscanf(fp, "%lf", &A[i][j]) != 1) {
                fprintf(stderr, "Error reading matrix from file.\n");
                free_matrix(A, N);
                fclose(fp);
                exit(EXIT_FAILURE);
            }
        }
    }
    fclose(fp);
    return A;
}

void print_solution(const double *x, int N)
{
    printf("Solution Vector:\n");
    for (int i = 0; i < N; ++i)
        printf("x[%d] = %lf\n", i, x[i]);
}

/* ============================================================================
   2. GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
   ============================================================================ */

void gaussian_elimination(double **A, double *x, int N)
{
    for (int i = 0; i < N; ++i) {
        /* Pivoting */
        int maxRow = i;
        for (int k = i + 1; k < N; ++k) {
            if (fabs(A[k][i]) > fabs(A[maxRow][i]))
                maxRow = k;
        }
        swap_matrix_rows(A, i, maxRow, N + 1);

        if (fabs(A[i][i]) < 1e-12) {
            fprintf(stderr, "Error: Singular matrix in Gaussian elimination.\n");
            exit(EXIT_FAILURE);
        }

        /* Elimination */
        for (int k = i + 1; k < N; ++k) {
            double factor = A[k][i] / A[i][i];
            for (int j = i; j <= N; ++j)
                A[k][j] -= factor * A[i][j];
        }
    }

    /* Back substitution */
    for (int i = N - 1; i >= 0; --i) {
        double sum = A[i][N];
        for (int j = i + 1; j < N; ++j)
            sum -= A[i][j] * x[j];
        x[i] = sum / A[i][i];
    }
}

void lu_decomposition(double **A, double *x, int N)
{
    /* Separate coefficient matrix and RHS vector */
    double **U = alloc_matrix(N, N);
    double **L = alloc_matrix(N, N);
    double *b  = (double *)malloc(N * sizeof(double));
    if (!b) { perror("gaussian_elimination: malloc b"); return; }

    for (int i = 0; i < N; ++i) {
        b[i] = A[i][N];
        for (int j = 0; j < N; ++j) {
            U[i][j] = A[i][j];
            L[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    /* LU decomposition with partial pivoting */
    for (int k = 0; k < N; ++k) {
        int maxRow = k;
        for (int i = k + 1; i < N; ++i) {
            if (fabs(U[i][k]) > fabs(U[maxRow][k]))
                maxRow = i;
        }

        /* Swap rows in U, L (only first k cols), and b */
        swap_matrix_rows(U, k, maxRow, N);
        swap_matrix_rows(L, k, maxRow, k);
        double temp = b[k];
        b[k] = b[maxRow];
        b[maxRow] = temp;

        if (fabs(U[k][k]) < 1e-12) {
            fprintf(stderr, "Error: Singular matrix in LU decomposition.\n");
            free_matrix(L, N);
            free_matrix(U, N);
            free(b);
            exit(EXIT_FAILURE);
        }

        for (int i = k + 1; i < N; ++i) {
            double factor = U[i][k] / U[k][k];
            L[i][k] = factor;
            for (int j = k; j < N; ++j)
                U[i][j] -= factor * U[k][j];
        }
    }

    /* Forward substitution: L*y = b */
    double *y = (double *)malloc(N * sizeof(double));
    if (!y) { perror("gaussian_elimination: malloc y"); free(b); return; }
    for (int i = 0; i < N; ++i) {
        y[i] = b[i];
        for (int j = 0; j < i; ++j)
            y[i] -= L[i][j] * y[j];
    }

    /* Back substitution: U*x = y */
    for (int i = N - 1; i >= 0; --i) {
        x[i] = y[i];
        for (int j = i + 1; j < N; ++j)
            x[i] -= U[i][j] * x[j];
        x[i] /= U[i][i];
    }

    free(b);
    free(y);
    free_matrix(L, N);
    free_matrix(U, N);
}

/* ============================================================================
   4. NATURAL CUBIC SPLINE INTERPOLATION
   ============================================================================ */

double cubic_spline_interpolate(const double *x, const double *y, int M, double xp)
{
    if (M < 3) {
        fprintf(stderr, "Error: Need at least 3 points for cubic spline.\n");
        exit(EXIT_FAILURE);
    }

    /* Step 1: Compute intervals and slopes */
    double *h = (double *)malloc((M - 1) * sizeof(double));
    if (!h) { perror("cubic_spline_interpolate: malloc h"); return NAN; }
    double *alpha = (double *)calloc(M, sizeof(double));
    for (int i = 0; i < M - 1; ++i)
        h[i] = x[i + 1] - x[i];

    for (int i = 1; i < M - 1; ++i)
        alpha[i] = (3.0 / h[i]) * (y[i + 1] - y[i]) -
                   (3.0 / h[i - 1]) * (y[i] - y[i - 1]);

    /* Step 2: Build tridiagonal system for c (second derivatives) */
    double **A = alloc_matrix(M, M + 1);
    for (int i = 1; i < M - 1; ++i) {
        A[i][i - 1] = h[i - 1];
        A[i][i]     = 2.0 * (h[i - 1] + h[i]);
        A[i][i + 1] = h[i];
        A[i][M]     = alpha[i];
    }

    /* Natural boundary conditions */
    A[0][0]     = 1.0;
    A[0][M]     = 0.0;
    A[M - 1][M - 1] = 1.0;
    A[M - 1][M]     = 0.0;

    /* Step 3: Solve for c using LU decomposition */
    double *c = (double *)calloc(M, sizeof(double));
    lu_decomposition(A, c, M);

    /* Step 4: Compute b and d coefficients */
    double *b = (double *)malloc((M - 1) * sizeof(double));
    if (!b) { perror("cubic_spline_interpolate: malloc b2"); free(h); free(alpha); free(c); return NAN; }
    double *d = (double *)malloc((M - 1) * sizeof(double));
    if (!d) { perror("cubic_spline_interpolate: malloc d"); free(h); free(alpha); free(c); free(b); return NAN; }
    for (int i = 0; i < M - 1; ++i) {
        b[i] = (y[i + 1] - y[i]) / h[i] - h[i] * (2.0 * c[i] + c[i + 1]) / 3.0;
        d[i] = (c[i + 1] - c[i]) / (3.0 * h[i]);
    }

    /* Step 5: Find the interval where xp lies */
    int i = 0;
    if (xp <= x[0]) i = 0;
    else if (xp >= x[M - 1]) i = M - 2;
    else {
        for (i = 0; i < M - 1; ++i) {
            if (xp >= x[i] && xp <= x[i + 1])
                break;
        }
    }

    double dx = xp - x[i];
    double yp = y[i] + b[i] * dx + c[i] * dx * dx + d[i] * dx * dx * dx;

    /* Clean up */
    free(h);
    free(alpha);
    free_matrix(A, M);
    free(c);
    free(b);
    free(d);

    return yp;
}

