#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <errno.h>
#include <windows.h>

// External assembly function implementation (NASM)
// Performs DAXPY using SSE2 instructions for vectorization
extern void daxpy_asm(int n, double a, double *x, double *y, double *z);

// Configuration constants
#define NUM_RUNS   30              // Number of timing runs per kernel for statistical reliability
#define EPSILON    1e-9            // Tolerance for floating-point comparison in correctness check
#define MANUAL_ENTRY_MAX 20        // Maximum vector size for manual data entry (prevents buffer overflow)

// Default test vector sizes: 1M, 16M, 268M elements
static const long long N_SIZES[] = { 1LL << 20, 1LL << 24, 1LL << 28 };
static const int NUM_SIZES = (int)(sizeof(N_SIZES) / sizeof(N_SIZES[0]));

/**
 * C reference implementation of DAXPY
 * Simple loop performing the element-wise operation
 *
 * @param n     Number of elements in vectors
 * @param a     Scalar multiplier
 * @param x     Input vector X (size n)
 * @param y     Input vector Y (size n)
 * @param z     Output vector Z (size n) - result of Z[i] = a * X[i] + Y[i]
 */
void daxpy_c(int n, double a, double *x, double *y, double *z) {
    for (int i = 0; i < n; i++) {
        z[i] = a * x[i] + y[i];
    }
}

/**
 * Allocate memory for a double-precision vector
 * Exits with error message if allocation fails
 *
 * @param n     Number of double elements to allocate
 * @return      Pointer to allocated memory block
 */
static double *alloc_vec(long long n) {
    double *v = (double *)malloc(sizeof(double) * (size_t)n);
    if (!v) {
        fprintf(stderr, "Error: malloc failed for n = %lld (out of memory)\n", n);
        exit(EXIT_FAILURE);
    }
    return v;
}

/**
 * Fill a vector with random values in range [0, 10)
 * Uses rand() for simplicity (not cryptographic quality needed)
 *
 * @param v     Pointer to vector to fill
 * @param n     Number of elements to fill
 */
static void fill_random(double *v, long long n) {
    for (long long i = 0; i < n; i++) {
        v[i] = ((double)rand() / (double)RAND_MAX) * 10.0; 
    }
}

/**
 * Compare two vectors element-wise with tolerance
 * Returns 1 if all elements match within EPSILON, 0 otherwise
 * Prints first mismatch for debugging
 *
 * @param ref   Reference vector (expected values)
 * @param test  Test vector (values to verify)
 * @param n     Number of elements to compare
 * @return      1 if all elements match within tolerance, 0 otherwise
 */
static int correctness_check(const double *ref, const double *test, long long n) {
    for (long long i = 0; i < n; i++) {
        if (fabs(ref[i] - test[i]) > EPSILON) {
            printf("  Mismatch at index %lld: C=%.10f  ASM=%.10f\n",
                   i, ref[i], test[i]);
            return 0;
        }
    }
    return 1;
}

/**
 * Utility function to print first 10 elements (or fewer for small vectors)
 * Useful for quick verification of results
 *
 * @param label Descriptive label for output
 * @param z     Vector to print
 * @param n     Total number of elements in vector
 */
static void print_first10(const char *label, const double *z, long long n) {
    long long lim = (n < 10) ? n : 10;
    printf("%s first %lld elements: ", label, lim);
    for (long long i = 0; i < lim; i++) {
        printf("%.4f ", z[i]);
    }
    printf("\n");
}

/**
 * High-resolution timestamp function using Windows Performance Counter
 * Returns current time in milliseconds with microsecond resolution
 *
 * @return      Current time in milliseconds since system start
 */
static double now_ms(void) {
    static LARGE_INTEGER freq;
    static int have_freq = 0;
    LARGE_INTEGER t;
    if (!have_freq) {
        QueryPerformanceFrequency(&freq);  // Get counter frequency once
        have_freq = 1;
    }
    QueryPerformanceCounter(&t);
    return (double)t.QuadPart * 1000.0 / (double)freq.QuadPart;
}

// Function pointer type for DAXPY kernels (supports both C and ASM)
typedef void (*daxpy_fn)(int, double, double *, double *, double *);

// Structure to hold comprehensive timing statistics
typedef struct {
    double min;          // Minimum execution time
    double median;       // Median execution time (robust to outliers)
    double max;          // Maximum execution time
    double raw_mean;     // Arithmetic mean of all runs
    double trimmed_mean; // Mean after removing extremes (more robust)
} KernelTiming;

#define TRIM_COUNT 3   // Number of extreme values to trim from each end

/**
 * Time a kernel function over multiple runs and compute statistics
 * Sorts times for median and trimmed mean calculations
 * Returns comprehensive timing statistics structure
 *
 * @param label     Descriptive label for output (e.g., "C" or "ASM")
 * @param fn        Function pointer to the kernel to time
 * @param n         Number of elements in vectors
 * @param a         Scalar multiplier
 * @param x         Input vector X
 * @param y         Input vector Y
 * @param z         Output vector Z (will be overwritten)
 * @param runs      Number of timing runs to perform
 * @return          KernelTiming structure containing timing statistics
 */
static KernelTiming time_kernel(const char *label, daxpy_fn fn, int n, double a,
                                 double *x, double *y, double *z, int runs) {
    double *times = (double *)malloc(sizeof(double) * runs);
    double total_ms = 0.0;
    
    // Perform all timing runs
    for (int r = 0; r < runs; r++) {
        double t0 = now_ms();
        fn(n, a, x, y, z);      // Execute the kernel
        double t1 = now_ms();
        times[r] = t1 - t0;
        total_ms += times[r];
    }

    // Sort times using insertion sort (efficient for small arrays)
    for (int i = 1; i < runs; i++) {
        double key = times[i];
        int j = i - 1;
        while (j >= 0 && times[j] > key) { times[j + 1] = times[j]; j--; }
        times[j + 1] = key;
    }

    // Print all run times for transparency
    printf("  [%s] per-run times (ms): ", label);
    for (int r = 0; r < runs; r++) printf("%.1f ", times[r]);
    printf("\n");

    // Compute statistics
    KernelTiming t;
    t.min      = times[0];
    t.median   = (runs % 2 == 0) ? (times[runs / 2 - 1] + times[runs / 2]) / 2.0
                                  : times[runs / 2];
    t.max      = times[runs - 1];
    t.raw_mean = total_ms / runs;

    // Trimmed mean: remove 3 smallest and 3 largest values
    int trim = (runs > 2 * TRIM_COUNT) ? TRIM_COUNT : 0; 
    double trimmed_sum = 0.0;
    for (int r = trim; r < runs - trim; r++) trimmed_sum += times[r];
    t.trimmed_mean = trimmed_sum / (runs - 2 * trim);

    // Print summary statistics
    printf("  [%s] min=%.2f  median=%.2f  max=%.2f  raw_mean=%.2f  trimmed_mean(drop %d each end)=%.2f\n",
           label, t.min, t.median, t.max, t.raw_mean, trim, t.trimmed_mean);

    free(times);
    return t;
}

/**
 * Run benchmark for a specific vector size
 * Allocates vectors, times both kernels, checks correctness
 * Optionally uses user-provided data instead of random values
 *
 * @param n_ll          Vector size in long long (will be cast to int)
 * @param has_fixed_a   Boolean flag (0=random A, 1=use fixed_a)
 * @param fixed_a       Fixed scalar value to use if has_fixed_a is 1
 * @param manual_x      Pre-filled X vector (NULL for random)
 * @param manual_y      Pre-filled Y vector (NULL for random)
 * @return              1 if correctness check passes, 0 if it fails
 */
static int run_one_size(long long n_ll, int has_fixed_a, double fixed_a,
                         const double *manual_x, const double *manual_y) {
    // Validate: ensure size is positive and fits in int (for kernel signature)
    if (n_ll <= 0) {
        fprintf(stderr, "Skipping n=%lld: vector size must be positive\n", n_ll);
        return 0;
    }
    if (n_ll > 2147483647LL) {
        fprintf(stderr, "Skipping n=%lld: exceeds int range for kernel signature\n", n_ll);
        return 0;
    }
    int n = (int)n_ll;
    
    // Use fixed A if provided, otherwise random
    double a = has_fixed_a ? fixed_a : ((double)rand() / (double)RAND_MAX) * 10.0;

    printf("------------------------------------------------------------\n");
    printf("n = %d (2^%.2f)\n", n, log((double)n) / log(2.0));
    printf("------------------------------------------------------------\n");

    // Allocate input and output vectors
    double *x     = alloc_vec(n);
    double *y     = alloc_vec(n);
    double *z_c   = alloc_vec(n);   // Output from C kernel
    double *z_asm = alloc_vec(n);   // Output from ASM kernel

    // Initialize X and Y (use manual values if provided, otherwise random)
    if (manual_x != NULL) {
        memcpy(x, manual_x, sizeof(double) * (size_t)n);
    } else {
        fill_random(x, n);
    }
    if (manual_y != NULL) {
        memcpy(y, manual_y, sizeof(double) * (size_t)n);
    } else {
        fill_random(y, n);
    }

    // Time the C reference implementation
    printf("[C kernel]\n");
    printf("  A = %.6f\n", a);
    KernelTiming t_c = time_kernel("C", daxpy_c, n, a, x, y, z_c, NUM_RUNS);
    print_first10("  Z", z_c, n);
    printf("  Average execution time over %d runs: %.6f ms  (median: %.6f ms, trimmed mean: %.6f ms)\n",
           NUM_RUNS, t_c.raw_mean, t_c.median, t_c.trimmed_mean);

    // Time the assembly implementation
    printf("[x86-64 ASM kernel]\n");
    KernelTiming t_asm = time_kernel("ASM", daxpy_asm, n, a, x, y, z_asm, NUM_RUNS);
    print_first10("  Z", z_asm, n);
    printf("  Average execution time over %d runs: %.6f ms  (median: %.6f ms, trimmed mean: %.6f ms)\n",
           NUM_RUNS, t_asm.raw_mean, t_asm.median, t_asm.trimmed_mean);

    // Verify correctness by comparing with reference implementation
    int ok = correctness_check(z_c, z_asm, n);
    printf("  Correctness check (ASM vs C): %s\n", ok ? "PASSED" : "FAILED");

    // Report speedup factor (C average / ASM average)
    printf("  Speedup (C average / ASM average): %.3fx\n\n", t_c.raw_mean / t_asm.raw_mean);

    // Free allocated memory
    free(x);
    free(y);
    free(z_c);
    free(z_asm);

    return ok;
}

/**
 * Parse a line of space-separated numbers into a double array
 * Returns 1 on success, 0 on failure
 *
 * @param line      Input string containing space-separated numbers
 * @param out       Output array to fill with parsed values
 * @param n         Expected number of values
 * @return          1 if exactly n valid numbers parsed, 0 otherwise
 */
static int parse_vector_line(char *line, double *out, int n) {
    int count = 0;
    char *tok = strtok(line, " \t");
    while (tok != NULL) {
        if (count >= n) return 0;  // Too many values
        char *end;
        double val = strtod(tok, &end);
        if (end == tok || *end != '\0') return 0;  // Invalid number
        if (isnan(val) || isinf(val)) return 0;   // NaN or infinity
        out[count++] = val;
        tok = strtok(NULL, " \t");
    }
    return (count == n) ? 1 : 0;  // Must have exactly n values
}

/**
 * Parse a vector size that can be either:
 * - A plain integer (e.g., 1048576)
 * - A power of two shorthand (e.g., 2^20)
 * Note: On Windows cmd, 2^20 must be quoted as "2^20" because ^ is escape char
 *
 * @param s     Input string to parse
 * @return      Parsed vector size if valid, -1 on error
 */
static long long parse_size_arg(const char *s) {
    // Skip leading whitespace
    while (*s == ' ' || *s == '\t') s++;
    if (*s == '\0') return -1;

    // Check for power-of-two notation (e.g., "2^20")
    char *caret = strchr(s, '^');
    if (caret != NULL) {
        // Parse base and exponent
        char base_str[32];
        size_t base_len = (size_t)(caret - s);
        if (base_len == 0 || base_len >= sizeof(base_str)) return -1;
        memcpy(base_str, s, base_len);
        base_str[base_len] = '\0';

        char *base_end, *exp_end;
        long long base = strtoll(base_str, &base_end, 10);
        long long exp  = strtoll(caret + 1, &exp_end, 10);
        if (*base_end != '\0' || base_end == base_str) return -1;
        if (*exp_end != '\0' || exp_end == caret + 1) return -1;
        if (exp < 0 || exp > 62 || base <= 0) return -1;

        // Compute base^exp with overflow checking
        long long result = 1;
        for (long long i = 0; i < exp; i++) {
            if (result > (2147483647LL / (base == 0 ? 1 : base))) return -1; 
            result *= base;
        }
        return result;
    }

    // Parse plain integer with error checking
    char *end;
    errno = 0;
    long long val = strtoll(s, &end, 10);
    if (end == s || *end != '\0') return -1;   
    if (errno == ERANGE) return -1;            // Overflow/underflow
    return val;
}

/**
 * Parse a double value from string, checking for validity
 *
 * @param s     Input string to parse
 * @param out   Pointer to store parsed value
 * @return      1 if valid finite double parsed, 0 otherwise
 */
static int parse_double_arg(const char *s, double *out) {
    while (*s == ' ' || *s == '\t') s++;
    if (*s == '\0') return 0;
    char *end;
    double val = strtod(s, &end);
    if (end == s || *end != '\0') return 0;  // No conversion or trailing chars
    if (isnan(val) || isinf(val)) return 0;  // NaN or infinity
    *out = val;
    return 1;
}

/**
 * Main program entry point
 * Supports three modes:
 * 1. Command-line: ./program [size] [A] - run single test with optional A
 * 2. Interactive: run with user input for size, A, and manual vectors
 * 3. Default: run the standard sweep (2^20, 2^24, 2^28)
 *
 * @param argc      Number of command-line arguments
 * @param argv      Array of command-line argument strings
 * @return          EXIT_SUCCESS on success, EXIT_FAILURE on error
 */
int main(int argc, char *argv[]) {
    // Seed random number generator
    srand((unsigned int)time(NULL));

    // Print program banner
    printf("==================================================================\n");
    printf(" DAXPY Benchmark:  Z[i] = A * X[i] + Y[i]\n");
    printf(" Kernels: C (reference)  vs  x86-64 NASM (scalar SSE2)\n");
    printf(" Runs per size: %d (averaged)\n", NUM_RUNS);
    printf("==================================================================\n\n");

    // Validate command line arguments
    if (argc > 3) {
        fprintf(stderr, "Usage: %s [vector_size] [A]\n", argv[0]);
        fprintf(stderr, "  vector_size accepts a plain integer (e.g. 1048576)\n");
        fprintf(stderr, "  or a power-of-two shorthand (e.g. 2^20 -- NOTE: on\n");
        fprintf(stderr, "  Windows cmd.exe, quote it: \"2^20\", since ^ is cmd's\n");
        fprintf(stderr, "  own escape character and will otherwise mangle it).\n");
        fprintf(stderr, "  A is an optional scalar; omit it for a random A.\n");
        fprintf(stderr, "  Omit both to run the default sweep: n = 2^20, 2^24, 2^28.\n");
        return EXIT_FAILURE;
    }

    // Mode 1: Command-line specified size
    if (argc >= 2) {
        long long n = parse_size_arg(argv[1]);
        if (n <= 0) {
            fprintf(stderr, "Error: '%s' is not a valid positive vector size.\n", argv[1]);
            return EXIT_FAILURE;
        }
        double a = 0.0;
        int has_a = 0;
        if (argc == 3) {
            if (!parse_double_arg(argv[2], &a)) {
                fprintf(stderr, "Error: '%s' is not a valid finite scalar A.\n", argv[2]);
                return EXIT_FAILURE;
            }
            has_a = 1;
        }
        int ran = run_one_size(n, has_a, a, NULL, NULL); 
        return ran ? 0 : EXIT_FAILURE;
    }

    // Mode 2: Interactive mode
    printf("Enter vector size n (plain integer or e.g. 2^26),\n");
    printf("or press Enter to run the default sweep (2^20, 2^24, 2^28): ");
    fflush(stdout);

    char line[256];
    if (fgets(line, sizeof(line), stdin) != NULL) {
        // Clean up input (remove newline and trailing whitespace)
        size_t len = strlen(line);
        while (len > 0 && (line[len - 1] == '\n' || line[len - 1] == '\r' ||
                            line[len - 1] == ' '  || line[len - 1] == '\t')) {
            line[--len] = '\0';
        }

        if (len > 0) {
            long long n = parse_size_arg(line);
            if (n <= 0) {
                fprintf(stderr, "Error: '%s' is not a valid positive vector size.\n", line);
                return EXIT_FAILURE;
            }

            // Ask for scalar A
            printf("Enter scalar A, or press Enter for a random A: ");
            fflush(stdout);
            char aline[256];
            double a = 0.0;
            int has_a = 0;
            if (fgets(aline, sizeof(aline), stdin) != NULL) {
                size_t alen = strlen(aline);
                while (alen > 0 && (aline[alen - 1] == '\n' || aline[alen - 1] == '\r' ||
                                     aline[alen - 1] == ' '  || aline[alen - 1] == '\t')) {
                    aline[--alen] = '\0';
                }
                if (alen > 0) {
                    if (!parse_double_arg(aline, &a)) {
                        fprintf(stderr, "Error: '%s' is not a valid finite scalar A.\n", aline);
                        return EXIT_FAILURE;
                    }
                    has_a = 1;
                }
            }

            // For small vectors, allow manual data entry
            double *manual_x = NULL, *manual_y = NULL;
            if (n <= MANUAL_ENTRY_MAX) {
                manual_x = (double *)malloc(sizeof(double) * (size_t)n);
                manual_y = (double *)malloc(sizeof(double) * (size_t)n);

                printf("Enter %lld values for X, separated by spaces, or press Enter to randomize: ", n);
                fflush(stdout);
                char xline[2048];
                int got_x = 0;
                if (fgets(xline, sizeof(xline), stdin) != NULL) {
                    xline[strcspn(xline, "\r\n")] = '\0';
                    if (strlen(xline) > 0) {
                        if (!parse_vector_line(xline, manual_x, (int)n)) {
                            fprintf(stderr, "Error: expected exactly %lld valid numbers for X.\n", n);
                            return EXIT_FAILURE;
                        }
                        got_x = 1;
                    }
                }

                printf("Enter %lld values for Y, separated by spaces, or press Enter to randomize: ", n);
                fflush(stdout);
                char yline[2048];
                int got_y = 0;
                if (fgets(yline, sizeof(yline), stdin) != NULL) {
                    yline[strcspn(yline, "\r\n")] = '\0';
                    if (strlen(yline) > 0) {
                        if (!parse_vector_line(yline, manual_y, (int)n)) {
                            fprintf(stderr, "Error: expected exactly %lld valid numbers for Y.\n", n);
                            return EXIT_FAILURE;
                        }
                        got_y = 1;
                    }
                }

                // Free unused manual arrays if not provided
                if (!got_x) { free(manual_x); manual_x = NULL; }
                if (!got_y) { free(manual_y); manual_y = NULL; }
            }

            printf("\n");
            int ran = run_one_size(n, has_a, a, manual_x, manual_y);
            free(manual_x);
            free(manual_y);
            return ran ? 0 : EXIT_FAILURE;
        }
    }

    // Mode 3: Default sweep (no input provided)
    printf("\n(no size entered -- running default sweep: n = 2^20, 2^24, 2^28)\n\n");
    int all_ok = 1;
    for (int s = 0; s < NUM_SIZES; s++) {
        if (!run_one_size(N_SIZES[s], 0, 0.0, NULL, NULL)) all_ok = 0;
    }

    return all_ok ? 0 : EXIT_FAILURE;
}