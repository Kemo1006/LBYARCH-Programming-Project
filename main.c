#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <windows.h>

/* x86-64 NASM kernel (see daxpy.asm) */
extern void daxpy_asm(int n, double a, double *x, double *y, double *z);

#define NUM_RUNS   30
#define EPSILON    1e-9

/* Vector sizes to benchmark: 2^20, 2^24, 2^28 (see note above re: 2^30) */
static const long long N_SIZES[] = { 1LL << 20, 1LL << 24, 1LL << 28 };
static const int NUM_SIZES = (int)(sizeof(N_SIZES) / sizeof(N_SIZES[0]));

/* C reference kernel                                                     */
void daxpy_c(int n, double a, double *x, double *y, double *z)
{
    for (int i = 0; i < n; i++) {
        z[i] = a * x[i] + y[i];
    }
}

/* Helpers                                                                 */
static double *alloc_vec(long long n)
{
    double *v = (double *)malloc(sizeof(double) * (size_t)n);
    if (!v) {
        fprintf(stderr, "Error: malloc failed for n = %lld (out of memory)\n", n);
        exit(EXIT_FAILURE);
    }
    return v;
}

static void fill_random(double *v, long long n)
{
    for (long long i = 0; i < n; i++) {
        v[i] = ((double)rand() / (double)RAND_MAX) * 10.0; /* [0, 10) */
    }
}

/* Returns 1 if all elements match within EPSILON, else 0 */
static int correctness_check(const double *ref, const double *test, long long n)
{
    for (long long i = 0; i < n; i++) {
        if (fabs(ref[i] - test[i]) > EPSILON) {
            printf("  Mismatch at index %lld: C=%.10f  ASM=%.10f\n",
                   i, ref[i], test[i]);
            return 0;
        }
    }
    return 1;
}

static void print_first10(const char *label, const double *z, long long n)
{
    long long lim = (n < 10) ? n : 10;
    printf("%s first %lld elements: ", label, lim);
    for (long long i = 0; i < lim; i++) {
        printf("%.4f ", z[i]);
    }
    printf("\n");
}

/* High resolution timer (Windows) */
static double now_ms(void)
{
    static LARGE_INTEGER freq;
    static int have_freq = 0;
    LARGE_INTEGER t;
    if (!have_freq) {
        QueryPerformanceFrequency(&freq);
        have_freq = 1;
    }
    QueryPerformanceCounter(&t);
    return (double)t.QuadPart * 1000.0 / (double)freq.QuadPart;
}

/* Benchmark driver for one kernel */
typedef void (*daxpy_fn)(int, double, double *, double *, double *);

/* Summary statistics for one kernel's set of timed runs.
 *
 * raw_mean is the simple average of all runs. On large workloads (e.g.
 * n = 2^28, where several GB of memory are live at once) occasional runs
 * can be hit by OS-level effects unrelated to the kernel itself (page
 * faults, working-set trimming, background paging) causing rare but huge
 * outlier times. A plain mean is heavily skewed by these; median and
 * trimmed_mean are more robust and are the recommended headline figures. */
typedef struct {
    double min;
    double median;
    double max;
    double raw_mean;      /* simple average of all runs (sensitive to outliers) */
    double trimmed_mean;  /* average after dropping the top/bottom TRIM_COUNT runs */
} KernelTiming;

#define TRIM_COUNT 3  /* drop this many of the lowest and highest runs before averaging */

/* Runs the kernel `runs` times, times each call individually, prints the
 * full distribution, and returns summary statistics (min/median/max/mean/
 * trimmed mean). */
static KernelTiming time_kernel(const char *label, daxpy_fn fn, int n, double a,
                                 double *x, double *y, double *z, int runs)
{
    double *times = (double *)malloc(sizeof(double) * runs);
    double total_ms = 0.0;
    for (int r = 0; r < runs; r++) {
        double t0 = now_ms();
        fn(n, a, x, y, z);
        double t1 = now_ms();
        times[r] = t1 - t0;
        total_ms += times[r];
    }

    /* simple insertion sort for median/min/max/trimmed mean, runs is small (30) */
    for (int i = 1; i < runs; i++) {
        double key = times[i];
        int j = i - 1;
        while (j >= 0 && times[j] > key) { times[j + 1] = times[j]; j--; }
        times[j + 1] = key;
    }

    printf("  [%s] per-run times (ms): ", label);
    for (int r = 0; r < runs; r++) printf("%.1f ", times[r]);
    printf("\n");

    KernelTiming t;
    t.min      = times[0];
    t.median   = (runs % 2 == 0) ? (times[runs / 2 - 1] + times[runs / 2]) / 2.0
                                  : times[runs / 2];
    t.max      = times[runs - 1];
    t.raw_mean = total_ms / runs;

    int trim = (runs > 2 * TRIM_COUNT) ? TRIM_COUNT : 0; /* skip trimming if too few runs */
    double trimmed_sum = 0.0;
    for (int r = trim; r < runs - trim; r++) trimmed_sum += times[r];
    t.trimmed_mean = trimmed_sum / (runs - 2 * trim);

    printf("  [%s] min=%.2f  median=%.2f  max=%.2f  raw_mean=%.2f  trimmed_mean(drop %d each end)=%.2f\n",
           label, t.min, t.median, t.max, t.raw_mean, trim, t.trimmed_mean);

    free(times);
    return t;
}

/* Main */
int main(void)
{
    srand(12345); /* fixed seed for reproducibility */

    printf("==================================================================\n");
    printf(" DAXPY Benchmark:  Z[i] = A * X[i] + Y[i]\n");
    printf(" Kernels: C (reference)  vs  x86-64 NASM (scalar SSE2)\n");
    printf(" Runs per size: %d (averaged)\n", NUM_RUNS);
    printf("==================================================================\n\n");

    for (int s = 0; s < NUM_SIZES; s++) {
        long long n_ll = N_SIZES[s];
        if (n_ll > 2147483647LL) {
            fprintf(stderr, "Skipping n=%lld: exceeds int range for kernel signature\n", n_ll);
            continue;
        }
        int n = (int)n_ll;
        double a = ((double)rand() / (double)RAND_MAX) * 10.0;

        printf("------------------------------------------------------------\n");
        printf("n = %d (2^%d)\n", n, (int)(log((double)n) / log(2.0) + 0.5));
        printf("------------------------------------------------------------\n");

        double *x     = alloc_vec(n);
        double *y     = alloc_vec(n);
        double *z_c   = alloc_vec(n);
        double *z_asm = alloc_vec(n);

        fill_random(x, n);
        fill_random(y, n);

        /* C kernel */
        printf("[C kernel]\n");
        printf("  A = %.6f\n", a);
        KernelTiming t_c = time_kernel("C", daxpy_c, n, a, x, y, z_c, NUM_RUNS);
        print_first10("  Z", z_c, n);
        /* Average is the headline figure required by the spec. Median and
         * trimmed mean are also reported (see per-run line above) as
         * supplementary detail -- useful context when a run is hit by an
         * unrelated OS-level pause (e.g. paging at very large n), since a
         * single such outlier can pull the plain average up noticeably. */
        printf("  Average execution time over %d runs: %.6f ms  (median: %.6f ms, trimmed mean: %.6f ms)\n",
               NUM_RUNS, t_c.raw_mean, t_c.median, t_c.trimmed_mean);

        /* x86-64 ASM kernel */
        printf("[x86-64 ASM kernel]\n");
        KernelTiming t_asm = time_kernel("ASM", daxpy_asm, n, a, x, y, z_asm, NUM_RUNS);
        print_first10("  Z", z_asm, n);
        printf("  Average execution time over %d runs: %.6f ms  (median: %.6f ms, trimmed mean: %.6f ms)\n",
               NUM_RUNS, t_asm.raw_mean, t_asm.median, t_asm.trimmed_mean);

        /* correctness check: ASM output vs C */
        int ok = correctness_check(z_c, z_asm, n);
        printf("  Correctness check (ASM vs C): %s\n", ok ? "PASSED" : "FAILED");

        /* Speedup computed from the average, as required by the spec. */
        printf("  Speedup (C average / ASM average): %.3fx\n\n", t_c.raw_mean / t_asm.raw_mean);

        free(x);
        free(y);
        free(z_c);
        free(z_asm);
    }

    return 0;
}