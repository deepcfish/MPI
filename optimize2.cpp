#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define MIN(a,b) ((a)<(b)?(a):(b))
#define LL long long

int main(int argc, char* argv[])
{
    LL count;              // Local prime count
    double elapsed_time;   // Parallel execution time
    LL first;              // Index of first multiple
    LL global_count = 0;   // Global prime count
    LL high_value;         // Highest value on this proc
    LL i;
    int id;                // Process ID number
    LL low_value;          // Lowest value on this proc
    char* marked;          // Portion of 2,...,'n'
    LL n;                  // Sieving from 2, ..., 'n'
    int p;                 // Number of processes
    LL proc0_size;         // Size of proc 0's subarray
    LL size;               // Elements in 'marked'

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime();

    if (argc != 2) {
        if (!id) printf("Command line: %s <m>\n", argv[0]);
        MPI_Finalize();
        exit(1);
    }

    n = atoll(argv[1]);

    low_value = 2 + id * (n - 1) / p;
    high_value = 1 + (id + 1) * (n - 1) / p;
    if (low_value % 2 == 0) low_value++;
    if (high_value % 2 == 0) high_value--;
    size = (high_value - low_value) / 2 + 1;

    proc0_size = (n - 1) / p;
    if ((2 + proc0_size) < (int)sqrt((double)n)) {
        if (!id) printf("Too many processes\n");
        MPI_Finalize();
        exit(1);
    }

    marked = (char*)malloc(size);
    if (marked == NULL) {
        printf("Cannot allocate enough memory\n");
        MPI_Finalize();
        exit(1);
    }

    for (i = 0; i < size; i++) marked[i] = 0;

    // ==== New: root process finds all primes ≤ sqrt(n) ====
    LL sqrt_n = (LL)sqrt((double)n);
    int* small_marked = NULL;
    int* small_primes = NULL;
    int small_size = (sqrt_n - 1) / 2 + 1; // only odd numbers

    int small_count = 0;
    if (id == 0) {
        small_marked = (int*)calloc(small_size, sizeof(int));
        for (int i = 0; i < small_size; i++) small_marked[i] = 0;
        for (int i = 0; i < small_size; i++) {
            if (!small_marked[i]) {
                int prime = 2 * i + 3;
                for (LL j = (prime * prime - 3) / 2; j < small_size; j += prime) {
                    small_marked[j] = 1;
                }
            }
        }

        // count and collect the primes
        for (int i = 0; i < small_size; i++) {
            if (!small_marked[i]) small_count++;
        }

        small_primes = (int*)malloc(small_count * sizeof(int));
        int idx = 0;
        for (int i = 0; i < small_size; i++) {
            if (!small_marked[i]) {
                small_primes[idx++] = 2 * i + 3;
            }
        }
    }

    // ==== Broadcast the number of small primes, then the list ====
    MPI_Bcast(&small_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (id != 0) small_primes = (int*)malloc(small_count * sizeof(int));
    MPI_Bcast(small_primes, small_count, MPI_INT, 0, MPI_COMM_WORLD);

    // ==== Use these small primes to mark composite numbers ====
    for (int j = 0; j < small_count; j++) {
        int prime = small_primes[j];

        // Find the first multiple of prime in [low_value, high_value]
        LL first;
        if (prime * prime > low_value)
            first = prime * prime;
        else {
            LL rem = low_value % prime;
            first = (rem == 0) ? low_value : (low_value + prime - rem);
        }
        if (first % 2 == 0) first += prime;

        for (i = first; i <= high_value; i += 2 * prime) {
            marked[(i - low_value) / 2] = 1;
        }
    }

    // ==== Count the number of unmarked (prime) numbers ====
    count = 0;
    for (i = 0; i < size; i++)
        if (!marked[i]) count++;

    // ==== Reduce and output ====
    MPI_Reduce(&count, &global_count, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    elapsed_time += MPI_Wtime();

    if (!id) {
        printf("There are %lld primes less than or equal to %lld\n", global_count + 1, n);
        printf("SIEVE (%d processes) took %10.6f seconds.\n", p, elapsed_time);
    }

    MPI_Finalize();
    return 0;
}
