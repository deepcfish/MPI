#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define MIN(a,b)  ((a)<(b)?(a):(b))
#define LL long long

int main(int argc, char* argv[]) {
    LL count = 0, global_count = 0;
    double elapsed_time;
    LL low_index, high_index, size;
    LL i, prime, index, first;
    int id, p;
    char *marked;
    LL n;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime();

    if (argc != 2) {
        if (!id) printf("Command line: %s <n>\n", argv[0]);
        MPI_Finalize();
        exit(1);
    }

    n = atoll(argv[1]);
    if (n < 2) {
        if (!id) printf("There are no primes less than or equal to %lld\n", n);
        MPI_Finalize();
        return 0;
    }

    // Only consider odd numbers: 3, 5, 7, ..., n
    LL odd_count = (n - 1) / 2; // exclude even numbers, index 0 maps to 3
    low_index = id * odd_count / p;
    high_index = (id + 1) * odd_count / p - 1;
    size = high_index - low_index + 1;

    marked = (char*) malloc(size);
    if (marked == NULL) {
        printf("Process %d: Cannot allocate memory\n", id);
        MPI_Finalize();
        exit(1);
    }
    for (i = 0; i < size; i++) marked[i] = 0;

    if (!id) index = 0;
    prime = 3;

    do {
        // prime index in odd space: (prime - 3)/2
        LL prime_index = (prime - 3) / 2;

        if (prime * prime > (2 * low_index + 3))
            first = (prime * prime - 3) / 2 - low_index;
        else {
            LL low_value = 2 * low_index + 3;
            if (low_value % prime == 0)
                first = 0;
            else
                first = prime - (low_value % prime);
            if (((2 * low_index + 3 + first) % 2) == 0) first++;
            first = (first + 1) / 2;
        }

        for (i = first; i < size; i += prime)
            marked[i] = 1;

        if (!id) {
            while (marked[++index]);
            prime = 2 * index + 3;
        }
        if (p > 1) MPI_Bcast(&prime, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);

    } while (prime * prime <= n);

    count = 0;
    for (i = 0; i < size; i++)
        if (!marked[i]) count++;

    MPI_Reduce(&count, &global_count, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    elapsed_time += MPI_Wtime();

    if (!id) {
        // Add 1 for the prime number 2
        global_count += 1;
        printf("There are %lld primes less than or equal to %lld\n", global_count, n);
        printf("Optimized SIEVE (%d processes): %10.6f seconds\n", p, elapsed_time);
    }

    free(marked);
    MPI_Finalize();
    return 0;
}