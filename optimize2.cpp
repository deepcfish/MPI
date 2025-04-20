#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define MIN(a,b)  ((a)<(b)?(a):(b))
#define LL long long

// Function to generate small primes up to sqrt(n) using serial sieve
void generate_small_primes(char *small_marked, LL sqrt_n, LL *small_primes, LL *count) {
    for (LL i = 0; i <= sqrt_n; ++i) small_marked[i] = 0;
    for (LL i = 2; i * i <= sqrt_n; ++i) {
        if (!small_marked[i]) {
            for (LL j = i * i; j <= sqrt_n; j += i)
                small_marked[j] = 1;
        }
    }
    *count = 0;
    for (LL i = 2; i <= sqrt_n; ++i) {
        if (!small_marked[i]) {
            small_primes[(*count)++] = i;
        }
    }
}

int main(int argc, char* argv[]) {
    LL    count;
    double elapsed_time;
    LL    global_count = 0;
    LL    high_value;
    LL    i;
    int   id;
    LL    low_value;
    char* marked;
    LL    n;
    int   p;
    LL    size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime();

    if (argc != 2) {
        if (!id) printf("Usage: %s <n>\n", argv[0]);
        MPI_Finalize();
        exit(1);
    }

    n = atoll(argv[1]);
    LL sqrt_n = (LL)sqrt((double)n);

    low_value = 2 + id * (n - 1) / p;
    high_value = 1 + (id + 1) * (n - 1) / p;
    size = high_value - low_value + 1;

    marked = (char*)malloc(size);
    for (i = 0; i < size; i++) marked[i] = 0;

    // Step 1: All processes generate small primes
    char *small_marked = (char*)malloc((sqrt_n + 1));
    LL *small_primes = (LL*)malloc((sqrt_n + 1) * sizeof(LL));
    LL small_count = 0;

    generate_small_primes(small_marked, sqrt_n, small_primes, &small_count);
    free(small_marked);

    // Step 2: Each process uses small primes to mark its own range
    for (LL j = 0; j < small_count; ++j) {
        LL prime = small_primes[j];
        LL first;
        if (prime * prime > low_value)
            first = prime * prime - low_value;
        else {
            if (!(low_value % prime)) first = 0;
            else first = prime - (low_value % prime);
        }
        for (i = first; i < size; i += prime) marked[i] = 1;
    }
    free(small_primes);

    count = 0;
    for (i = 0; i < size; i++)
        if (!marked[i]) count++;
    free(marked);

    MPI_Reduce(&count, &global_count, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    elapsed_time += MPI_Wtime();

    if (!id) {
        printf("There are %lld primes <= %lld\n", global_count, n);
        printf("SIEVE (no Bcast, %d processes) %10.6f seconds\n", p, elapsed_time);
    }

    MPI_Finalize();
    return 0;
}
