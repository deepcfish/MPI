#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define MIN(a,b)  ((a)<(b)?(a):(b))
#define LL long long

// Generate small odd primes up to sqrt(n)
void generate_small_primes(char *small_marked, LL sqrt_n, LL *small_primes, LL *count) {
    for (LL i = 0; i <= sqrt_n; ++i) small_marked[i] = 0;
    for (LL i = 3; i * i <= sqrt_n; i += 2) {
        if (!small_marked[i]) {
            for (LL j = i * i; j <= sqrt_n; j += 2 * i)
                small_marked[j] = 1;
        }
    }
    *count = 0;
    small_primes[(*count)++] = 2;
    for (LL i = 3; i <= sqrt_n; i += 2) {
        if (!small_marked[i]) {
            small_primes[(*count)++] = i;
        }
    }
}

int main(int argc, char* argv[]) {
    LL count;
    double elapsed_time;
    LL global_count = 0;
    LL high_index, low_index;
    LL i, j;
    int id;
    char *marked;
    LL n;
    int p;
    LL size;
    const LL BLOCK_SIZE = 32768;

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

    // Only store odd numbers > 2: index i represents number 2*i + 3
    LL total_odd_count = (n - 1) / 2;
    low_index = id * total_odd_count / p;
    high_index = ((id + 1) * total_odd_count / p) - 1;
    size = high_index - low_index + 1;

    marked = (char*)malloc(size);
    for (i = 0; i < size; i++) marked[i] = 0;

    char *small_marked = (char*)malloc((sqrt_n + 1));
    LL *small_primes = (LL*)malloc((sqrt_n + 1) * sizeof(LL));
    LL small_count = 0;
    generate_small_primes(small_marked, sqrt_n, small_primes, &small_count);
    free(small_marked);

    // Perform blocked sieve
    for (LL block_low = 0; block_low < size; block_low += BLOCK_SIZE) {
        LL block_high = MIN(block_low + BLOCK_SIZE, size);
        for (j = 0; j < small_count; ++j) {
            LL prime = small_primes[j];
            LL prime2 = 2 * prime;
            LL first_val = 2 * low_index + 3 + 2 * block_low;
            LL first;

            if (prime * prime > first_val)
                first = (prime * prime - 3) / 2 - low_index;
            else {
                LL r = first_val % prime;
                if (r == 0) first = block_low;
                else {
                    LL delta = prime - r;
                    if ((first_val + delta) % 2 == 0) delta += prime;
                    first = (first_val + delta - 3) / 2 - low_index;
                }
            }
            for (i = first; i < block_high; i += prime2) marked[i] = 1;
        }
    }

    free(small_primes);
    count = 0;
    for (i = 0; i < size; i++)
        if (!marked[i]) count++;
    free(marked);

    if (2 <= n) count++; // account for prime number 2

    MPI_Reduce(&count, &global_count, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    elapsed_time += MPI_Wtime();

    if (!id) {
        printf("There are %lld primes <= %lld\n", global_count, n);
        printf("SIEVE (no evens, no Bcast, %d processes) %10.6f seconds\n", p, elapsed_time);
        printf("the process id",id);
    }
    fflush(stdout);
    MPI_Finalize();
    return 0;
}
