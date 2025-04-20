#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define LL long long

void simple_sieve(LL limit, LL **primes, LL *count) {
    char *is_prime = malloc((limit + 1) * sizeof(char));
    for (LL i = 0; i <= limit; i++) is_prime[i] = 1;
    is_prime[0] = is_prime[1] = 0;

    for (LL i = 2; i * i <= limit; i++) {
        if (is_prime[i]) {
            for (LL j = i * i; j <= limit; j += i)
                is_prime[j] = 0;
        }
    }

    *count = 0;
    for (LL i = 2; i <= limit; i++)
        if (is_prime[i]) (*count)++;

    *primes = malloc((*count) * sizeof(LL));
    LL idx = 0;
    for (LL i = 2; i <= limit; i++)
        if (is_prime[i]) (*primes)[idx++] = i;

    free(is_prime);
}

int main(int argc, char* argv[]) {
    LL count = 0, global_count = 0;
    double elapsed_time;
    LL low_index, high_index, size;
    LL i, prime, first;
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

    // Only odd numbers > 2 are considered
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

    // Only process 0 computes the small prime list
    LL *small_primes = NULL;
    LL prime_count = 0;
    if (id == 0) {
        simple_sieve((LL)sqrt(n), &small_primes, &prime_count);
    }

    // Broadcast prime count and primes to all processes
    MPI_Bcast(&prime_count, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    if (id != 0) small_primes = malloc(prime_count * sizeof(LL));
    MPI_Bcast(small_primes, prime_count, MPI_LONG_LONG, 0, MPI_COMM_WORLD);

    for (LL k = 0; k < prime_count; k++) {
        prime = small_primes[k];

        // Mapping from index to value: value = 2 * idx + 3
        LL low_value = 2 * low_index + 3;

        if (prime * prime > low_value) {
            first = (prime * prime - 3) / 2 - low_index;
        } else {
            LL rem = low_value % prime;
            LL delta = (rem == 0) ? 0 : prime - rem;
            if ((low_value + delta) % 2 == 0) delta += prime;
            first = (delta + 1) / 2;
        }

        for (i = first; i < size; i += prime)
            marked[i] = 1;
    }

    count = 0;
    for (i = 0; i < size; i++)
        if (!marked[i]) count++;

    MPI_Reduce(&count, &global_count, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    elapsed_time += MPI_Wtime();

    if (!id) {
        global_count += 1; // Include prime number 2
        printf("There are %lld primes less than or equal to %lld\n", global_count, n);
        printf("Optimized SIEVE (%d processes): %10.6f seconds\n", p, elapsed_time);
    }

    free(marked);
    free(small_primes);
    MPI_Finalize();
    return 0;
}
