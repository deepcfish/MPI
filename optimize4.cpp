#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define MIN(a,b) ((a)<(b)?(a):(b))
#define LL long long

void generate_small_primes(char *small_marked, LL sqrt_n, LL *small_primes, LL *count) {
    LL i = 0;
    for (; i + 7 <= sqrt_n; i += 8) {
        small_marked[i] = 0;
        small_marked[i + 1] = 0;
        small_marked[i + 2] = 0;
        small_marked[i + 3] = 0;
        small_marked[i + 4] = 0;
        small_marked[i + 5] = 0;
        small_marked[i + 6] = 0;
        small_marked[i + 7] = 0;
    }
    for (; i <= sqrt_n; i++) {
        small_marked[i] = 0;
    }
    
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


    LL *small_primes = NULL;
    char *small_marked = NULL;
    LL small_count = 0;
    LL sqrt_n = (LL)sqrt((double)n);
    small_primes = (LL*)malloc((sqrt_n + 1) * sizeof(LL));
    if (id == 0) {
        small_marked = (char*)malloc((sqrt_n + 1) * sizeof(char));
        generate_small_primes(small_marked, sqrt_n, small_primes, &small_count);
        free(small_marked);
    }

    MPI_Bcast(&small_count, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(small_primes, small_count, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    count=size;

    if (posix_memalign((void**)&marked, 64, size) != 0) {
        printf("Cannot allocate aligned memory\n");
        MPI_Finalize();
        exit(1);
    }

    LL BLOCK_SIZE = 524288; // cache 
    for (LL block_start = 0; block_start < size; block_start += BLOCK_SIZE) {
        LL block_end = MIN(block_start + BLOCK_SIZE, size);
        
        // 初始化标记数组
        LL i = block_start;
        for (; i + 7 < block_end; i += 8) {
            marked[i] = 0;
            marked[i + 1] = 0;
            marked[i + 2] = 0;
            marked[i + 3] = 0;
            marked[i + 4] = 0;
            marked[i + 5] = 0;
            marked[i + 6] = 0;
            marked[i + 7] = 0;
        }
        for (; i < block_end; i++) {
            marked[i] = 0;
        }
        
        // 执行素数筛选过程
        for (int j = 1; j < small_count; j++) {
            const LL prime = small_primes[j];
            LL start_index;
            LL value = low_value + block_start * 2;
            
            if (prime * prime > value)
                start_index = (prime * prime - low_value) / 2;
            else {
                LL rem = (value) % prime;
                LL first = (rem == 0) ? value : (value + prime - rem);
                if (first % 2 == 0) first += prime;
                start_index = (first - low_value) / 2;
            }
            LL step = 2 * prime;
            LL k = start_index
            for (; k < block_end; k += step) {
                if (!marked[k]) {
                    marked[k] = 1;
                    count-=2;
                }
            }
            if(k-block_end>prime) count++;
        }
    }
    
    

    //for (i = 0; i < size; i++)
        //if (!marked[i]) count++;

    MPI_Reduce(&count, &global_count, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    elapsed_time += MPI_Wtime();

    if (!id) {
        //printf("%lld",small_count);
        printf("There are %lld primes less than or equal to %lld\n", global_count+1, n);
        printf("SIEVE (%d processes) took %10.6f seconds.\n", p, elapsed_time);
    }

    MPI_Finalize();
    return 0;
}