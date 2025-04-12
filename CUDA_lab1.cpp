#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <cuda.h>
#include <cuda_runtime.h>

#define N 10000000
#define MAX_ERR 1e-6

// CUDA kernel for element-wise vector addition
__global__ void vector_add(float *out, float *a, float *b, int n) {
    //int index = threadIdx.x;
    //int stride = blockDim.x;

    // Each thread processes multiple elements with stride access pattern
    //for (int i = index; i < n; i += stride) {
        //out[i] = a[i] + b[i];
    //}
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if(tid<n){
        out[tid]=a[tid]+b[tid];
    }
}

int main() {
    float *a, *b, *out;
    float *d_a, *d_b, *d_out;

    // Allocate memory on the host (CPU)
    a   = (float *)malloc(sizeof(float) * N);
    b   = (float *)malloc(sizeof(float) * N);
    out = (float *)malloc(sizeof(float) * N);

    // Initialize host arrays
    for (int i = 0; i < N; i++) {
        a[i] = 1.0f;
        b[i] = 2.0f;
    }

    // Allocate memory on the device (GPU)
    cudaMalloc((void **)&d_a, sizeof(float) * N);
    cudaMalloc((void **)&d_b, sizeof(float) * N);
    cudaMalloc((void **)&d_out, sizeof(float) * N);

    // Copy input data from host to device
    cudaMemcpy(d_a, a, sizeof(float) * N, cudaMemcpyHostToDevice);
    cudaMemcpy(d_b, b, sizeof(float) * N, cudaMemcpyHostToDevice);

    int block_size=256;
    int grid_size=(N+block_size)/block_size;
    // Launch the kernel with 1 block and 256 threads
    vector_add<<<grid_size, block_size>>>(d_out, d_a, d_b, N);

    // Copy result back from device to host
    cudaMemcpy(out, d_out, sizeof(float) * N, cudaMemcpyDeviceToHost);

    // Verify results
    for (int i = 0; i < N; i++) {
        assert(fabs(out[i] - a[i] - b[i]) < MAX_ERR);
    }

    printf("PASSED\n");

    // Free device memory
    cudaFree(d_a);
    cudaFree(d_b);
    cudaFree(d_out);

    // Free host memory
    free(a);
    free(b);
    free(out);
}
