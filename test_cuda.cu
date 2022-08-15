#include <random>
#include <cuda.h>
#include <sys/time.h>
#include <iostream>

#define M 4096
#define N 4096
#define K 4096
#define THREAD_NUM 1024

void mat_init (int *mat_a, int *mat_b, int *mat_res, int m, int n, int k) {
    srand(0);
    // mat_a
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            mat_a[i * m + j] = random() % 10;
        }
    }
    // mat_b
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < k; j++) {
            mat_b[i * n + j] = random() % 10;
        }
    }
    // mat_res
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < k; j++) {
            mat_res[i * m + j] = 0;
        }
    }
}

void cpu_mat_mul (int *mat_a, int *mat_b, int *mat_res, int m, int n, int k) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < k; j++) {
            int tmp = 0;
            for (int q = 0; q < n; q++) {
                tmp += mat_a[i * n + q] * mat_b[q * k + j];
            }
            mat_res[i * m + j] = tmp;
        }
    }
}

__global__ void gpu_mat_mul (int *mat_a, int *mat_b, int *mat_res, int m, int n, int k) {
    int tid = threadIdx.x;
    int bid = blockIdx.x;
    int idx = bid * THREAD_NUM + tid;

    int row = idx / m;
    int col = idx % m;

    if (row < m && col < n) {
        int tmp = 0;
        for (int q = 0; q < n; q++) {
            tmp += mat_a[row * n + q] * mat_b[q * k + col];
        }
        mat_res[row * m + col] = tmp;
    }
}

void print_mat(int *mat, int row, int col) {
    std::cout << std::endl;
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            std::cout << mat[i * col + j] <<" ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

int main() {
    int *a, *b, *res;
    int m = M, n = N, k = K;

    // 分配CPU内存
    a = (int *)malloc(m * n * sizeof(int));
    b = (int *)malloc(n * k * sizeof(int));
    res = (int *)malloc(m * k * sizeof(int));
    // 矩阵初始化
    mat_init(a, b, res, m, n, k);
    // print_mat(a, m, n);
    // print_mat(b, n, k);

    struct timeval cpu_start, cpu_end, gpu__start, gpu_end;
    double cpu_time, gpu_time;

    // cpu_mat_mul
    gettimeofday(&cpu_start, nullptr);
    cpu_mat_mul(a, b, res, m, n, k);
    gettimeofday(&cpu_end, nullptr);
    cpu_time = (cpu_end.tv_sec*1000000 + cpu_end.tv_usec) - (cpu_start.tv_sec*1000000 + cpu_start.tv_usec); //um

    // 打印结果
    // print_mat(res, m, k);

    // gpu_mat_mul
    //从主机拷贝数据到设备
    int *dev_a, *dev_b, *dev_res;
    cudaMalloc((void **)&dev_a, m * n * sizeof(int));
    cudaMalloc((void **)&dev_b, n * k * sizeof(int));
    cudaMalloc((void **)&dev_res, m * k * sizeof(int));
    cudaMemcpy(dev_a, a, m * n * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_b, b, n * k * sizeof(int), cudaMemcpyHostToDevice);

    int block_num = (m * n - 1) / THREAD_NUM + 1;

    gettimeofday(&gpu__start, nullptr);
    gpu_mat_mul<<<block_num, THREAD_NUM>>>(dev_a, dev_b,dev_res, m, n, k);
    cudaMemcpy(res, dev_res, m * k * sizeof(int), cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    gettimeofday(&gpu_end, nullptr);
    gpu_time = (gpu_end.tv_sec*1000000 + gpu_end.tv_usec) - (gpu__start.tv_sec*1000000 + gpu__start.tv_usec); //um

    // 打印结果
    // print_mat(res, m, k);

    // 输出结果
    std::cout << "矩阵规模：[" << m  <<"][" << n << "] * [" << n << "][" << k << "]\n";
    std::cout << "cpu_time = " << cpu_time/1000000 << " s\n" << "gpu_time = " << gpu_time/1000000 << " s\n" << "加速比 = " << cpu_time/gpu_time << std::endl;

    free(a);
    free(b);
    free(res);
    cudaFree(dev_a);
    cudaFree(dev_b);
    cudaFree(dev_res);

    return 0;
}