#include <Matrix.h>
#include <mpi.h>
#include <stdlib.h> 
#include <iostream>
#include <time.h>
#include <omp.h>
#include <fstream>

void print_mat(int* M, int row, int col) {
    std::cout << std::endl;
    for (int i = 0; i <row; i++) {
        for (int j = 0; j < col; j++) {
            std::cout << M[i * col + j] <<" ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void Exp(Matrix &A, Matrix &B, int thrd, std::fstream &outfile) {

    int r_1 = A.get_row();
    int c_1 = A.get_col();
    int r_2 = B.get_row();
    int c_2 = B.get_col();

    if (r_2 != c_1) {
        std::cout << "error input, r_2 != c_1" <<std::endl;
        exit(EXIT_FAILURE);
    }

    // A.display();
    // B.display();

    int my_rank;
    int numprocs;
    int line;
    int *mat_a, *mat_b, *buffer_a, *ans, *res;
    double start, end, proc_t; // proc_t 表示多进程耗时

    // MPI环境初始化
    MPI_Init(nullptr,nullptr);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    // std::cout << "当前进程号：" << my_rank << std::endl;

    int r_1_s = r_1;
    if (r_1 % numprocs != 0) {
        r_1 -= r_1 % numprocs;
        r_1 += numprocs;
    }

    line = r_1 / numprocs; //每个进程处理的行数
    // std::cout << "line = " << line << std::endl;

    mat_a = (int *)malloc(r_1 * c_1 * sizeof(int));
    mat_b = (int *)malloc(r_2 * c_2 * sizeof(int));
    buffer_a = (int *)malloc(line * c_1 * sizeof(int));
    ans = (int *)malloc(line * c_2 * sizeof(int));
    res = (int *)malloc(r_1 * c_2 * sizeof(int));

    if (mat_a == nullptr || buffer_a == nullptr || ans == nullptr || res == nullptr) {
        std::cout << "no enough memory" <<std::endl;
        exit(EXIT_FAILURE);
    }

    if (my_rank == 0) { // 0号进程负责初始化 mat_a, mat_b; 并开始计时
        //将 A 矩阵拉平为 mat_a
        for (int i = 0; i < r_1; i++) {
            for (int j = 0; j < c_1; j++) {
                if (i < r_1_s){
                    mat_a[i * c_1 + j] = A.get_element(i, j);
                } else {
                    mat_a[i * c_1 + j] = 0;
                }  
            }
        }
        //将 B 矩阵拉平为 mat_b
        for (int i = 0; i < r_2; i++) {
            for (int j = 0; j < c_2; j++) {
                mat_b[i * c_2 + j] = B.get_element(i, j);
            }
        }

        // print_mat(mat_a, r_1, c_1);
        // B.display();

        //开始计时
        start = MPI_Wtime();
    }

    MPI_Scatter(mat_a, line * c_1, MPI_INT, buffer_a, line * c_1, MPI_INT, 0, MPI_COMM_WORLD); // 数据分发(mat_a)

    // std::cout << "buffer_a:\n";
    // print_mat(buffer_a, line, c_1);

    MPI_Bcast(mat_b, r_2 * c_2, MPI_INT, 0, MPI_COMM_WORLD); //数据广播(mat_b)

    int tmp = 0;
    #pragma omp parallel for num_threads(thrd) // 对最外层for并行
    for (int i = 0; i < line; i++) {
        for (int j = 0; j < c_2; j++) {
            tmp = 0;
            for (int k = 0; k < c_1; k++) {
                tmp += buffer_a[i * c_1 + k] * mat_b[k * c_2 + j];
            }
            ans[i * c_2 + j] = tmp;

            // print ans
            // print_mat(ans, line, c_2); 
            // std::cout << "i = " << i <<" j = " << j << "的ans" << std::endl;
        }
    }

    MPI_Gather(ans, line * c_2, MPI_INT, res, line * c_2, MPI_INT, 0, MPI_COMM_WORLD); // 数据收集到0号进程

    if (my_rank == 0) { // 0号进程结束计时，并输出结果
        // print_mat(res, r_1, c_2); // 打印计算结果
        end = MPI_Wtime();
        proc_t = end - start; 
        std::cout << "proc_t: " << proc_t  << " s\n" <<std::endl;  
        outfile << std::endl << numprocs << "," << thrd << "," << proc_t;
    }
    
    // 释放申请的内存
    free(mat_a);
    free(mat_b);
    free(buffer_a);
    free(ans);
    free(res);

    //结束MPI环境
    MPI_Finalize();
}

int main(int argc, char* argv[]) {

    if (argc != 2) { //判断参数输入是否正确
        std::cout << "Two arguments are needed, please enter the thread number.\n";
        exit(EXIT_FAILURE);
    }

    int thrd = atoi(argv[1]); // thrd 接收传入的线程数
    int scale = 2048;
    Matrix A(scale, scale);
    Matrix B(scale, scale);

    std::fstream outfile;
    outfile.open("data_mpi.csv",std::ios::out | std::ios::app);

    Exp(A, B, thrd, outfile);

    outfile.close();
    
    return 0;
}