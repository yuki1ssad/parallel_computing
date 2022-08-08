#include <Matrix.h>
#include <vector>
#include <sys/time.h>
#include <iostream>
#include <fstream>


void Exp_One(int scale, std::ofstream &outfile) { //计算量 A*A A*A*A
    Matrix A(scale, scale);
    struct timeval start, end, para_start, para_end;
    double t, para_t;

    // A*A
    std::cout << "\n计算 A*A, scale = " << scale << std::endl;
    //串行
    gettimeofday(&start, nullptr);
    Matrix res0 = A * A;
    gettimeofday(&end, nullptr);
    t = (end.tv_sec*1000000 + end.tv_usec) - (start.tv_sec*1000000 + start.tv_usec); //um

    //并行
    for (int thrd = 2; thrd <= 32; thrd += 2) {
        gettimeofday(&para_start, nullptr);
        Matrix para_res0 = A.para_mul(A, thrd);
        gettimeofday(&para_end, nullptr);
        para_t = (para_end.tv_sec*1000000 + para_end.tv_usec) - (para_start.tv_sec*1000000 + para_start.tv_usec); //um

        std::cout << "     t = " << t/1000000 << "s\n" << "para_t = " << para_t/1000000 << "s\n" << thrd << " 线程的加速比： " << t/para_t << std::endl;
        outfile << "A*A" << "," << thrd << "," << t/para_t << std::endl;
    }

    // A*A*A
    std::cout << "\n\n计算 A*A*A, scale = " << scale << std::endl;
    //串行
    gettimeofday(&start, nullptr);
    Matrix res1 = A * A * A;
    gettimeofday(&end, nullptr);
    t = (end.tv_sec*1000000 + end.tv_usec) - (start.tv_sec*1000000 + start.tv_usec); //um

    //并行
    for (int thrd = 2; thrd <= 32; thrd += 2) {
        gettimeofday(&para_start, nullptr);
        Matrix para_res1 = A.para_mul(A, thrd).para_mul(A, thrd);
        gettimeofday(&para_end, nullptr);
        para_t = (para_end.tv_sec*1000000 + para_end.tv_usec) - (para_start.tv_sec*1000000 + para_start.tv_usec); //um

        std::cout << "     t = " << t/1000000 << "s\n" << "para_t = " << para_t/1000000 << "s\n" << thrd << " 线程的加速比： " << t/para_t << std::endl;
        outfile << "A*A*A" << "," << thrd << "," << t/para_t << std::endl;
    }

    // A*A*A*A
    std::cout << "\n\n计算 A*A*A*A, scale = " << scale << std::endl;
    //串行
    gettimeofday(&start, nullptr);
    Matrix res2 = A * A * A;
    gettimeofday(&end, nullptr);
    t = (end.tv_sec*1000000 + end.tv_usec) - (start.tv_sec*1000000 + start.tv_usec); //um

    //并行
    for (int thrd = 2; thrd <= 32; thrd += 2) {
        gettimeofday(&para_start, nullptr);
        Matrix para_res2 = A.para_mul(A, thrd).para_mul(A, thrd);
        gettimeofday(&para_end, nullptr);
        para_t = (para_end.tv_sec*1000000 + para_end.tv_usec) - (para_start.tv_sec*1000000 + para_start.tv_usec); //um

        std::cout << "     t = " << t/1000000 << "s\n" << "para_t = " << para_t/1000000 << "s\n" << thrd << " 线程的加速比： " << t/para_t << std::endl;
        outfile << "A*A*A*A" << "," << thrd << "," << t/para_t << std::endl;
    }
}

void Exp_Two(std::ofstream &outfile) { // 计算规模 scale = {128, 256, 512, 1024, 2048}
    int scale_array[] = {128, 256, 512, 1024, 2048};
    struct timeval start, end, para_start, para_end;
    double t, para_t;
    for (int &s : scale_array) {
        std::cout << "\nSCALE = " << s << std::endl;
        Matrix A(s,s);
        //串行
        gettimeofday(&start, nullptr);
        Matrix res = A*A;
        gettimeofday(&end, nullptr);
        t = (end.tv_sec*1000000 + end.tv_usec) - (start.tv_sec*1000000 + start.tv_usec); //um

        //并行
        for (int thrd = 2; thrd <= 32; thrd += 2) {
            gettimeofday(&para_start, nullptr);
            Matrix para_res = A.para_mul(A, thrd);
            gettimeofday(&para_end, nullptr);
            para_t = (para_end.tv_sec*1000000 + para_end.tv_usec) - (para_start.tv_sec*1000000 + para_start.tv_usec); //um

            std::cout << "     t = " << t/1000000 << "s\n" << "para_t = " << para_t/1000000 << "s\n" << thrd << " 线程的加速比： " << t/para_t << std::endl;
            outfile << s << "," << thrd << "," << t/para_t << std::endl;
        }
        


    }

}

int main() {
// 计算量 A*A A*A*A scale=512

    // std::ofstream outfile;
    // outfile.open("data_omp_1.csv");
    // outfile << "QUANTITY,THREAD,RATIO" << std::endl;

    // Exp_One(512, outfile);

    // outfile.close();

// 计算规模 scale = {128, 256, 512, 1024, 2048}

    std::ofstream outfile;
    outfile.open("data_omp_2.csv");
    outfile << "SCALE,THREAD,RATIO" << std::endl;

    Exp_Two(outfile);

    outfile.close();


    return 0;
}