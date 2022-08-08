#include <iostream> 
#include <stdlib.h> 
#include <time.h>  
#include <vector>
#include <string.h>
#include <omp.h>
#include <Matrix.h>

Matrix::Matrix(int r, int c) {
    _row = r;
    _col = c;
    _ptr = (int**)malloc(_row * sizeof(int*));
    for (int i = 0; i < _row; i++) {
        _ptr[i] = (int*)malloc(_col * sizeof(int));
    }

    srand (time(NULL));
    for (int i = 0; i < _row; i++) {
        for (int j = 0; j < _col; j++) {
            _ptr[i][j] = rand() % 10;
        }
    }
}

Matrix::Matrix(std::vector<std::vector<int>>& vec) {
    this->_row = vec.size();
    this->_col = vec[0].size();
    _ptr = (int**)malloc(_row * sizeof(int*));
    for (int i = 0; i < _row; i++) {
        _ptr[i] = (int*)malloc(_col * sizeof(int));
    }

    for (int i = 0; i < _row; i++) {
        memcpy(_ptr[i], vec[i].data(),  sizeof(int) * _col);
    }
}

Matrix::~Matrix() {
    for (int i = 0; i < _row; i++) {
        free(_ptr[i]);
    }
    free(_ptr);
}

int Matrix::get_element(int r, int c) const {
    return _ptr[r][c];
}

void Matrix::set_element(int r, int c, int val) {
    _ptr[r][c] = val;
}

void Matrix::display() {
    std::cout << std::endl;
    for (int i = 0; i < _row; i++) {
        for (int j = 0; j < _col; j++) {
            std::cout << _ptr[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

Matrix Matrix::operator*(const Matrix& a) {

    Matrix m(_row, a._col);
    for (int i = 0; i < this->_row; i++) {
        for (int j = 0; j < a._col; j++) {
            int sum = 0;
            for (int k = 0; k < a._row; k++) {
                sum += _ptr[i][k] * a.get_element(k, j);
            }
            m.set_element(i, j, sum);
        }
    }
    return m;
}

Matrix Matrix::para_mul(const Matrix& a, int trd) {
    Matrix m(_row, a._col);

    omp_set_num_threads(trd);

    std::cout << "trd = " << trd << std::endl;
    #pragma omp parallel for

    for (int i = 0; i < this->_row; i++) {
        for (int j = 0; j < a._col; j++) {
            int sum = 0;
            for (int k = 0; k < a._row; k++) {
                sum += _ptr[i][k] * a.get_element(k, j);
            }
            m.set_element(i, j, sum);
        }
    }
    return m;
}
