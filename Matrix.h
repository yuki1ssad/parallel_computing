#include <vector>

class Matrix {
private:
    int _row;
    int _col;
public:
    int** _ptr;
    Matrix(int r, int c);
    Matrix(std::vector<std::vector<int>>& vec);
    ~Matrix();
    Matrix operator*(const Matrix& a);
    Matrix para_mul(const Matrix& a, int trd);
    int get_element(int r, int c) const;
    void set_element(int r, int c, int val);
    void display();
};