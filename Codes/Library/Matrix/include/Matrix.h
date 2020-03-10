#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <cassert>

template <typename T>
class Matrix
{
private:
    unsigned int rows;
    unsigned int cols;
    std::vector<std::vector<T>> mat;

public:
    //Constructor and destructor
    Matrix(unsigned int rows, unsigned int cols, const T &value);
    Matrix(const Matrix<T> &rhs);
    ~Matrix();

    //Operator Overloading Matrix-Matrix Operation
    Matrix<T>& operator=(const Matrix<T> &rhs);
    Matrix<T> operator+(const Matrix<T> &rhs);
    Matrix<T> operator-(const Matrix<T> &rhs);
    Matrix<T> operator*(const Matrix<T> &rhs);
    Matrix<T> transpose();

    //Operator Overloading Matrix-Scalar operation: //TODO Optimise
    Matrix<T> operator+(const T& rhs);
    Matrix<T> operator-(const T& rhs);
    Matrix<T> operator*(const T& rhs);
    Matrix<T> operator/(const T& rhs);




    //Access the individual element
    T& operator()(const unsigned int &i, const unsigned int &j);
    const T& operator()(const unsigned int &i, const unsigned int &j) const;

    //Access row and col size
    unsigned int get_rows() const;
    unsigned int get_cols() const;
};

#include "../src/Matrix.cc"

#endif
