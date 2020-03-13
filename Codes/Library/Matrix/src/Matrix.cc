#include "../include/Matrix.h"

#ifndef MATRIX_CC
#define MATRIX_CC

//Constructor
template <typename T>
Matrix<T>::Matrix(unsigned int rows, unsigned int cols, const T &value) : rows(rows),
                                                                          cols(cols)
{
    mat.resize(rows);
    for (unsigned int i = 0; i < mat.size(); i++)
        mat[i].resize(cols, value);
}

//Copy Constructor
template <typename T>
Matrix<T>::Matrix(const Matrix<T> &rhs) : rows(rhs.rows),
                                          cols(rhs.cols),
                                          mat(rhs.mat)
{
}

//Destructor
template <typename T>
Matrix<T>::~Matrix() {}

//Operator Overloading

//Assignment opearator = TODO Optimize
template <typename T>
Matrix<T> &Matrix<T>::operator=(const Matrix<T> &rhs)
{
    if (&rhs == this)
        return *this;

    assert(this->rows == rhs.rows);
    assert(this->cols == rhs.cols);
    for (unsigned i = 0; i < this->rows; i++)
    {
        for (unsigned j = 0; j < this->cols; j++)
        {
            mat[i][j] = rhs.mat[i][j];
        }
    }

    return *this;
}

//Matrix Addition operator +
template <typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T> &rhs)
{
    assert(this->rows == rhs.rows);
    assert(this->cols == rhs.cols);

    Matrix<T> result(rows, cols, 0.0);
    for (unsigned i = 0; i < rows; i++)
    {
        for (unsigned j = 0; j < cols; j++)
        {
            result(i, j) = this->mat[i][j] + rhs(i, j);
        }
    }

    return result;
}

//Matrix Subtraction operator -
template <typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T> &rhs)
{
    assert(this->rows == rhs.rows);
    assert(this->cols == rhs.cols);

    Matrix<T> result(rows, cols, 0.0);
    for (unsigned i = 0; i < rows; i++)
    {
        for (unsigned j = 0; j < cols; j++)
        {
            result(i, j) = this->mat[i][j] - rhs(i, j);
        }
    }

    return result;
}

//Matrix Multiplication opearator *
template <typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T> &rhs) const
{
    assert(this->cols == rhs.rows);

    Matrix result(this->rows, rhs.get_cols(), 0.0);
    for (unsigned i = 0; i < this->rows; i++)
    {
        for (unsigned j = 0; j < rhs.get_cols(); j++)
        {
            for (unsigned k = 0; k < this->cols; k++)
            {
                result(i, j) += this->mat[i][k] * rhs(k, j);
            }
        }
    }

    return result;
}

// Calculate a transpose of this matrix
template <typename T>
Matrix<T> Matrix<T>::transpose() const
{
    Matrix result(rows, cols, 0.0);

    for (unsigned i = 0; i < rows; i++)
    {
        for (unsigned j = 0; j < cols; j++)
        {
            result(i, j) = this->mat[j][i];
        }
    }

    return result;
}

//Operator Overloading Matrix-Scalar operation:
// Matrix/scalar addition
template <typename T>
Matrix<T> Matrix<T>::operator+(const T &rhs)
{
    Matrix result(rows, cols, 0.0);

    for (unsigned i = 0; i < rows; i++)
    {
        for (unsigned j = 0; j < cols; j++)
        {
            result(i, j) = this->mat[i][j] + rhs;
        }
    }

    return result;
}

// Matrix/scalar subtraction
template <typename T>
Matrix<T> Matrix<T>::operator-(const T &rhs)
{
    Matrix result(rows, cols, 0.0);

    for (unsigned i = 0; i < rows; i++)
    {
        for (unsigned j = 0; j < cols; j++)
        {
            result(i, j) = this->mat[i][j] - rhs;
        }
    }

    return result;
}

// Matrix/scalar multiplication
template <typename T>
Matrix<T> Matrix<T>::operator*(const T &rhs)
{
    Matrix result(rows, cols, 0.0);

    for (unsigned i = 0; i < rows; i++)
    {
        for (unsigned j = 0; j < cols; j++)
        {
            result(i, j) = this->mat[i][j] * rhs;
        }
    }

    return result;
}

// Matrix/scalar division
template <typename T>
Matrix<T> Matrix<T>::operator/(const T &rhs)
{
    Matrix result(rows, cols, 0.0);

    for (unsigned i = 0; i < rows; i++)
    {
        for (unsigned j = 0; j < cols; j++)
        {
            result(i, j) = this->mat[i][j] / rhs;
        }
    }

    return result;
}

//Access the individual element
template <typename T>
T &Matrix<T>::operator()(const unsigned int &i, const unsigned int &j)
{
    return this->mat[i][j];
}

template <typename T>
const T &Matrix<T>::operator()(const unsigned int &i, const unsigned int &j) const
{
    return this->mat[i][j];
}

//Access row and col size

template <typename T>
unsigned int Matrix<T>::get_rows() const
{
    return this->rows;
}

template <typename T>
unsigned int Matrix<T>::get_cols() const
{
    return this->cols;
}

#endif
