#ifndef MATRIX_H
#define MATRIX_H

#include <cassert>
#include "Vector.h"

class Matrix
{
public:
    Matrix(int nrows, int ncols);
    Matrix(const Matrix& OtherMatrix);
    ~Matrix();
    int Rows() const;
    int Cols() const;
    double& operator()(int i,int j);
    Matrix operator*(const Matrix& OtherMatrix) const;
     
    

protected:
    int nrows, ncols;
    double **A;
};

#endif