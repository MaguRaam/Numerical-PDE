#include "../include/Matrix.h"

//Constructor:
Matrix::Matrix(int nrows, int ncols):nrows(nrows),ncols(ncols){
    assert(nrows > 0); assert(ncols > 0);
    A = new double*[nrows];
    for (int i = 0;i<nrows;i++){
        A[i] = new double[ncols];
        for (int j = 0;j<ncols;j++)
            A[i][j] = 0.0;
    }

}

//Copy Constructor:
Matrix::Matrix(const Matrix& OtherMatrix):nrows(OtherMatrix.nrows),ncols(OtherMatrix.ncols){
    assert(nrows > 0); assert(ncols > 0);
    A = new double*[nrows];
    for (int i = 0;i<nrows;i++){
        A[i] = new double[ncols];
        for (int j = 0;j<ncols;j++)
            A[i][j] = OtherMatrix.A[i][j];
    }
}


//Destructor:
Matrix::~Matrix(){
    for (int i = 0;i<nrows;i++){
        delete [] A[i];
    }
    delete [] A;
}

//size:
int Matrix::Rows() const{
    return nrows;
}


int Matrix::Cols() const{
    return ncols;
}

//access elements ():
double& Matrix::operator()(int i,int j){
    assert(i > -1); assert(i < nrows);
    assert(j > -1); assert(j < ncols);
    return A[i][j];
}

//Matrix-Matrix multiplication:
Matrix Matrix::operator*(const Matrix& OtherMatrix) const{
    assert(ncols == OtherMatrix.nrows);
    Matrix mat(nrows,OtherMatrix.ncols);
    for (int i = 0;i<nrows;i++){
        for (int j = 0;j<OtherMatrix.ncols;j++){
            for (int k=0;k<ncols;k++)
                mat(i,j) += A[i][k]*OtherMatrix.A[k][j];
        }
    }

    return mat;
}

 