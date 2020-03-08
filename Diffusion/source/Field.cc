#include "../include/Field.h"

Field::Field(int Nx,int Ny,double Lx,double Ly) :
                       Matrix(Nx+1,Ny+1),
                       Grid(Nx,Ny,Lx,Ly)
{
}

Field::Field(const Field& OtherField):Matrix(OtherField),Grid(OtherField)
{}


void Field::WriteDataTecplot(std::string folderpath,std::string filename){
    std::ofstream File(folderpath + filename, std::ios::out);
    File.flags( std::ios::dec | std::ios::scientific );
    File.precision(16) ;
    if(!File)
        std::cerr<<"Error opening file";
    File<<"TITLE = "+ filename<<"\n"<<"VARIABLES = X, Y," + filename<<"\n";
    File<<"Zone T = psi I = "<<Nx+1<<"J = "<<Ny+1<<"\n";
    for (int i = 0;i<=Nx;i++){
        for (int j = 0;j<=Ny;j++)
            File<<(*x)(i)<<"\t"<<(*y)(j)<<"\t"<<A[i][j]<<"\t"<<"\n";
    }
    File.close();
}

void Field::SetFunction(double (*func) (double x,double y,double t) , double t ){
    for (int i = 0;i<=Nx;i++){
        for (int j = 0;j<=Ny;j++)
            A[i][j] = func( (*x)(i) , (*y)(j),t);
    }
}


//TODO:Check Matrix Size:

//Matrix-Matrix multiplication:
Field Field::operator*(const Field& OtherField) const{
    assert(ncols == OtherField.nrows);
    Field mat(nrows-1,OtherField.ncols-1,Lx,Ly);
    for (int i = 0;i<nrows;i++){
        for (int j = 0;j<OtherField.ncols;j++){
            for (int k=0;k<ncols;k++)
                mat(i,j) += A[i][k]*OtherField.A[k][j];
        }
    }

    return mat;
} 