#ifndef HELMHOLTZ_H
#define HELMHOLTZ_H

#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <chrono>
 

using namespace std;
class helmholtz{
    public:
        helmholtz(double L,int N);
        ~helmholtz();
        void SetGrid();
        void SetExactSolution();
        void SetRightHandSide();
        void SetEigenValues();
        void SetPandPinv();
        void ComputeFtilde();
        void ComputeUtilde();
        void ComputeU();
        void L2Error();

        void Solve();
        

    private:
        double L;
        int N;
        double h;
        double  *x,*y;

        double  **u,**uexact,**P,**Pinverse;                
        double  **b,**btilde,**temp;
        double  *eigenvalue;
        double  l2error, time;

        void CreateMatrix(double**& A);
        void DeleteMatrix(double**& A);
        void WriteDataTecplot2D(double**& A,string folderpath,string filename);
        void WriteDataTecplot1D(double*& b,string folderpath,string filename);

};

#endif

