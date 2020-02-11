#ifndef FIELD_H
#define FIELD_H

#include <iostream>
#include <fstream>
#include <cassert>
#include <string>
#include <cmath>

using namespace std;

class field{
    public:
        field(double lx,double ly,int nx,int ny);
        field(const field& v);
        ~field();
        void set_value(const double value);
        void set_function(double (*func)(const double&,const double&));
        void write_matplolib(string folderpath, string filename);
        double hx()const;
        double hy()const;
        int    nx()const;
        int    ny()const;
        
        //operators:
        double& operator()(const int j,const int i);
        field&  operator=(const field&);

        
        friend double L2(field& u);
        //friend double Linfy(const field& u);

    private:
        double** u;
        double Lx,Ly,dx,dy;
        int    Nx,Ny;


};

#endif
