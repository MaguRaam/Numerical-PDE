#ifndef FIELDHEADERDEF
#define FIELDHEADERDEF

#include <iostream>
#include <cassert>
#include <cmath>
#include <fstream>

class field
{
public:
    field(double Lengthx, double Lengthy, int ngridx, int ngridy);
    field(const field &other_field);
    void set_value(double value);
    void set_function(double (*func)(double, double, double), double t);
    ~field();
    double L2norm();
     

    //access private variables:
    int nptsx();
    int nptsy();
    double lengthx();
    double lengthy();
    double grid_sizex();
    double grid_sizey();

    //operators:
    double &operator()(int i, int j);
    field &operator=(const field& other_field);
    field operator-(const field& other_field);

private:
    double **u;
    double Lx, Ly;
    int nx, ny;
    double hx, hy;
};

#endif