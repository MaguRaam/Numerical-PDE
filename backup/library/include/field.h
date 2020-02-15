//field object:

// stores the data in a uniform N×N computational grid.

#ifndef FIELD_H
#define FIELD_H

#include <iostream>
#include <fstream>
#include <cassert>
#include <string>
#include <cmath>

using namespace std;

class field
{
public:
    field(double L, int N);
    ~field();
    void write_output(string folderpath, string filename);
    void set_function(double (*func)(const double &, const double &));
    double length();
    double grid_size();
    int npts();
    double L2norm()const;
    double read(int j,int i)const;
    
    

    double& operator()(int j,int i);

private:
    double L, h;
    int N;
    double **u;
};

#endif