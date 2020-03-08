#ifndef GRID_H
#define GRID_H

#include "Vector.h"
#include <cassert>

class Grid
{
protected:
    int Nx, Ny;
    double Lx, Ly, hx, hy;
    Vector *x;
    Vector *y;

public:
    Grid(int Nx, int Ny, double Lx, double Ly);
    Grid(const Grid& OtherGrid);
    ~Grid();
    double X(int i) const;
    double Y(int j) const;
    int  nx() const;
    int  ny() const;
    double dx() const;
    double dy() const; 
    double lx() const;
    double ly() const;  
};

#endif