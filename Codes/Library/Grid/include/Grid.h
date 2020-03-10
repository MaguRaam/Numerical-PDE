#ifndef GRID_H
#define GRID_H

#include <vector>
#include "Matrix.h"

//2D Square Grid:

class Grid
{
private:
    unsigned int N_;
    double L_,h_;
    std::vector<double> x_;

public:
    Grid(unsigned int N,double L);
    const double& x(const unsigned int& i) const;
    const double& y(const unsigned int& j) const;
    const unsigned int&    N() const;
    const double& h() const;
    const double& L() const;
    friend Matrix<double> SetFunction(double (*func)(const double& x,const double& y),const Grid& g);
 
};


#endif
