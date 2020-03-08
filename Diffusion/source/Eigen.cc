#include "../include/Eigen.h"

using namespace Laplace;

Field Laplace::EigenVector(const Grid &g)
{

    assert(g.nx() == g.ny());
    Field P(g.nx(), g.ny(),g.lx(),g.ly());
    for (int i = 1; i < g.nx(); i++)
    {
        for (int j = 1; j < g.ny(); j++)
        {
            P(i, j) = sin(i * M_PI * g.X(j));
        }
    }

    return P;
}

Field Laplace::EigenVectorInverse(const Grid &g)
{

    assert(g.nx() == g.ny());
    Field PInverse(g.nx(), g.ny(),g.lx(),g.ly());
    for (int i = 1; i < g.nx(); i++)
    {
        for (int j = 1; j < g.ny(); j++)
        {
            PInverse(i, j) = 2.0 * g.dx() * sin(j * M_PI * g.X(i));
        }
    }

    return PInverse;
}

Vector Laplace::EigenValueX(const Grid &g){
    Vector eigenvalue(g.nx());
    for (int i=1;i<g.nx();i++)
        eigenvalue(i) = -4.0*sin(M_PI*0.5*g.X(i))*sin(M_PI*0.5*g.X(i));
    return eigenvalue;
}

 
