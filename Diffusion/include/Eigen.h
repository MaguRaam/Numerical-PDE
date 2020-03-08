#ifndef EIGEN_H
#define EIGEN_H

#include "Field.h"
#include "Vector.h"
#include <cmath>
#include <cassert>

namespace Laplace
{
    Field EigenVector(const Grid& g);
    Field EigenVectorInverse(const Grid& g);
    Vector EigenValueX(const Grid &g);
} // namespace Laplace



#endif