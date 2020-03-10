#ifndef EIGEN_H
#define EIGEN_H

#include "Matrix.h"
#include "Grid.h"
#include <vector>
#include <cmath>

namespace Laplace
{
Matrix<double> EigenVector(const Grid &g);
std::vector<double> EigenValue(const Grid &g);

} // namespace Laplace

namespace MatrixDiagonalization
{
Matrix<double> &Solve(Matrix<double> &btilde, const std::vector<double> &eigenvalue, const double &value);
}
#endif