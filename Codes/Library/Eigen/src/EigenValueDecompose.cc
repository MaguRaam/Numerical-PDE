#include "../include/EigenValueDecompose.h"

namespace Laplace
{
Matrix<double> EigenVector(const Grid &g)
{
    Matrix<double> eigenvector(g.N() + 1, g.N() + 1, 0.0);
    for (unsigned int i = 1; i < g.N(); i++)
    {
        for (unsigned int j = 1; j < g.N(); j++)
            eigenvector(i, j) = sin(i * M_PI * g.x(j));
    }
    return eigenvector;
}

std::vector<double> EigenValue(const Grid &g)
{
    std::vector<double> eigenvalue(g.N() + 1, 0.0);
    for (unsigned int i = 1; i < g.N(); i++)
        eigenvalue[i] = -4.0 * sin(M_PI * 0.5 * g.x(i)) * sin(M_PI * 0.5 * g.x(i));

    return eigenvalue;
}

} // namespace Laplace

namespace MatrixDiagonalization
{
Matrix<double> &Solve(Matrix<double> &btilde, const std::vector<double> &eigenvalue, const double &value)
{
    for (unsigned int i = 1; i < eigenvalue.size() - 1; i++)
    {
        for (unsigned int j = 1; j < eigenvalue.size() - 1; j++)
            btilde(i, j) /= (eigenvalue[i] + eigenvalue[j] + value);
    }
    return btilde;
}
} // namespace MatrixDiagonalization
