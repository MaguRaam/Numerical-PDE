#include <iostream>
#include <cmath>
#include "EigenValueDecompose.h"
#include "Matrix.h"
#include "Grid.h"
#include "Time.h"
#include "Plot.h"

double exact_solution(const double &x, const double &y)
{
    return sin(10.0 * M_PI * x) * sin(10.0 * M_PI * y);
}

double right_handside(const double &x, const double &y)
{
    double value = -200.0 * M_PI * M_PI * sin(10.0 * M_PI * x) * sin(10.0 * M_PI * y);
    value += 1000.0 * (sin(10.0 * M_PI * x) * sin(10.0 * M_PI * y));
    return value;
}

int main()
{

    Grid g(100, 1.0);
    Matrix<double> uexact = SetFunction(exact_solution, g);
    Matrix<double> b = SetFunction(right_handside, g) * (g.h() * g.h());

    Matrix<double> P = Laplace::EigenVector(g);
    Matrix<double> Pinverse = P.transpose() * (2.0 * g.h());
    std::vector<double> Lambda = Laplace::EigenValue(g);

    Matrix<double> temp = Pinverse * b;
    Matrix<double> btilde = temp * P;

    btilde = MatrixDiagonalization::Solve(btilde, Lambda, 1000 * g.h() * g.h());

    temp = P * btilde;
    Matrix<double> u = temp * Pinverse;

    std::cout << "L2error = " << L2(uexact, u, g) << "\n";

    return 0;
}
