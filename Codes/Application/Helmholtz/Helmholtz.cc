#include <iostream>
#include <cmath>
#include "EigenValueDecompose.h"
#include "Matrix.h"
#include "Grid.h"
#include "Plot.h"

typedef Matrix<double> ScalarField;

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
    //Grid
    Grid g(200, 1.0);

    //Exact Solution
    ScalarField uexact = SetFunction(exact_solution, g);

    //Right Hand side
    ScalarField b = SetFunction(right_handside, g) * (g.h() * g.h());

    //Eigen Vector of Laplace operator
    Matrix<double> P = Laplace::EigenVector(g);
    Matrix<double> Pinverse = P.transpose() * (2.0 * g.h());

    //Eigen value of Laplace Operator
    std::vector<double> Lambda = Laplace::EigenValue(g);

    //Solve Linear system using Matrix Diagonalization

    //Compute P^{-1} F Q^{-T}
    ScalarField btilde = Pinverse * b * P;

    //Compute Utilde
    btilde = MatrixDiagonalization::Solve(btilde, Lambda, 1000 * g.h() * g.h());

    //Compute U = P Utilde Q^T
    ScalarField u = P * btilde * Pinverse;

    std::cout << "L2error = " << L2(uexact, u, g) << "\n";

    return 0;
}
