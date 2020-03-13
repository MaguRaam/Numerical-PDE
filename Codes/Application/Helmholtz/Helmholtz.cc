#include <iostream>
#include <cmath>
#include "Eigen.h"
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
    const Grid g(512, 1.0);
    const double h = g.h();

    //Exact Solution
    const ScalarField Uexact = SetFunction(exact_solution, g);

    //Right Hand side
    const ScalarField F = SetFunction(right_handside, g) * h * h;

    //Eigen Vector of Laplace operator
    const Matrix<double> P = Laplace::EigenVector(g);
    const Matrix<double> Pinverse = P.transpose() * (2.0 * h); //Normalize:

    //Eigen value of Laplace Operator
    const std::vector<double> Lambda = Laplace::EigenValue(g);

    //Solve Linear system using Matrix Diagonalization

    //Compute P^{-1} F Q^{-T}
    ScalarField Ftilde = Pinverse * F * P;

    //Compute Utilde
    ScalarField &Utilde = Ftilde;
    for (unsigned int m = 1; m < g.N(); m++)
    {
        for (unsigned int n = 1; n < g.N(); n++)
            Utilde(m, n) /= (Lambda[m] + Lambda[n] +  1000 * h * h);
    }

    //Compute U = P Utilde Q^T
    ScalarField U = P * Utilde * Pinverse;

    std::cout << "L2error = " << L2(Uexact, U, g) << "\n";

    return 0;
}
