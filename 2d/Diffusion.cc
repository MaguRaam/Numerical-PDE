//Solves Diffusion Equation using Matrix Diagonalization

//Schemes:
//Explicit Euler  alpha = 0
//Implicit Euler  alpha = 1
//Crank-Nicholson alpha = 0.5

#include <iostream>
#include <cmath>
#include "Eigen.h"
#include "Matrix.h"
#include "Grid.h"
#include "Time.h"
#include "Plot.h"

using std::cout;
using std::vector;
typedef Matrix<double> ScalarField;

double exact_solution(const double &x, const double &y, const double &t)
{
    return (exp(-8.0 * M_PI * M_PI * t) - ((cos(2.0 * M_PI * t) - 4.0 * M_PI * sin(2.0 * M_PI * t)) / ((2.0 * M_PI) + (32.0 * M_PI * M_PI * M_PI)))) * (sin(2.0 * M_PI * x) * sin(2.0 * M_PI * y));
}

double right_handside(const double &x, const double &y, const double &t)
{
    return sin(2.0 * M_PI * x) * sin(2.0 * M_PI * y) * sin(2.0 * M_PI * t);
}

int main()
{
    //space:
    const Grid G(128, 1.0);

    
    //time:
    const Time T(0.5, 0.00390625);

    ScalarField U(G.N() + 1, G.N() + 1, 0.0), Uexact(U), B(U), Bnew(U);
    ScalarField Utilde(U), Btilde(U);
    
    //Eigen Vector and Eigen Value of Laplace operator
    const Matrix<double> P = Laplace::EigenVector(G);
    const Matrix<double> Pinverse = P.transpose() * (2.0 * G.h()); //Normalize:
    const vector<double> Lambda = Laplace::EigenValue(G);

    double t = 0.0;
    const double value = pow(G.h(), 2.0) / (0.5 * T.dt());

    for (unsigned int i = 0; i < T.Nt(); i++)
    {
        cout << "Time = " << t << "\n";

        //Compute Btilde from B:
        SetFunction(B, right_handside, t, G);
        SetFunction(Bnew, right_handside, t + T.dt(), G);
        Btilde = Pinverse * (B + Bnew) * P;

        //Solve:
        for (unsigned int m = 1; m < G.N(); m++)
        {
            for (unsigned int n = 1; n < G.N(); n++)
            {
                Utilde(m, n) *= (value + Lambda[m] + Lambda[n]);
                Utilde(m, n) += pow(G.h(), 2.0) * Btilde(m, n);
                Utilde(m, n) /= (value - Lambda[m] - Lambda[n]);
            }
        }

        t += T.dt();
    }

    //Compute U from Utilde;
    U = P * Utilde * Pinverse;

    //Compute Uexact:
    SetFunction(Uexact, exact_solution, t, G);

    //Write data:
    WriteDataMatplot(Uexact, G, "../plot/", "Uexact.dat");
    WriteDataMatplot(U, G, "../plot/", "U.dat");
    WriteDataMatplot(B, G, "../plot/", "B.dat");

    //Error:
    cout << "L2 error = " << L2(Uexact, U, G) << "\n";
    cout << "Linfty error = " << Linfty(Uexact, U, G) << "\n";

    return 0;
}
