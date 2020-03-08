#include <iostream>
#include <cmath>
#include <chrono>
#include "../include/Grid.h"
#include "../include/Eigen.h"
#include "../include/Field.h"

double ExactSolution(double x, double y, double t)
{
    return sin(10.0 * M_PI * x) * sin(10.0 * M_PI * y) * t;
}

double RightHandSide(double x, double y, double t)
{
    double value = -200.0 * M_PI * M_PI * sin(10.0 * M_PI * x) * sin(10.0 * M_PI * y);
    value += 1000.0 * (sin(10.0 * M_PI * x) * sin(10.0 * M_PI * y));
    return value * t;
}

int main()
{

    const Grid g(20, 20, 1.0, 1.0);
    const Field P = Laplace::EigenVector(g);
    const Field PInv = Laplace::EigenVectorInverse(g);
    const Vector LambdaX = Laplace::EigenValueX(g);

    //Exact Solution:
    Field uexact(20, 20, 1.0, 1.0);
    uexact.SetFunction(ExactSolution, 1.0);

    //Right Hand Side:
    Field F(uexact);
    F.SetFunction(RightHandSide, 1.0);
    //TODO add h2


    //Compute Ftilde = P^{-1} F Q^{-T}
    Field Temp   = PInv*F;
    Field Ftilde = Temp*P;

    //Solve:
    

     

    return 0;
}