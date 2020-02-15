#include "field.h"

double rhs(const double& x,const double& y){
    int N=5,p=1,q=1;
    double norm = (pow(2.0/N,0.5))*(pow(2.0/N,0.5));
    return (norm*sin(p*M_PI*x)*sin(q*M_PI*y));
}



int main()
{
    //grid:
    double L=1.0;
    int N = 5;
    
    field b(L,N);
    return 0;
}