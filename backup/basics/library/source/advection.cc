#include <iostream>
#include <cmath>
#include <vector>
#include "field.h"


using namespace std;


double analytic(double x,double y,double t){
    double c = 1.0;
    double xprime = x - c*t;
    double yprime = y - c*t;
    if (xprime>0.5 && yprime>0.5 && xprime<1.0 && yprime<1.0)
        return 2.0;
    else 
        return 1.0;
} 

 




int main(){

    //grid:
    int    nx = 10 ,ny = 10;
    double Lx = 1.0, Ly = 1.0;
    double hx = Lx/(nx-1),hy = Ly/(ny-1);


    //time:
    double Nt = 20;
    double dt = 0.01;

    //wave spped:
    double c =1.0;

    //initialize field:
    field u(Lx,Ly,nx,ny);
    u.set_function(analytic,0.0);

    //temporary field:
    field un(u);
    double t=0.0;

    for (int it=1;it<Nt;it++){
        t+=dt;
        un = u;
        for (int i=1;i<nx;i++)
            for (int j=1;j<ny;j++){
                u(i,j) = un(i,j) - (c*dt/hx*(un(i,j)-un(i-1,j))) - (c*dt/hy*(un(i,j)-un(i,j-1)));
            }
        
    }

    //exact solution computed at time t:
    field uexact(u);
    uexact.set_value(0.0);
    uexact.set_function(analytic,t);

    //compute error:
    field error(u);
    error.set_value(0.0);
    error = u - uexact;

    //L2 error:
    cout<<"L2 error = "<<error.L2norm()<<"\n";
    

    return 0;
}
