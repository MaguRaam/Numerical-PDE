
#include "../include/field.h"


 



//right handside 
double rhs(const double& x,const double& y){
   return (1000.0 -200.0*pow(M_PI,2))*sin(10.0*M_PI*x)*cos(10.0*M_PI*y);
}

 
int main(){

    //grid:
    double Lx = 1.0, Ly = 1.0;
    int nx = 100, ny = 100;
     
    //create and initialize b:
    field b(Lx,Ly,nx,ny);
    b.set_function(rhs);
    b.write_matplolib("plot/","b.dat"); 
    
    

    


    return 0;
}