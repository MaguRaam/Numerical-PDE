#include "../include/poisson.h"


poisson::poisson(double l, int n):L(l),N(n),h(L/N){
    b       =   new field(L,N);
    u       =   new field(L,N);
    un       =   new field(L,N);
   
}

poisson::~poisson(){
    delete b;
    delete u;
    delete un;
  
}

 void poisson::initialize_b(int n){
     
     for (int j=1;j<N;j++){
         for (int i=1;i<N;i++)
            (*b)(j,i) = -2.0*n*n*M_PI*M_PI*sin(n*M_PI*i*h)*sin(n*M_PI*j*h);
     }
     b->write_output("../plot/","b.dat");
 }



void poisson::exact_solution(int n){
    for (int j=1;j<N;j++){
        for (int i=1;i<N;i++)
            (*un)(j,i) = sin(n*M_PI*i*h)*sin(n*M_PI*j*h);
    }
    un->write_output("../plot/","uexact.dat");
}
 

  


 int poisson::solve(){
     int iter = 0;
     
     while (iter <= 1000){
         iter++;
         (*un) = (*u);
         for (int j=1;j<N;j++){
             for (int i=1;i<N;i++)
                (*u)(j,i) = 0.25*( (-(b->read(j,i))*h*h)  + ( un->read(j,i+1) + un->read(j,i-1) + un->read(j+1,i) + un->read(j-1,i)   )   );
         }

     }
     u->write_output("../plot/","u.dat");
     
     return iter;
 }