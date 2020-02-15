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
     
 }