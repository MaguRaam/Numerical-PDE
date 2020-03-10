#include "../include/helmholtz.h"



void convergence(int N){
    double L = 1.0;
    helmholtz problem(L,N);
    problem.Solve();   
}


int main(){

 int N = 200; 
 convergence(N);
 



    return 0;
}


