#include "../include/helmholtz.h"



void convergence(int N){
    double L = 1.0;
    helmholtz problem(L,N);
    problem.Solve();  
}


int main(){

 int N = 16;
 for (int i = 1;i<6;i++){
     N = 2*N;
     convergence(N);
 }   



    return 0;
}


