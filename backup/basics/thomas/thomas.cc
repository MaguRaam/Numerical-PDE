#include <iostream>
#include <cmath>
#include <fstream>
//solves laplace equation in 1d:

using namespace std;

void Thomas(int N,double* b,double* a,double* c,double* x,double* q){
    
    //lower,upper and diagonal
    double* l = new double[N];
    double* u = new double[N];
    double* d = new double[N];

    double* y = new double[N];

    //Step 1 : LU decomposition A = LU:
    d[0] = a[0]; u[0] = c[0];

    for (int i =1;i<N-2;i++){
        l[i] = b[i]/d[i];
        d[i+1] = a[i+1] - l[i]*u[i];
        u[i+1] = c[i+1]; 

    }

    l[N-2] = b[N-2]/d[N-2];
    d[N-1] = a[N-1] - l[N-2]*u[N-2];

    //Step 2 : Forwad substitution Ly = q:
    y[0] = q[0];
    for (int i=1;i<N;i++)
        y[i] = q[i] - l[i-1]*y[i-1];

    //Step 3 : Backward substitution Ux = y:
    x[N-1] = y[N-1]/d[N-1];
    for(int i=N-2;i>=0;i--)
        x[i] = (y[i] - u[i]*x[i+1])/d[i];


    delete[] l;
    delete[] u;
    delete[] d;
    delete[] y;
}

double write_data(double* x,double h,int N){
    ofstream output("result.dat");
    for (int i=0;i<N;i++)
        output<<i*h<<"\t"<<x[i]<<"\n";
    output.close();
}


int main(){
    
     

    return 0;
}