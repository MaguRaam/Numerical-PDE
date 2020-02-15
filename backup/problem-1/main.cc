


//Solves tridiagonal system that arise from finite diffrencing of Laplacian using SVD:
// See https://en.wikipedia.org/wiki/Eigenvalues_and_eigenvectors_of_the_second_derivative for eigevectors of Discrete Laplacian:

#include <iostream>
#include <cmath>
#include <vector>
#include <cassert>
#include <fstream>
#include <iomanip>
#include <string>

using namespace std;



void inititalize_b(vector<double>& b,int n){
    for(int i=0;i<=n;i++)
        b[i]=-sin(2*M_PI*i/n);
}
 
void fourier_transform(const vector<double>& b,vector<double>& b_tilde,int n){
    //Loop over modes:
    double norm = pow(2.0/n,0.5);
    for (int j=1;j<=n-1;j++){
        for (int i=1;i<=n-1;i++){
            b_tilde[j] += norm*sin(i*j*M_PI/n)*b[i];
        }
    }
}

void solve(vector<double>& x_tilde,const vector<double>& b_tilde,int n,double h){
    double temp = -4.0/(h*h);
    for(int j=1;j<=n-1;j++)
        x_tilde[j] = b_tilde[j]/(temp*pow(sin(M_PI*j/(2*n)),2));

}

void inverse_fourier_transform(vector<double>& x,const vector<double>& x_tilde,int n){
    double norm = pow(2.0/n,0.5);
    //Loop over grid points:
    for(int i=1;i<=n-1;i++){
        for (int j=1;j<=n-1;j++){
            x[i] +=norm*sin(i*j*M_PI/n)*x_tilde[j];
        }
    }
}

double L2error(const vector<double>& x,double h,int n){
    double error=0.0;
    for (int i = 0; i < x.size(); i++)
        error += pow( (x[i] - sin(2.0*M_PI*i/n))*h ,2);

    return pow(error,0.5);
}

double write_data(const vector<double>& u,double h,string filename){
    ofstream output(filename);
    output << setprecision(10) << setiosflags(ios::scientific);
    for (int i=0;i<u.size();i++)
        output<<u[i]<<"\n";
    output.close();
}

void Laplace(int n){
    //grid:
    double L = 1.0;
    double h = L/n;

    //initialize solution vector:
    vector<double> x(n+1,0.0);

    //inititalize b vector
    vector<double> b(n+1,0.0);
    inititalize_b(b,n);
    write_data(b,h,"plot/b.dat");
    
    //fourier transform:
    vector<double> b_tilde(n+1,0.0);
    fourier_transform(b,b_tilde,n);
    write_data(x,h,"plot/u.dat");

    //solve for x_tilde in fourier domain:
    vector<double> x_tilde(n+1,0.0);
    solve(x_tilde,b_tilde,n,h);

    //inverse fourier transform of the solution:
    inverse_fourier_transform(x,x_tilde,n);

    //compute error:
    cout<< setprecision(10) << setiosflags(ios::scientific);
    cout<<n<<"\t\t\t"<<L2error(x,h,n)<<"\n";

    //output:
    write_data(x,h,"plot/u.dat");

}

int main(){
    int n=50;  
    Laplace(n);
    

    return 0;
}
