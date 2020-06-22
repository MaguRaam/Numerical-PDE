#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>

using namespace std;

double func(double x,double t,double nu){
    double value=0.0;
    value = -2*nu*(-(-8*t + 2*x)*exp(-pow((-4*t + x),2)/(4*nu*(t + 1)))/(4*nu*(t + 1)) - 
    (-8*t + 2*x - 4*M_PI)*exp(-pow((-4*t + x - 2*M_PI),2)/(4*nu*(t + 1)))/(4*nu*(t + 1)))/(exp(-pow((-4*t + x - 2*M_PI),2)/(4*nu*(t + 1))) 
    + exp(-pow((-4*t + x),2)/(4*nu*(t + 1)))) + 4;
    return value;
}

double u_exact(vector<double>& u,double t,double h,double nu){
    double x=0.0;
    for (int i=0;i<u.size();i++){
        x = i*h;
        u[i] = func(x,t,nu);
    }
}

double write_data(const vector<double>& u,double h){
    ofstream output("result.dat");
    for (int i=0;i<u.size();i++)
        output<<i*h<<"\t"<<u[i]<<"\n";
    output.close();
}

double norm(const vector<double>& u,double p){
    double value = 0.0;
    for(int i=0;i<u.size();i++)
        value+=pow(u[i],p);
    return pow(value,1.0/p);
}

void advection_diffusion(int N){
    //grid:
     double L = 2*M_PI;
     double h = L/(N-1);

     //stability:
     double nu = 0.07;

     //time:
     double dt = nu*h;
     int    Nt = 100;
     double t = 0.0;

     //initialize:
     vector<double> u(N,0.0);
     vector<double> un(N,0.0);
     u_exact(u,t,h,nu);

     //start time loop:
     for (int i=0;i<Nt;i++){
         t+=dt;
         un = u;
         for (int j=1;j<N-1;j++){
             u[j] = un[j] - un[j]*(dt/h)*(un[j]-un[j-1]) + (nu*dt/h*h)*(un[j+1] - 2 * un[j] + un[j-1]);           
         }
         u[0] = un[0] - un[0]*(dt/h)*(un[0]-un[N-1]) + (nu*dt/h*h)*(un[1] - 2 * un[0] + un[N-1]);
         u[N-1] = un[N-1] - un[N-1]*(dt/h)*(un[N-1]-un[N-2]) + (nu*dt/h*h)*(un[N-2] - 2 * un[N-1] + un[0]);

     }

     //compute exact solution at time t:
     vector<double> u_analytic(N,0.0);
     u_exact(u_analytic,t,h,nu);

     //error:
     vector<double> error(N,0.0);
     for (int i=0;i<N;i++)
        error[i] = (u_analytic[i] - u[i])*h;
     
     cout<<N<<"\t\t"<<norm(error,2.0)<<"\n";

     //write data:
     write_data(u,h);
}


int main(){

    int Ngrid = 8;
    for (int iter=1;iter<10;iter++){
        Ngrid*=2;
        advection_diffusion(Ngrid);
    }
      
    return 0;
}