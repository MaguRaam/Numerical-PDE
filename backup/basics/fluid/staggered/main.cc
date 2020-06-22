
//Simple fluid code for droplet falling:



#include <iostream>
#include <cmath>
#include <fstream>
#include <cassert>
#include <vector>

using namespace std;

void linspace(vector<double>& x,double dx,double Lx){
    
    //Interior cells:
    int nx = x.size()-2;
    for (int i = 1; i <= nx; i++)
        x[i] = (i-0.5)*dx;

    //Ghost cells:
    x[nx+1] =  Lx + (0.5*dx);
    x[0]    = -0.5*dx;
} 

//allocate memory for the fields:
double** create_field(int Ny,int Nx,double value){
   //create memory:
    double** u = new double*[Ny];
    for (int j=0;j<Ny;j++)
        u[j] = new double[Nx];
    
    for (int j=0;j<Ny;j++){
        for (int i=0;i<Nx;i++)
            u[j][i] = value;
    }
    
    return u;
        
} 

//delete memory:
void delete_field(double** u,int Ny){
    for (int j=0;j<Ny;j++)
        delete [] u[j];
    delete [] u;
}

//set droplet density:
void set_density(double** rho,int nx,int ny,const vector<double>& x,const vector<double>& y,double rho2){

    //droplet geometry and location:
    double rad = 0.15, xc = 0.5, yc = 0.7;

    for (int j = 1; j <= ny; j++)
    {
        for (int i = 1; i <= nx; i++)
        {
            if (( pow(x[i] - xc,2) + pow(y[j] - yc,2)   )<rad*rad){
                rho[j][i] = rho2;
            }
        }
        
    }
    

}

//writes data in matlplotlib format:
void plot_field(double** U,const vector<double>& X,const vector<double>& Y,int Ny,int Nx){
    ofstream x,y,u;

    //open files:
    x.open("plot/x.dat");
    y.open("plot/y.dat");
    u.open("plot/u.dat");

    //assert file is open:
    assert(x.is_open());
    assert(y.is_open());
    assert(u.is_open());

    //set control on output data
    x.setf(std::ios::scientific);
    x.setf(std::ios::showpos);
    x.precision(13);

    y.setf(std::ios::scientific);
    y.setf(std::ios::showpos);
    y.precision(13);

    u.setf(std::ios::scientific);
    u.setf(std::ios::showpos);
    u.precision(13);

    for (int i=0; i<Nx; i++)  {
            x << X[i] <<"\n"; 
        }

    for (int j=0; j<Ny; j++)  {
            y<< Y[j] <<"\n"; 
    }
    for (int j=0; j<Ny; j++) {
        for (int i=0; i<Nx; i++) {
            u << U[j][i] << " " ;
        }
        
        u<<"\n"; 
    }

    //close fiels:
    x.close();y.close();u.close();

}



int main(int argc, char const *argv[])
{
    //Interior grid:
    double Lx = 1.0,   Ly = 1.0;
    int    nx = 320,    ny = 320;
    double dx = Lx/nx, dy = Ly/ny;

    //x and y coord (Interior + Ghost):
    vector<double> x(nx+2,0.0);
    linspace(x,dx,Lx);

    vector<double> y(ny+2,0.0);
    linspace(y,dy,Ly);

    //Time:
    double dt = 0.00125;
    double t  = 0.0;
    int    Nt = 100;
    
    //set density and viscosity (fluid + droplet):
    double rho1 = 1.0, rho2 = 2.0, mu = 0.01;
    
    //set droplet density:
    double** rho = create_field(ny+2,nx+2,rho1);
    set_density(rho,nx,ny,x,y,rho2);

    //plot density field:
    plot_field(rho,x,y,ny+2,nx+2);


    //Inititalize fields:
    double** u = create_field(ny+2,nx+1,0.0);
    double** v = create_field(ny+1,nx+2,0.0);
    double** p = create_field(ny+2,nx+2,0.0);
    
    //start time loop:    
    for (int it = 0; it < Nt; it++)
    {
        t+=dt;
        cout<<"time = "<<t<<"\n";
    }
    

     
    //delete field:
    delete_field(u,ny);
    delete_field(v,ny);
    delete_field(p,ny);
    delete_field(rho,ny);


    return 0;
}
