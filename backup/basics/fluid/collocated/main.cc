//Lid driven cavity:

#include <iostream>
#include <cmath>
#include <fstream>
#include <cassert>
#include <vector>

using namespace std;

//create x and y array:
double* linspace(double xo,double xl,int N){
    double* x = new double[N];
    double h = (xl - xo)/(N-1);

    for (int i=0;i<N;i++)
        x[i] = xo + i*h;
    return x;
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

//copy field:
void copy_field(double** ucopy,double** uoriginal,int Ny,int Nx){
    for (int j=0;j<Ny;j++){
        for (int i=0;i<Nx;i++)
            ucopy[j][i] = uoriginal[j][i];
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

//solve pressure poisson equation:
void poisson(double** p,double** pn, double** b,double dx,double dy,int ny,int nx,double* x,double* y,int iter){
    
    for (int i=0;i<iter;i++){
        copy_field(pn,p,ny,nx);
        for (int j=1;j<ny-1;j++){
            for (int i=1;i<nx-1;i++)
                 p[j][i] = (  (pn[j][i+1] + pn[j][i-1])*dy*dy + (pn[j+1][i] + pn[j-1][i])*dx*dx - (b[j][i]*dx*dx*dy*dy) )/(2* (dx*dx + dy*dy));
        }

    }

 
}








int main(int argc, char const *argv[])
{
    //grid:
    double Lx = 2.0,   Ly = 2.0;
    int    nx = 41,    ny = 41;
    double dx = Lx/(nx-1), dy = Ly/(ny-1);

    double* x = linspace(0,Lx,nx);
    double* y = linspace(0,Ly,ny); 
   

    //Time:
    double dt = 0.001;
    double t  = 0.0;
    int    nt = 500;
    int    nit = 50;

    //Lid velocity:
    double Ulid = 1.0;

    //fluid property:
    double rho = 1.0;
    double nu  = 0.1;

    //Reynolds and cfl:
    double Re =  Ulid*Lx/nu;
    double cfl = Ulid*dt/dx;

    cout<<"\n\nReynold's number = "<<Re<<"\n";
    cout<<"Cfl number = "<<cfl<<"\n\n";
    

    //Inititalize fields:
    double** u = create_field(ny,nx,0.0);
    double** v = create_field(ny,nx,0.0);
    double** p = create_field(ny,nx,0.0);
    double** b = create_field(ny,nx,0.0);

    //temporary fields:
    double** pn = create_field(ny,nx,0.0);

    cout<<"Start time loop"<<"\n\n";
    
    //start time loop:    
    for (int it = 0; it < nt; it++)
    {
        t+=dt;
        

        if (it%10 == 0)
            cout<<"time = "<<t<<"\n";
    }
    

     
    //delete field:
    delete_field(u,ny);
    delete_field(v,ny);
    delete_field(p,ny);
    delete_field(b,ny);

    //delete temporary field:
    delete_field(pn,ny);


    return 0;
}
