//2d diffusion:

#include <iostream>
#include <string>
#include <cmath>
#include <cassert>
#include <vector>
#include <fstream>
 

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

//initialize field:
void initialize_field(double** u,double dx,double dy){
    for (int j = int(0.5/dy);j<int(1.0/dy);j++){
        for (int i = int(0.5/dx);i<int(1.0/dx);i++)
            u[j][i] = 2.0;
    }

}

//writes data in matlplotlib format:
void plot_field(double** U,double* X,double* Y,int Ny,int Nx){
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

 

//copy field:
void copy_field(double** ucopy,double** uoriginal,int Ny,int Nx){
    for (int j=0;j<Ny;j++){
        for (int i=0;i<Nx;i++)
            ucopy[j][i] = uoriginal[j][i];
    }
}



//set boundary condition:
void set_boundary(double** u,int ny,double dy){
    for (int j=0;j<ny;j++){
        u[j][0] = 0.0;             //set u = 0      at  x = 0
        u[j][ny-1] = dy*j;          //set u = y      at  x = 2.0
        u[0][j] = u[1][j];         //set du/dy = 0  at  y = 0
        u[ny-1][j] = u[ny-2][j];   //set du/dy = 0  at  y = 2.0
    }
}

//solve laplace:
void laplace(double** u,double dx,double dy,int ny,int nx,int iter){
     
    double** un = create_field(ny,nx,0.0);

    //dummy time for plotting purpose:
    double t=0.0;

    for (int i=0;i<iter;i++){
        copy_field(un,u,ny,nx);
        for (int j=1;j<ny-1;j++){
            for (int i=1;i<nx-1;i++)
                u[j][i] = ( dy*dy*(un[j][i+1] + un[j][i-1])  +  dx*dx*(un[j+1][i] + un[j-1][i]) )/(2*(dx*dx + dy*dy) ); 
        }
        //set boundary:
        set_boundary(u,ny,dy);


    }

    delete_field(un,Ny);
}

int main(){

    //grid:
    double Lx = 2.0, Ly = 2.0;
    int nx = 31, ny = 31;
    double dx = Lx/(nx-1), dy = Ly/(ny-1);
    double* x = linspace(0,Lx,nx);
    double* y = linspace(0,Ly,ny);

    //initialize field:
    double** u = create_field(ny,nx,0.0);
    set_boundary(u,ny,dy);
     
    //solve laplace:
    laplace(u,dx,dy,ny,nx,100000);

    //plot:
    plot_field(u,x,y,ny,nx);

    //delete memory:
    delete_field(u,ny);
    delete [] x;
    delete [] y;  

    return 0;
}