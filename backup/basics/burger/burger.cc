//2d Nonlinear advection:

#include <iostream>
#include <cmath>
#include <fstream>
#include <cassert>
#include <string>

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

//write data in tecplot:
void plot_tecplot(unsigned int i,const unsigned int digits,int nx,int ny,double* x,double* y,double** u) const{
    ofstream output;
    const string filename = "plot/u.dat";
    output.open(filename);
    output.flags( std::ios::dec | std::ios::scientific );
    output.precision(6);

    assert(output.is_open());

    output<< "TITLE = \"Linear Convection 2D\" " << "\n"
          << "VARIABLES = \"x\", \"y\", \"U\" " <<  "\n";
    output<< "Zone I = " << ny << " J = " << nx <<  "\n";
    
    for (unsigned int i=0; i < nx; i++) {
		    for (unsigned int j=0; j < ny; j++) {
                output << x[i] << "\t" << y[j] << "\t" << u[i][j] << "\n";
		    }
	    }
    output.close();

}



int main(){

    //grid:
    double Lx = 2.0, Ly = 2.0;
    int nx = 41, ny = 41;
    double dx = Lx/(nx-1), dy = Ly/(ny-1);
    double* x = linspace(0,Lx,nx);
    double* y = linspace(0,Ly,ny);

    //stability criteria:
    double sigma = .0009;

    //diffusion coefficient:
    double nu = 0.01; 

    //time:
    double dt = sigma * dx * dy / nu;
    double T = 0.027;

    //initialize field:

    //initialize u:
    double** u = create_field(ny,nx,1.0);
    double** un = create_field(ny,nx,0.0);
    initialize_field(u,dx,dy);

    //initialize v:
    double** v = create_field(ny,nx,1.0);
    double** vn = create_field(ny,nx,0.0);
    initialize_field(v,dx,dy);

    //start time loop:
    double t =0.0;
    int    iter=0;
    while (t<T){
        t+=dt;
        iter++;
        copy_field(un,u,ny,nx);
        copy_field(vn,v,ny,nx);
        for (int j=1;j<ny-1;j++){
            for (int i=1;i<nx-1;i++){
                u[j][i] = un[j][i] - (un[i][j]*dt/dx)*(un[j][i] - un[j][i-1]) - (vn[i][j]*dt/dy)*(un[j][i] - un[j-1][i])
                          + (nu*dt/(dx*dx))*(un[j][i+1] - 2*un[j][i] + un[j][i-1]) 
                          + (nu*dt/(dy*dy))*(un[j+1][i] - 2*un[j][i] + un[j-1][i]);


                v[j][i] = vn[j][i] - (un[i][j]*dt/dx)*(vn[j][i] - vn[j][i-1]) - (vn[i][j]*dt/dy)*(vn[j][i] - vn[j-1][i])
                          + (nu*dt/(dx*dx))*(vn[j][i+1] - 2*vn[j][i] + vn[j][i-1]) 
                          + (nu*dt/(dy*dy))*(vn[j+1][i] - 2*vn[j][i] + vn[j-1][i]);
            }                
        }
        if (iter%10 == 0){
            cout<<iter<<"  time = "<<t<<"\n";
            plot_tecplot(iter,3,nx,ny,x,y,u);
        }
        
    }

    //plot field 
    plot_field(u,x,y,ny,nx);


    //delete field:
    delete_field(u,ny);
    delete_field(v,ny);
    delete [] x;
    delete [] y;
    return 0;
}