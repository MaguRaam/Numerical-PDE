//2d poisson:

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

//copy field:
void copy_field(double** ucopy,double** uoriginal,int Ny,int Nx){
    for (int j=0;j<Ny;j++){
        for (int i=0;i<Nx;i++)
            ucopy[j][i] = uoriginal[j][i];
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

std::string int_to_string (unsigned int value, const unsigned int digits) {
    std::string lc_string = std::to_string(value);  
    
    if (lc_string.size() < digits) {
        // We have to add the padding zeroes in front of the number
        const unsigned int padding_position = (lc_string[0] == '-')
                                                ?
                                                1
                                                :
                                                0;

        const std::string padding(digits - lc_string.size(), '0');
        lc_string.insert(padding_position, padding);
    }
         
    return lc_string;
}




//plot vtk:
void plot_vtk(double time, unsigned int i, const unsigned int digits,double** U,double* X,double* Y,int Ny,int Nx){
    std::ofstream vtk;
    const std::string filename = "vtk/plot_" + int_to_string (i, digits) + ".vtk";
    vtk.open (filename);
    vtk.flags( std::ios::dec | std::ios::scientific );
    vtk.precision(6);

    assert(vtk.is_open());
    vtk << "# vtk DataFile Version 2.0" << "\n"; 
    vtk << "Poisson" << "\n";
    vtk << "ASCII" << "\n";
    vtk << "\nDATASET STRUCTURED_GRID" << "\n"; 
    vtk << "\nFIELD FieldData 1" << "\n"; 
    vtk << "TIME 1 1 double" << "\n"; 
    vtk << time << "\n"; 
    vtk << "\nDIMENSIONS " <<  Ny << " " << Nx << " " << 1 << "\n"; 
    vtk << "POINTS " << Nx*Ny <<  " double" << "\n";
    vtk << "\n";

    for (unsigned int j = 0; j < Ny; j++) {
        for (unsigned int i = 0; i < Nx; i++) {
            vtk << Y[j] << " " << X[i] << " " << 0.0 << "\n"; 
        }
    }
    vtk << "\nPOINT_DATA " << Nx*Ny << "\n";
    vtk << "\nSCALARS U double 1" << "\n"; 
    vtk << "LOOKUP_TABLE default" << "\n"; 
    vtk << "\n"; 

    for (unsigned int j = 0; j < Ny; j++) {
        for (unsigned int i = 0; i < Nx; i++) {
            vtk << U[j][i] << "\n"; 
        }
    }

   vtk.close();  
    
}

//solve poisson equation:
void poisson(double** p,double** b,double dx,double dy,int ny,int nx,double* x,double* y,int iter){
    double** pn = create_field(ny,nx,0.0);

    double t=0.0,dt = 0.01;

    for (int i=0;i<iter;i++){
        t+=dt;
        copy_field(pn,p,ny,nx);
        for (int j=1;j<ny-1;j++){
            for (int i=1;i<nx-1;i++)
                 p[j][i] = (  (pn[j][i+1] + pn[j][i-1])*dy*dy + (pn[j+1][i] + pn[j-1][i])*dx*dx - (b[j][i]*dx*dx*dy*dy) )/(2* (dx*dx + dy*dy));
        }

        
        plot_vtk(t,iter,1,p,x,y,ny,nx);
    }

    delete_field(pn,ny);
}




int main(){

    //grid:
    double Lx = 2.0, Ly = 1.0;
    int nx = 50, ny = 50;
    double dx = Lx/(nx-1), dy = Ly/(ny-1);
    double* x = linspace(0,Lx,nx);
    double* y = linspace(0,Ly,ny);    

    //initialize pressure:
    double** p = create_field(ny,nx,0.0);

    //initialize right handside:
    double** b = create_field(ny,nx,0.0);
    b[int(ny/4)][int(nx/4)] = 100;
    b[int(3*ny/4)][int(3*nx/4)] = -100;

    //solve poisson equation:
    poisson(p,b,dx,dy,ny,nx,x,y,100);

    //plot:
    plot_field(p,x,y,ny,nx);

    //delete memory:
    delete_field(p,ny);
    delete [] x;
    delete [] y;  

    return 0;
}