
#include <iostream>
#include <cmath>
#include <fstream>
#include <cassert>
#include <string>

using namespace std;

//Auxillary functions:

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
void set_b(double** u,double dx,double dy,int Nx,int Ny){
    int p=25,q=25;
    double norm = (pow(2.0/Nx,0.5))*(pow(2.0/Ny,0.5));
     for (int i=0;i<=Nx;i++){
         for (int j=0;j<=Ny;j++)
            u[j][i] = (norm*sin(p*M_PI*i*dx)*sin(q*M_PI*j*dy));
     }
}

//writes data in matlplotlib format:
void plot_field(double** U,double dx,double dy,int Ny,int Nx,string filename){
    ofstream x,y,u;

    //open files:
    x.open("plot/x.dat");
    y.open("plot/y.dat");
    u.open(filename);

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
            x << i*dx <<"\n"; 
        }

    for (int j=0; j<Ny; j++)  {
            y<< j*dy <<"\n"; 
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

//projects given field on given eigen vector: 
double project_on_fourier_basis(double** b,int p,int q,int nx,int ny,double hx,double hy){
    double value = 0.0;

    double norm = (pow(2.0/nx,0.5))*(pow(2.0/ny,0.5));
    for (int j=1;j<=ny-1;j++){
        for (int i=1;i<=nx-1;i++)
            value += (norm*sin(p*M_PI*i*hx)*sin(q*M_PI*j*hy)*b[j][i]);
    }
     
    return value; 
}

//eigenvalue:
double eigen_value(int q,int p,int nx,int ny,double hx,double hy){
    double value = 0.0;
    value = -( (4.0/hx*hx)*(pow(sin(M_PI*q/2*nx),2))    ) - ( (4.0/hy*hy)*(pow(sin(M_PI*p/2*ny),2)) );
    return value;   
}



//Important functions:

//STEP-1: b~ = P'b
void fourier_tranform(double** b_tilde,double** b,int ny,int nx,double hx,double hy){
    for (int q=1;q<=ny-1;q++){
        for (int p=1;p<=nx-1;p++)
            b_tilde[q][p] =  project_on_fourier_basis(b,p,q,nx,ny,hx,hy);
    }   
}

//STEP-2:  xi~ = bi~/lambdai
void solve(double** x_tilde,double** b_tilde,int nx,int ny,double hx,double hy){
    for (int q=1;q<=ny-1;q++){
        for (int p=1;p<=nx-1;p++)
            x_tilde[q][p] =  b_tilde[q][p]/eigen_value(p,q,nx,ny,hx,hy);
    }
}

//STEP-3:  x = Px~
void inverse_fourier_transform(double** x,double** x_tilde,int nx,int ny,double hx,double hy){
    
    double norm = (pow(2.0/nx,0.5))*(pow(2.0/ny,0.5));
    //Loop over grid points
    for (int j=1;j<=ny-1;j++){
        for (int i=1;i<=nx-1;i++){

            //Loop over all eigen states
            for (int q=1;q<=ny-1;q++){
                for (int p=1;p<=nx-1;p++){
                    x[j][i] += x_tilde[q][p]*( norm*sin(p*M_PI*i*hx)*sin(q*M_PI*j*hy) );

                }
             
            }  


        }
            
    }
}


int main(){
    //grid:
    double Lx = 1.0, Ly = 1.0;
    int nx = 50, ny = 50;
    double dx = Lx/nx, dy = Ly/ny;

    //create and initialize rhs:
    double** b = create_field(ny+1,nx+1,0.0);
    set_b(b,dx,dy,nx,ny);
    plot_field(b,dx,dy,ny,nx,"plot/b.dat"); 


    //STEP-1: Fourier transform:
    double** b_tilde = create_field(ny+1,nx+1,0.0);
    fourier_tranform(b_tilde,b,ny,nx,dx,dy); 
    plot_field(b_tilde,dx,dy,ny,nx,"plot/b_tilde.dat"); 


    //STEP-2: Solve:
    double** x_tilde = create_field(ny+1,nx+1,0.0);
    solve(x_tilde,b_tilde,nx,ny,dx,dy);
    plot_field(x_tilde,dx,dy,ny,nx,"plot/u_tilde.dat"); 


    /*//STEP-3: Inverse Fourier transform:
    double** x = create_field(ny+1,nx+1,0.0);
    inverse_fourier_transform(x,x_tilde,nx,ny,dx,dy);
    plot_field(x,dx,dy,ny,nx,"plot/u.dat");*/ 





    //delete field:
    delete_field(b,ny); 
    delete_field(b_tilde,ny); 
    delete_field(x_tilde,ny);
    //delete_field(x,ny);
    return 0;
}
