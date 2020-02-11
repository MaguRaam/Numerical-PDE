
#include "field.h"

//right handside 
double rhs(const double& x,const double& y){
   int n = 2;
   return -2.0*n*n*M_PI*M_PI*sin(n*M_PI*x)*sin(n*M_PI*y);
}


class poisson{
    public:
        poisson(double Lx,double Ly,int nx,int ny);
        ~poisson();
        void initialize_b();
        void solve();

    private:
        field* u;
        field* un;
        field* b;
        void jacobi();
          
};

poisson::poisson(double Lx,double Ly,int nx,int ny){
    u  = new field(Lx,Ly,nx,ny);
    un = new field(*u);
    b  = new field(*u);
  


      
}

void poisson::initialize_b(){
    b->set_function(rhs);
    b->write_matplolib("plot/","b.dat");
}

void poisson::jacobi(){
    for (int j=1;j<u->ny();j++){
        for (int i=1;i<u->nx();i++)
            (*u)(j,i) = 0.25*((*un)(j,i+1) + (*un)(j,i-1) + (*un)(j+1,i) + (*un)(j-1,i) - (*b)(j,i)*(pow(u->hx(),2.0)));
    }
}

void poisson::solve(){
    for (int iter = 1;iter<50;iter++){
        *un = *u;
        jacobi();
        cout<<"iter = "<<iter<<"\n";
    }
    //write data:
    u->write_matplolib("plot/","u.dat");
}


poisson::~poisson(){
    delete u;
    delete un;
    delete b;
}

 
int main(){

    //grid:
    double Lx = 1.0, Ly = 1.0;
    int nx = 100, ny = 100;
     
    poisson problem(Lx,Ly,nx,ny);
    problem.initialize_b();
     
    problem.solve();


    return 0;
}
