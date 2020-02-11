
#include "field.h"

class poisson{
    public:
        poisson(double Lx,double Ly,int nx,int ny);
        void initialize_b();
        void solve();
        ~poisson();
    private:
        field* u_;
        field* un_;
        field* b_;
        double residue;

        void jacobi();

         
};

poisson::poisson(double Lx,double Ly,int nx,int ny){
    u_  = new field(Lx,Ly,nx,ny);
    un_ = new field(*u_);
    b_  = new field(*u_);

    field u  = *u_;
    field un = *un_;
    field b  = *b_; 


      
}

poisson::~poisson(){
    delete u_;
    delete un_;
    delete b_;
}

//right handside 
double rhs(const double& x,const double& y){
   int n = 2;
   return -2.0*n*n*M_PI*M_PI*sin(n*M_PI*x)*sin(n*M_PI*y);
}

//initialize b:
void poisson::initialize_b(){
    b.set_function(rhs);
    b.write_matplolib("plot/","b.dat");
}

//jacobi:
void poisson::jacobi(){
    for (int j=1;j<n;j++){
        for (int i=1;i<n;i++)
            u(j,i) = 0.25*(un(j-1,i) + un(j+1,i) + un(j,i+1) + un(j,i-1) - b(j,i)*(dx*dx) );
    }
}


//solve:
void solve(){
    for (int iter=1;iter<100;iter++){
        un = u;


        cout<<"iteration = "<<iter<<"\n";
    }
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
