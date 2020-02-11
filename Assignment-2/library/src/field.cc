#include "field.h"


//constructor:
field::field(double lx,double ly,int nx,int ny):Lx(lx),Ly(ly),Nx(nx),Ny(ny){
    dx = Lx/nx;
    dy = Ly/ny;

    //create memory:
    u = new double*[Ny+1];
    for (int j=0;j<Ny+1;j++)
        u[j] = new double[Nx+1];
    
    //initialize with 0:
    for (int j=0;j<Ny+1;j++){
        for (int i=0;i<Nx+1;i++)
            u[j][i] = 0.0;
    }
}

//copy constructor:
field::field(const field& v):Lx(v.Lx),Ly(v.Ly),Nx(v.Nx),Ny(v.Ny),dx(v.dx),dy(v.dy){
    //create memory:
    u = new double*[Ny+1];
    for (int j=0;j<Ny+1;j++)
        u[j] = new double[Nx+1];

    //copy values:
    for (int j=0;j<Ny+1;j++){
        for (int i=0;i<Nx+1;i++)
            u[j][i] = v.u[j][i];
    }
}

//destructor:
field::~field(){
    for (int j=0;j<Ny+1;j++)
        delete [] u[j];
    delete [] u;
}

//set value:
void field::set_value(const double value){
    for (int j=0;j<Ny+1;j++){
        for (int i=0;i<Nx+1;i++)
            u[j][i] = value;
    }
}

//set function:
void field::set_function(double (*func)(const double&,const double&)){
    for (int i=0;i<=Nx;i++){
         for (int j=0;j<=Ny;j++){
            u[j][i] = func(i*dx,j*dy); 
         }
    }
}

//access elements u[j,i]:
double& field::operator()(const int j,const int i){
    assert(i>=0 && i<=Nx);
    assert(j>=0 && j<=Ny);
    return u[j][i];
}

//asssignment operator u = v:
field& field::operator=(const field& v){
    assert(Lx == v.Lx);
    assert(Ly == v.Ly);
    assert(Nx == v.Nx);
    assert(Ny == v.Ny);

    for (int i=0;i<=Nx;i++){
         for (int j=0;j<=Ny;j++){
            u[j][i] = v.u[j][i]; 
         }
    }


    return *this;
}

//grid size:
double field::hx()const{
	return dx;
}


double field::hy()const{
	return dy;
}

//no of grid points:
int field::nx()const{
	return Nx;
}

int field::ny()const{
	return Ny;
}



//norm:
double L2(field& u){
    double value = 0.0;
    for (int i=0;i<=u.Nx;i++){
         for (int j=0;j<=u.Ny;j++){
            value += pow(u(j,i),2.0)*u.dx*u.dy;
         }
    }
    return sqrt(value);
}


//write data in matplotlib format:
void field::write_matplolib(string folderpath, string filename){
    ofstream x,y,u_;
	
	 string x_filename = folderpath + "x.dat"; 
     string y_filename = folderpath + "y.dat";
     string u_filename = folderpath + filename;

    //open files:
    x.open(x_filename.c_str());
    y.open(y_filename.c_str());
    u_.open(u_filename.c_str());

    //assert file is open:
    assert(x.is_open());
    assert(y.is_open());
    assert(u_.is_open());

    //set control on output data
    x.setf(std::ios::scientific);
    x.setf(std::ios::showpos);
    x.precision(13);

    y.setf(std::ios::scientific);
    y.setf(std::ios::showpos);
    y.precision(13);

    u_.setf(std::ios::scientific);
    u_.setf(std::ios::showpos);
    u_.precision(13);
    
    for (int i=0; i<Nx+1; i++)  {
            x << i*dx <<"\n"; 
        }

    for (int j=0; j<Ny+1; j++)  {
            y<< j*dy <<"\n"; 
    }
    for (int j=0; j<Ny+1; j++) {
        for (int i=0; i<Nx+1; i++) {
            u_ << u[j][i] << " " ;
        }
        
        u_<<"\n"; 
    }

    //close fiels:
    x.close();y.close();u_.close();
}
