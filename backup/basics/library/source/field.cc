#include "field.h"

using namespace std;

//override default constructor : dynamicaly allocate memory and initialize to zeros:
field::field(double Lengthx, double Lengthy, int ngridx, int ngridy) : Lx(Lengthx), Ly(Lengthy), nx(ngridx), ny(ngridy)
{

    //compute grid size:
    hx = Lx / (nx - 1);
    hy = Ly / (ny - 1);

    u = new double *[nx];
    for (int i = 0; i < nx; i++)
        u[i] = new double[ny];

    //initialize to zeros:
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            u[i][j] = 0.0;
        }
    }
}

//copy constructor:
//Note : use copy constructor for multiple fields:
field::field(const field &other_field)
{

    Lx = other_field.Lx;
    Ly = other_field.Ly;
    nx = other_field.nx;
    ny = other_field.ny;
    nx = other_field.nx;
    ny = other_field.ny;

    u = new double *[nx];
    for (int i = 0; i < nx; i++)
        u[i] = new double[ny];

    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            u[i][j] = other_field.u[i][j];
        }
    }
}

//set field values to be constant:
void field::set_value(double value)
{
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            u[i][j] = value;
        }
    }
}

//set field to given analytical function:
void field::set_function(double (*func)(double, double, double), double t)
{
    double x = 0.0, y = 0.0;

    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {

            x = (i * hx);
            y = (j * hy);

            //field value of the cell:
            u[i][j] = func(x, y, t);
        }
    }
}

//access and modify elements of the field:
double &field::operator()(int i, int j)
{
    assert(i > -1 && i < nx);
    assert(j > -1 && j < ny);
    return u[i][j];
}

//= operator:
field &field::operator=(const field &other_field)
{
    assert(nx = other_field.nx);
    assert(ny = other_field.ny);
    assert(Lx = other_field.Lx);
    assert(Ly = other_field.Ly);

    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
            u[i][j] = other_field.u[i][j];
    }

    return *this;
}

//binary - operator:
field field::operator-(const field& other_field){
    assert(nx==other_field.nx);
    assert(ny==other_field.ny);
    field v(Lx,Ly,nx,ny);
    for (int i=0;i<nx;i++){
        for (int j=0;j<ny;j++){
            v.u[i][j]=u[i][j]+other_field.u[i][j];
        }
    }
    return v;
}

double field::L2norm(){
    double value = 0.0;
    for (int i=0;i<nx;i++)
        for (int j=0;j<ny;j++){
            value+=(pow(u[i][j],2)*hx*hy);
        }
    return pow(value,0.5);
}

double field::grid_sizex()
{
    return hx;
}
double field::grid_sizey()
{
    return hy;
}
double field::lengthx()
{
    return Lx;
}
double field::lengthy()
{
    return Ly;
}

int field::nptsx()
{
    return nx;
}

int field::nptsy()
{
    return ny;
}

 

//destructor:
field::~field()
{
    for (int i = 0; i < nx; i++)
        delete[] u[i];
    delete[] u;
}