#include "../include/field.h"

field::field(double L, int N) : L(L), N(N), h(L / N)
{

    //allocate memory:
    u = new double *[N];
    for (int i = 0; i <= N; i++)
        u[i] = new double[N];

    //initialize to zero:
    for (int j = 0; j <= N; j++)
    {
        for (int i = 0; i <= N; i++)
            u[j][i] = 0.0;
    }
}

field::~field()
{

    //deallocate memory:
    for (int i = 0; i <= N; i++)
        delete[] u[i];
    delete[] u;
}

void field::set_function(double (*func)(const double &, const double &))
{
    for (int j = 0; j <= N; j++)
    {
        for (int i = 0; i <= N; i++)
            u[j][i] = func(i * h, j * h);
    }
}

double field::length()
{
    return L;
}

int field::npts()
{
    return N;
}

double field::grid_size()
{
    return h;
}

double &field::operator()(int j, int i)
{
    return u[j][i];
}

void field::write_output(string folderpath, string filename)
{
    ofstream x, y, u_;
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

    for (int i = 0; i <= N; i++)
        x << i * h << "\n";
    for (int j = 0; j <= N; j++)
        y << j * h << "\n";

    for (int j = 0; j <= N; j++)
    {
        for (int i = 0; i <= N; i++)
        {
            u_ << u[j][i] << " ";
        }

        u_ << "\n";
    }

    //close fiels:
    x.close();
    y.close();
    u_.close();
}

double field::L2norm() const
{
    double value = 0.0;
    for (int j = 0; j <= N; j++)
    {
        for (int i = 0; i <= N; i++)
            value += pow(u[j][i], 2.0);
    }
    return sqrt(value);
}

double field::read(int j, int i) const
{
    return u[j][i];
}

void field::zero()
{
    for (int j = 0; j <= N; j++)
    {
        for (int i = 0; i <= N; i++)
            u[j][i] = 0.0;
    }
}
