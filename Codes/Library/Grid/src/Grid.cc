#include "../include/Grid.h"

 

Grid::Grid(unsigned int N, double L) : N_(N), L_(L), h_(L / double(N))
{
    x_.resize(N_ + 1, 0.0);
    for (unsigned int i = 0; i <= N_; i++)
        x_[i] = i * h_;
}

const double &Grid::x(const unsigned int &i) const
{
    return this->x_[i];
}

const double &Grid::y(const unsigned int &j) const
{
    return this->x_[j];
}

const unsigned int &Grid::N() const
{
    return N_;
}

const double &Grid::h() const
{
    return this->h_;
}

const double &Grid::L() const
{
    return this->L_;
}

Matrix<double> SetFunction(double (*func)(const double &x, const double &y), const Grid &g)
{
    Matrix<double> u(g.N() + 1, g.N() + 1, 0.0);
    for (unsigned int i = 0; i <= g.N(); i++)
    {
        for (unsigned int j = 0; j <= g.N(); j++)
            u(i, j) = func(g.x(i), g.y(j));
    }
    return u;
}
