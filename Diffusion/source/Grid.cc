#include "../include/Grid.h"

Grid::Grid(int Nx, int Ny, double Lx, double Ly) : Nx(Nx),
                                                   Ny(Ny),
                                                   Lx(Lx),
                                                   Ly(Ly),
                                                   hx(double(Lx / Nx)),
                                                   hy(double(Ly / Ny))
{
    x = new Vector(Nx + 1);
    y = new Vector(Ny + 1);

    for (int i = 0; i <= Nx; i++)
        (*x)(i) = i * hx;
    for (int j = 0; j <= Ny; j++)
        (*y)(j) = j * hy;
}

Grid::Grid(const Grid &OtherGrid) : Nx(OtherGrid.Nx),
                                    Ny(OtherGrid.Ny),
                                    Lx(OtherGrid.Lx),
                                    Ly(OtherGrid.Ly),
                                    hx(OtherGrid.hx),
                                    hy(OtherGrid.hy)
{
    x = new Vector(Nx + 1);
    y = new Vector(Ny + 1);

    for (int i = 0; i <= Nx; i++)
        (*x)(i) = OtherGrid.X(i);
    for (int j = 0; j <= Ny; j++)
        (*y)(j) = OtherGrid.Y(j);
}

double Grid::X(int i) const
{
    assert(i > -1);
    assert(i <= Nx);
    return (*x)(i);
}

double Grid::Y(int j) const
{
    assert(j > -1);
    assert(j <= Ny);
    return (*y)(j);
}

int Grid::nx() const
{
    return Nx;
}

int Grid::ny() const
{
    return Ny;
}

double Grid::dx() const
{
    return hx;
}
double Grid::dy() const
{
    return hy;
}

double Grid::lx() const
{
    return Lx;
}
double Grid::ly() const
{
    return Ly;
}

Grid::~Grid()
{
    delete x;
    delete y;
}