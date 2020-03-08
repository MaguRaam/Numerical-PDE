#ifndef FIELD_H
#define FIELD_H

#include "Matrix.h"
#include "Grid.h"
#include <iostream>
#include <fstream>
#include <string>

//The Field Class is inherited from matrix and grid class:

class Field : public Matrix, public Grid
{

public:
    Field(int Nx, int Ny, double Lx, double Ly);
    Field(const Field &OtherField);
    void WriteDataTecplot(std::string folderpath, std::string filename);
    void SetFunction(double (*func)(double x, double y, double t), double t);
    Field operator*(const Field& OtherField) const;
};

#endif