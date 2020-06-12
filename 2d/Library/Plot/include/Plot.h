#ifndef PLOT_H
#define PLOT_H

#include "Grid.h"
#include "Matrix.h"
#include <string>
#include <fstream>
#include <iostream>
#include <cmath>

void WriteDataTecplot(const Matrix<double>& A,const Grid& g,const std::string folderpath,const std::string filename);
void WriteDataMatplot(const Matrix<double>& A,const Grid& g,const std::string folderpath,const std::string filename);
double L2(const Matrix<double>& uexact,const Matrix<double>& u,const Grid& g); 
double Linfty(const Matrix<double>& uexact,const Matrix<double>& u,const Grid& g);


#endif
