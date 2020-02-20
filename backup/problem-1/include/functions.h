#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>
#include <cmath>
#include <fstream>
#include <string>

double L2error(const vector<double>& u,const vector<double>& uexact);
void write_data(const vector<double>& u,const string filename);


#endif 
