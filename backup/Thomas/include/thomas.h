#ifndef THOMAS_H
#define THOMAS_H

#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

using namespace std;

void thomas(int N,const vector<double>& b,const vector<double>& a,
			    const vector<double>& c,vector<double>& x,const vector<double>& q);


#endif
