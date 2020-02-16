#include "../include/functions.h"

double L2error(const vector<double>& u,const vector<double>& uexact){
	double error = 0.0;
	int N = u.size()+1; 
	for (unsigned int i = 0;i<u.size();i++)
		error += pow( (u[i]-uexact[i] ),2.0)*(1.0/N);
	return sqrt(error);
}



void write_data(const vector<double>& u,const string filename){
    ofstream output(filename);
    for (unsigned int i=0;i<u.size();i++)
        output<<u[i]<<"\n";
    output.close();
}


