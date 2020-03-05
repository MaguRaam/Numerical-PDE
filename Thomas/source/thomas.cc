#include "../include/thomas.h"


void thomas(int N,const vector<double>& b,const vector<double>& a,
			    const vector<double>& c,vector<double>& x,const vector<double>& q){

	int i;
	vector<double> l(N,0.0);
	vector<double> u(N,0.0);
	vector<double> d(N,0.0);
	vector<double> y(N,0.0);
	
	// LU Decomposition:
	d[0] = a[0];
	u[0] = c[0];
	
	for(i=0;i<N-2;i++){
		l[i] = b[i]/d[i];
		d[i+1] = a[i+1] - l[i]*u[i];
		u[i+1] = c[i+1];
	}
	
	l[N-2] = b[N-2]/d[N-2];
	d[N-1] = a[N-1] - l[N-2]*u[N-2];
	
	//Forward substitution [L][y] = [q] :
	y[0] = q[0];
	for(i=1;i<N;i++)
		y[i] = q[i] - l[i-1]*y[i-1];
	
	/* Backward Substitution [U][x] = [y] */
	x[N-1] = y[N-1]/d[N-1];
	for(i=N-2;i>=0;i--)
		x[i] = (y[i] - u[i]*x[i+1])/d[i];
}