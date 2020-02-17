#include "../include/helmholtz.h"

double rhs(const double& x,const double& y){
	 return (1000.0-(200.0*M_PI*M_PI))*(sin(10.0*M_PI*x))*(cos(10.0*M_PI*y));
}

double exact_solution(const double& x,const double& y){	 
	return sin(10.0*M_PI*x)*cos(10.0*M_PI*y);
}

void convergence(int N){
	//grid:
	double L = 1.0;

	//solves helmholtz equation:
	helmholtz problem(L,N);
	problem.initialize_b(rhs);
	problem.compute_btilde();
	problem.compute_utilde();
	problem.compute_u();

	//exact solution:
    field uexact(L,N);
    uexact.set_function(exact_solution);
    uexact.write_output("../plot/","uexact.dat");

	cout<<N<<"\t\t\t"<<problem.L2error()<<"\t\t"<<problem.Linfyerror()<<"\n";
}



int main(){
	int N = 10;
	for (int i = 1;i<5;i++){
		N = 2*N;
		convergence(N);
	}
 
}
