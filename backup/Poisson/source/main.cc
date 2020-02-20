#include "../include/helmholtz.h"

double rhs(const double& x,const double& y){
   int N = 20;       //TODO
	double h = 1.0/N;
	int p = 10;
	int q = 10;
	double norm = (pow(2.0/N,0.5))*(pow(2.0/N,0.5));
	double eigenvalue = -( (4.0/(h*h))*(pow(sin((M_PI*q)/(2*N)),2))    ) - ( (4.0/(h*h))*(pow(sin((M_PI*p)/(2*N)),2)) );
	return eigenvalue*norm*sin(p*M_PI*x)*sin(q*M_PI*y);
}

double exact_solution(const double& x,const double& y){
    int N = 20;  //TODO
	int p = 10;
	int q = 10;
	double norm = (pow(2.0/N,0.5))*(pow(2.0/N,0.5));
	 
	return norm*sin(p*M_PI*x)*sin(q*M_PI*y);
}


int main(){
	
	//grid:
	double L = 1.0;
	int N = 20; //TODO

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

	cout<<"L2error = "<<problem.L2error()<<"\n";


	
	 
}
