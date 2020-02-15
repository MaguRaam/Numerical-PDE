#include "../include/helmholtz.h"

//Right handside of helmholtz equation:
double rhs(const double& x,const double& y){
	return (1000 - 200*M_PI*M_PI)*sin(10*M_PI*x)*sin(10*M_PI*y);  
}

int main(){

	//grid:
	double L = 1.0;
	int N = 20;

	helmholtz problem(L,N);
	problem.initialize_b(rhs);
	problem.set_bcs();

	//STEP-1:Transform b to b~ 
	problem.compute_btilde();

	//STEP-2:Solve for u~ = b~/lambda
	problem.compute_utilde();

	//STEP-3:Transform x~ to x
	problem.compute_u();


	return 0;
}






//problem.set_bcs();

	//STEP-1:Transform b to b~ 
	//problem.compute_btilde();

	//STEP-2:Solve for u~ = b~/lambda
	//problem.compute_utilde();

	//STEP-3:Transform x~ to x
	//problem.compute_u();
