#include "../include/poisson.h"
#include <vector>


int main(){

	vector<int> n{1,2,5,10,20};
	cout<<"wave number"<<"\t"<<"no of iteration"<<"\n";

	for (unsigned int i=0;i<n.size();i++){

		
		//grid:
		double L = 1.0;
		int N = 100;

		poisson problem(L,N);
		problem.initialize_b(n[i]);
		cout<<n[i]<<"\t\t\t"<<problem.solve()<<"\n";
		problem.exact_solution(n[i]);

	}
	return 0;
}
