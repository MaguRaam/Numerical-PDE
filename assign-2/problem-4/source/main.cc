#include "../include/poisson.h"
#include <vector>


int main(){

	vector<int> n{1,2,5,10,20};

	for (unsigned int i=0;i<n.size();i++){

		cout<<n[i]<<"\n";
		//grid:
		double L = 1.0;
		int N = 100;
		poisson problem(L,N);

	}
	return 0;
}
