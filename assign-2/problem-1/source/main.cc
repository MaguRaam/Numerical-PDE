#include "../include/thomas.h" 

double L2error(const vector<double>& u,const vector<double>& uexact){
	double error = 0.0;
	int N = u.size()+1; 
	for (unsigned int i = 0;i<u.size();i++)
		error += pow( (u[i]-uexact[i] ),2.0)*(1.0/N);
	return sqrt(error);
}



void write_data(const vector<double>& u,const string filename,double h,int beta){
    ofstream output(filename);
    output<<0.0<<"\t\t\t"<<0.0<<"\n";
    for (unsigned int i=0;i<u.size();i++)
        output<<(i+1)*h<<"\t\t\t"<<(1.0/(beta*beta*M_PI*M_PI))*u[i]<<"\n";
    output<<1.0<<"\t\t\t"<<0.0<<"\n";
    output.close();
}

double poisson(int N,int beta){
	//grid:
	double L = 1.0;
	double h = L/N;

	//initialize rhs vector q:
	vector<double> q(N-1,0.0);
	for (int i=1;i<N;i++)
		q[i-1] = -sin(beta*M_PI*i*h)*h*h;
	
	//build tridiagonal matrix:
	vector<double> a(N-1,-2.0);
	vector<double> b(N-1,1.0);
	vector<double> c(N-1,1.0);

	//solution vector:
	vector<double> x(N-1,0.0);
	thomas(N-1,b,a,c,x,q);

	//exact solution:
	vector<double> xexact(N-1,0.0);
	for (int i=1;i<N;i++)
		xexact[i-1] = (1.0/(beta*beta*M_PI*M_PI))*sin(beta*M_PI*i*h);
	
	//write solution:
	write_data(x,"../plot/u100.dat",h,beta);
	write_data(xexact,"../plot/uexact100.dat",h,beta);

	return L2error(x,xexact);
}



int main(){

	//ngpts:
	vector<int> N{25,50,100,200,400};

	//beta:
	vector<int> beta{1,10,100};


	//Part a:
	cout<<"ngpts"<<"\t\t\t"<<"L2error"<<"\n";
	 
    cout<<N[1]<<"\t\t\t"<<poisson(N[1],beta[2])<<"\n";
 



	 
	 


	 


 return 0;
}
