#include<iostream>
#include<cstdlib>
#include<fstream>
#include<cstring>
#include<math.h>

using namespace std ;

double const PI = 4.0*atan(1.0); // Value of PI

void Allocate_2D_R(double**& m, int d1, int d2) {
        m=new double* [d1];
        for (int i=0; i<d1; ++i) {
                m[i]=new double [d2];
                for (int j=0; j<d2; ++j)
                        m[i][j]=0.0;
        }
}

int main() {

	cout.flags( ios::dec | ios::scientific );
	cout.precision(5);

	double *x, *uExact, **P, **PINV, *PEig, *u, *RHS ;	 
	double **A, **B, *RHS_Tilde, *u_Tilde, *Temp ;
	double L2_Error = 0.0, Max_Error = 0.0, hx ;

	int Nx = 200, i, j, k ;

	x = new double[Nx+1] ;
	uExact = new double[Nx+1] ; u = new double[Nx+1] ; 
	Allocate_2D_R(P, Nx-1, Nx-1) ; Allocate_2D_R(PINV, Nx-1, Nx-1) ; PEig = new double[Nx-1] ;
	Allocate_2D_R(A, Nx-1, Nx-1) ; Allocate_2D_R(B, Nx-1, Nx-1) ;
	RHS = new double[Nx-1] ; Temp = new double[Nx-1] ; RHS_Tilde = new double[Nx-1] ; u_Tilde = new double[Nx-1] ;

	// set a uniform grid.
	for(i = 0 ; i <= Nx ; i++) x[i] = i/double(Nx) ;	
	hx = 1.0/double(Nx) ; 

	// set initial solution and Exact solution for computing error.
	for(i = 0 ; i <= Nx ; i++) {
		u[i] = 0.0 ; uExact[i] = sin(10.0*PI*x[i]) ;
	}

	// Set up eigenvalues and matrices
	for(i = 1 ; i < Nx ; i++) {
		PEig[i-1] = -4.0*sin(PI*0.5*x[i])*sin(PI*0.5*x[i]) ;
		for(j = 1 ; j < Nx ; j++) {
			P[i-1][j-1] = sin(i*PI*x[j]) ; 
			PINV[i-1][j-1] = 2.0*hx*sin(j*PI*x[i]) ;
			B[i-1][j-1] = 0.0 ;
		}
	}
	for(i = 1 ; i < Nx ; i++) {
		if(i == 1) {B[i-1][i-1] = -2.0 ; B[i-1][i] = 1.0; }
		else if(i == Nx-1) {B[i-1][i-1] = -2.0 ; B[i-1][i-2] = 1.0 ; }
		else {B[i-1][i-1] = -2.0 ; B[i-1][i-2] = 1.0 ; B[i-1][i] = 1.0 ;}
	}

	for(i = 1 ; i < Nx ; i++) {
		for(j = 1 ; j < Nx ; j++) {
			A[i-1][j-1] = 0.0 ;
			for(k = 1 ; k < Nx ; k++) A[i-1][j-1] += PINV[i-1][k-1]*P[k-1][j-1] ;
		//	if ( fabs(A[i-1][j-1] - P[i-1][j-1]*PEig[j-1]) > 1.0E-10) cout << i << "\t" << j << "\t" << fabs(A[i-1][j-1] - P[i-1][j-1]*PEig[j-1]) << endl ;
		//	cout << A[i-1][j-1] << "\t" ;
		}
	}
	cout << endl ;
//	for(i = 1 ; i < Nx ; i++) cout << PEig[i-1] << "\t" ;
	cout << endl ;

	// Set the right hand side
	for(i = 1 ; i < Nx ; i++) {
		RHS[i-1] = -100.0*PI*PI*sin(10.0*PI*x[i]) ;
		RHS[i-1] *= (hx*hx) ;
	}

	// compute P^{-1} F Q^{-T}
	for(i = 1 ; i < Nx ; i++) {
		RHS_Tilde[i-1] = 0.0 ;
		for(j = 1 ; j < Nx ; j++) {
			RHS_Tilde[i-1] += PINV[i-1][j-1]*RHS[j-1] ;
		}
	}
	
	// compute U_Tilde
	for(i = 1 ; i < Nx ; i++) {
		RHS_Tilde[i-1] = RHS_Tilde[i-1]/(PEig[i-1]) ;
	}
	for(i = 1 ; i < Nx ; i++) {
		u[i] = 0.0 ;
		for(j = 1 ; j < Nx ; j++) {
			u[i] += P[i-1][j-1]*RHS_Tilde[j-1] ;
		}
	}

	ofstream File("Output.dat", ios::out) ;
	File.flags( ios::dec | ios::scientific );
	File.precision(16) ;
	if(!File) {cerr<< "Error: Output file couldnot be opened.\n";}
	
	Max_Error = L2_Error = 0.0 ;

	for(i = 0 ; i <= Nx ; i++) {
		if( fabs(u[i] - uExact[i] ) > Max_Error) Max_Error = fabs( u[i] - uExact[i] );
		L2_Error += ( u[i] - uExact[i])*( u[i] - uExact[i] )/ ( ( Nx+1.0 ) ) ;

		File << x[i] << "\t" << u[i] << "\t" << uExact[i] << endl ;
	}
	L2_Error = sqrt(L2_Error) ;
	File.close() ;
	cout << "\n L2 : " << L2_Error << "\t Max : " << Max_Error <<  endl ;

	// Free all the memory and exit
	delete[] x ; 

	return 0;
}
