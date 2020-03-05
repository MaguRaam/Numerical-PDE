#include "../include/helmholtz.h"

helmholtz::helmholtz(double L,int N):L(L),N(N),h(L/double(N)){
    x = new double[N+1];
    y = new double[N+1];

    eigenvalue = new double[N+1];
    CreateMatrix(u);CreateMatrix(uexact);CreateMatrix(P);CreateMatrix(Pinverse);
    CreateMatrix(b);CreateMatrix(btilde);CreateMatrix(temp);


}

void helmholtz::CreateMatrix(double**& A){
    A = new double*[N+1];
    for (int i = 0;i<N+1;i++){
        A[i] = new double[N];
        for (int j = 0;j<N+1;j++)
            A[i][j] = 0.0;
    }
}

void helmholtz::DeleteMatrix(double**& A){
    for (int i = 0;i<N+1;i++){
        delete [] A[i];
    }
    delete [] A;
}

helmholtz::~helmholtz(){
    delete [] x;
    delete [] y;
    delete [] eigenvalue;

    DeleteMatrix(u);DeleteMatrix(uexact);DeleteMatrix(P);DeleteMatrix(Pinverse);
    DeleteMatrix(b);DeleteMatrix(btilde);DeleteMatrix(temp);

}

//create grid:
void helmholtz::SetGrid(){
    for (int i = 0;i<N+1;i++){
        x[i] = i*h; y[i] = i*h;
    }
}

//set exact solution:
void helmholtz::SetExactSolution(){
    for(int i = 0;i<=N;i++){
        for (int j=0;j<=N;j++)
            uexact[i][j] = sin(10.0*M_PI*x[i])*sin(10.0*M_PI*y[j]);
    }
    WriteDataTecplot2D(uexact,"../plot/","uexact.dat");
}

//set right hand side:
void helmholtz::SetRightHandSide(){
    for (int i = 1;i<N;i++){
        for (int j = 1;j<N;j++){
            b[i][j]  = -200.0*M_PI*M_PI*sin(10.0*M_PI*x[i])*sin(10.0*M_PI*y[j]);
            b[i][j] += 1000.0*( sin(10.0*M_PI*x[i])*sin(10.0*M_PI*y[j]) ) ;
            b[i][j] *= h*h; 
        }
    }
    WriteDataTecplot2D(b,"../plot/","b.dat");
}

//set eigenvalues:
void helmholtz::SetEigenValues(){
     
    for (int i=1;i<N;i++)
        eigenvalue[i] = -4.0*sin(M_PI*0.5*x[i])*sin(M_PI*0.5*x[i]);
}

//set P and Pinv matrix:
void helmholtz::SetPandPinv(){
    for (int i=1;i<N;i++){
        for (int j=1;j<N;j++){
            P[i][j] = sin(i*M_PI*x[j]);
            Pinverse[i][j] = 2.0*h*sin(j*M_PI*x[i]) ;
        }
    }
}

//Compute Ftilde = P^{-1} F Q^{-T}
void helmholtz::ComputeFtilde(){
    for (int i = 1;i < N;i++){
        for (int j = 1;j < N;j++){
            temp[i][j] = 0.0;
            for (int k = 1;k < N;k++)
                temp[i][j] += Pinverse[i][k]*b[k][j];
        }
    }

    for (int i = 1;i < N;i++){
        for (int j = 1;j < N;j++){
            btilde[i][j] = 0.0;
            for (int k = 1;k < N;k++)
                btilde[i][j] += temp[i][k]*P[k][j]; 
        }
    }

    WriteDataTecplot2D(btilde,"../plot/","btilde.dat");
}

//Compute Utilde = Ftilde/lambda:
void helmholtz::ComputeUtilde(){
    for (int i = 1;i < N;i++){
        for (int j = 1;j < N;j++){
            btilde[i][j] /= (eigenvalue[i] + eigenvalue[j] + 1000.0*h*h);
        }   
    }
    WriteDataTecplot2D(btilde,"../plot/","utilde.dat");
}

//Compute U = P \bar{U} Q^T
void helmholtz::ComputeU(){
    for (int i = 1;i < N;i++){
        for (int j = 1;j < N;j++){
            temp[i][j] = 0.0;
            for (int k = 1;k < N;k++)
                temp[i][j] += P[i][k]*btilde[k][j];
        }
    }

    for (int i = 1;i < N;i++){
        for (int j = 1;j < N;j++){
            u[i][j] = 0.0;
            for (int k = 1;k < N;k++)
                u[i][j] += temp[i][k]*Pinverse[k][j]; 
        }
    }
    WriteDataTecplot2D(u,"../plot/","u.dat");
}


//write data:
void helmholtz::WriteDataTecplot2D(double**& A,string folderpath,string filename){
    ofstream File(folderpath + filename, ios::out);
    File.flags( ios::dec | ios::scientific );
    File.precision(16) ;
    if(!File)
        cerr<<"Error opening file";
    File<<"TITLE = "+ filename<<"\n"<<"VARIABLES = X, Y," + filename<<"\n";
    File<<"Zone T = psi I = "<<N+1<<"J = "<<N+1<<"\n";
    for (int i = 0;i<=N;i++){
        for (int j = 0;j<=N;j++)
            File<<x[i]<<"\t"<<y[j]<<"\t"<<A[i][j]<<"\t"<<"\n";
    }
    File.close();

}

void helmholtz::WriteDataTecplot1D(double*& b,string folderpath,string filename){
    ofstream File(folderpath + filename, ios::out);
    File.flags( ios::dec | ios::scientific );
    File.precision(16) ;
    if(!File)
        cerr<<"Error opening file";
    for (int i = 0;i<=N;i++)
        File<<b[i]<<"\n";
}

void helmholtz::L2Error(){
    l2error = 0.0;
    for (int i = 0;i<=N;i++){
        for (int j = 0;j<=N;j++)
             l2error += pow(u[i][j] - uexact[i][j],2.0)*h*h;
    }
    l2error = sqrt(l2error);

    cout.flags( ios::dec | ios::scientific );
	cout.precision(5);
    cout<<"Ngpts = "<<N<<"\t\t\tL2Error = "<<l2error<<"\t\t\tTime = "<<time<<"\tsec\n";
}



void helmholtz::Solve(){
    auto start = chrono::steady_clock::now();
    SetGrid();
    SetExactSolution();
    SetRightHandSide();
    SetEigenValues();
    SetPandPinv();
    ComputeFtilde();
    ComputeUtilde();
    ComputeU();
    auto end = chrono::steady_clock::now();
    time = chrono::duration_cast<chrono::seconds>(end - start).count(); 
    L2Error();
}