#include "../include/Plot.h"

 void WriteDataTecplot(const Matrix<double>& A,const Grid& g,const std::string folderpath,const std::string filename){
    std::ofstream File(folderpath + filename, std::ios::out);
    File.flags( std::ios::dec | std::ios::scientific );
    File.precision(16) ;
    if(!File)
        std::cerr<<"Error opening file";
    File<<"TITLE = "+ filename<<"\n"<<"VARIABLES = X, Y," + filename<<"\n";
    File<<"Zone T = psi I = "<<g.N()+1<<"J = "<<g.N()+1<<"\n";
    for (unsigned int i = 0;i<=g.N();i++){
        for (unsigned int j = 0;j<=g.N();j++)
            File<<g.x(i)<<"\t"<<g.y(j)<<"\t"<<A(i,j)<<"\t"<<"\n";
    }
    File.close();
}

double L2(const Matrix<double>& uexact,const Matrix<double>& u,const Grid& g){
    double l2error = 0.0;
    for (unsigned int i = 0;i<=g.N();i++){
        for (unsigned int j = 0;j<=g.N();j++)
             l2error += pow(u(i,j) - uexact(i,j),2.0)*g.h()*g.h();
    }
    return sqrt(l2error);
} 


 