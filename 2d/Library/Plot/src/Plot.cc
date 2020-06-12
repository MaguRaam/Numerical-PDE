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


void WriteDataMatplot(const Matrix<double>& A,const Grid& g,const std::string folderpath,const std::string filename){
    std::ofstream mpl_x;
    std::ofstream mpl_y;
    std::ofstream mpl_u;

    mpl_x.open (folderpath + "/x.dat", std::ios::out);
    mpl_y.open (folderpath + "/y.dat", std::ios::out);
    mpl_u.open (folderpath + filename, std::ios::out);

    if ( !(mpl_x.is_open()) || !(mpl_y.is_open()) || !(mpl_u.is_open()) ) {
        
        std::cerr<<"Error opening file";
    }

    else {
        

        for (unsigned int i=0; i<g.N(); i++)  {
            mpl_x << g.x(i) << std::endl; 
        }

        for (unsigned int j=0; j<g.N(); j++)  {
            mpl_y << g.y(j) << std::endl; 
        }

        for (unsigned int i=0; i<g.N(); i++) {
            for (unsigned int j=0; j<g.N(); j++) {
                mpl_u << A(i,j) << " " ;
            }
            
            mpl_u << std::endl; 
        }
    }

    mpl_u.close();

}



double L2(const Matrix<double>& uexact,const Matrix<double>& u,const Grid& g){
    double l2error = 0.0;
    for (unsigned int i = 0;i<=g.N();i++){
        for (unsigned int j = 0;j<=g.N();j++)
             l2error += pow(u(i,j) - uexact(i,j),2.0)*g.h()*g.h();
    }
    return sqrt(l2error);
} 

double Linfty(const Matrix<double>& uexact,const Matrix<double>& u,const Grid& g){
    double maxerror = 0.0;
    for (unsigned int i = 0;i<=g.N();i++){
        for (unsigned int j = 0;j<=g.N();j++)
             if( fabs(u(i,j) - uexact(i,j) ) > maxerror) maxerror = fabs( u(i,j) -uexact(i,j) );
    }

    return maxerror;
 }
