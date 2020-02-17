#include "../include/helmholtz.h"

helmholtz::helmholtz(double l,int n):L(l),N(n),h(L/N){
    
    //allocate memory for fields:
    b       =   new field(L,N);
    btilde  =   new field(L,N);
    utilde  =   new field(L,N);
    u       =   new field(L,N);

}



helmholtz::~helmholtz(){
	//delete memory
    delete b;
    delete btilde;
    delete utilde;
    delete u;
}


//initialize rhs of helmholtz equation:
void helmholtz::set_bcs_b(){
    for (int i=1;i<N;i++){
        (*b)(1,i) -= sin(10*M_PI*i*h)*(1.0/(h*h));
        (*b)(N-1,i) -= sin(10*M_PI*i*h)*(1.0/(h*h));
    }
}


void helmholtz::initialize_b(double (*func)(const double &, const double &)){
    b->set_function(func);
    b->write_output("../plot/","b.dat");
    set_bcs_b();
    b->write_output("../plot/","bmodified.dat");
     
}

//eigenvalue: //TODO
double helmholtz::eigenvalue(int q,int p){
    return -( (4.0/(h*h))*(pow(sin((M_PI*q)/(2*N)),2.0))    ) - ( (4.0/(h*h))*(pow(sin((M_PI*p)/(2*N)),2.0)) )+ 1000.0;
}


//STEP-1:Compute b~
double helmholtz::project_field_on_eigenbasis(int q,int p){
   double value = 0.0;
   double norm = (pow(2.0/N,0.5))*(pow(2.0/N,0.5));

   for (int j=1;j<N;j++){
        for (int i=1;i<N;i++)
            value += (norm*sin(p*M_PI*i*h)*sin(q*M_PI*j*h)*(b->read(j,i)) );
    }  
 
   return value; 
}

void helmholtz::compute_btilde(){
    for (int q=1;q<N;q++){
        for (int p=1;p<N;p++)
            (*btilde)(q,p) = project_field_on_eigenbasis(q,p);
    }
    btilde->write_output("../plot/","btilde.dat");
}

//STEP-2: Solve u~ = b~/eigenvalue;
void helmholtz::compute_utilde(){

    for (int q=1;q<N;q++){
        for (int p=1;p<N;p++){
            (*utilde)(q,p) = (btilde->read(q,p))/eigenvalue(q,p);
        }
    }
    utilde->write_output("../plot/","utilde.dat");
}

//STEP-3: Compute u:
void helmholtz::set_bcs_u(){
    for (int i=1;i<N;i++){
        (*u)(0,i)    = sin(10*M_PI*i*h);
        (*u)(N,i)    = sin(10*M_PI*i*h);
    }
}


double helmholtz::expand_field_using_eigenbasis(int j,int i){
    double value = 0.0;
    double norm = (pow(2.0/N,0.5))*(pow(2.0/N,0.5));
    for (int q=1;q<N;q++){
                for (int p=1;p<N;p++){
                    value+= norm*sin(p*M_PI*i*h)*sin(q*M_PI*j*h)*(utilde->read(q,p));
                }
                  
            }
    return value;
}

void helmholtz::compute_u(){
    for (int j=1;j<N;j++){
        for (int i=1;i<N;i++){
            (*u)(j,i) = expand_field_using_eigenbasis(j,i);
        }
              
    }
    set_bcs_u();
    u->write_output("../plot/","u.dat");
     
}

double helmholtz::L2error(){
    double error = 0.0;
    for (int j=0;j<=N;j++){
        for (int i=0;i<=N;i++){
            error += pow( ( (u->read(j,i)) - (sin(10.0*M_PI*i*h)*cos(10.0*M_PI*j*h))),2.0)*h*h;
            //cout<<"error = "<<error<<"\n";
        }
    }

    return sqrt(error);
}

double helmholtz::Linfyerror(){
    double error = 0.0;
    double value = 0.0;
    for (int j=0;j<=N;j++){
        for (int i=0;i<=N;i++){
             value = abs(u->read(j,i) - (sin(10.0*M_PI*i*h)*cos(10.0*M_PI*j*h)) );
             if (value > error){
                 error = value;
             } 
        }
    }
    return error;
}





