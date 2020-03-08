#include "../include/Vector.h"

//Constructor:
Vector::Vector(int size):size(size){
    assert(size>0);
    Vec = new double[size];
    for (int i = 0;i<size;i++)
        Vec[i] = 0.0;
}

//Copy Constructor:
Vector::Vector(const Vector& OtherVec):size(OtherVec.Size()){
    Vec = new double[size];
    for (int i = 0;i<size;i++)
        Vec[i] = OtherVec.Vec[i];
}

//Destructor:
Vector::~Vector(){
    delete [] Vec;
}

//size:
int Vector::Size() const{
    return size;
}

//operator 0 based indexing ():
double& Vector::operator()(int i){
    assert(i > -1); assert(i < size);
    return Vec[i];
}

 