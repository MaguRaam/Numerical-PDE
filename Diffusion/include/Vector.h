#ifndef VECTOR_H
#define VECTOR_H


#include <cassert>

class Vector{
    public:
        Vector(int Size);
        Vector(const Vector& OtherVec);
        ~Vector();
        int Size() const;
        double& operator()(int i);
    private:
        int size;
        double* Vec;

};

#endif