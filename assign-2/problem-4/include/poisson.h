#ifndef POISSON_H
#define POISSON_H

#include "field.h"
class poisson{
    public:
        poisson(double l, int n);
        ~poisson();
        void initialize_b(int n);
        int solve();
        void exact_solution(int n);
         

    private:
        double L;
        int N;
        double h;
        field* b;
        field* u;
        field* un;
        

        

};



#endif