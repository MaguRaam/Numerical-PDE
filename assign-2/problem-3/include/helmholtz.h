#ifndef HELMHOLTZ_H
#define HELMHOLTZ_H

#include "field.h"

class helmholtz
{
public:
    helmholtz(double l, int n);
    ~helmholtz();
    void initialize_b(double (*func)(const double &, const double &));
    void compute_btilde();
    void compute_utilde();
    void compute_u();
    double L2error();
     


private:
    double L;
    int N;
    double h;
    field *b;
    field *btilde;
    field *&utilde = btilde;
    field *&u = b;

    double project_field_on_eigenbasis(int q, int p);
    double eigenvalue(int q, int p);
    double expand_field_using_eigenbasis(int j, int i);
    void set_bcs();
};

#endif
