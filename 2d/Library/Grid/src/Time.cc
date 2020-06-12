#include "../include/Time.h"

Time::Time(double T,double dt):T_(T),dt_(dt),Nt_(T/dt)
{}

 

const double& Time::T() const{
    return T_;
}

const double& Time::dt() const{
    return dt_;
}

const unsigned int&  Time::Nt() const{
    return Nt_;
}