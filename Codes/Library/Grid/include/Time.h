#ifndef TIME_H
#define TIME_H

class Time{
    private:
        double T_,dt_;
        unsigned int Nt_;
    public:
        Time(double T,double dt);
        const double& T() const;
        const double& dt() const;
        const unsigned int&    Nt() const;
};



#endif