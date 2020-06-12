#include <iostream>
#include <cassert>
#include <vector>
#include <cmath>
#include <algorithm>

template <typename Vector>
inline typename Vector::value_type Linfty(const Vector &v)
{
    auto w(v);
    for (auto &elem : w)
        elem = abs(elem);
    return *max_element(w.begin(), w.end());
}




int main()
{
    std::vector<float> v{1.0,-2.0,3.0,-4.0};
    std::cout<<Linfty(v)<<std::endl;
}

/*
//space:
    const size_t N = 100;
    const double L = 1.0, dx = L / N;

    //time:
    const double cfl = 0.5, dt = 0.5 * dx;
    const double T = 5;

    //initialize:
    std::vector<double> uold(N), unew(N);
    double lambda = 0.5;

    for (size_t i = 0; i < N; i++)
        uold[i] = cos((2.0 * M_PI / lambda) * i * dx);

    double t = 0;
    while (t < T - dt)
    {
        t += dt;
        for (size_t i = 0; i < N; i++)
            unew[i] = cfl * (uold[i - 1] - uold[i]) + uold[i];
        uold = unew;
    }
    std::cout << t << std::endl;
*/