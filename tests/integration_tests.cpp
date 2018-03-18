#include "integration_tests.hpp"
#include <cmath>
#include <iomanip>  // For setprecision
#include <iostream>

void rastrigin_func(double* x, double* f) /* Rastrigin's  */
{
    /**
     * This function has a global minimum of 0 at x = 0, this is the
     * one and only global minimum
     */
    int nx = 2;  // nx is the dimension of the problem
    f[0] = 0.0;
    for (int i = 0; i < nx; i++) {
        f[0] += (x[i] * x[i] - 10.0 * cos(2.0 * PI * x[i]) + 10.0);
    }
    std::cout << "The fitness is " << std::setprecision(3) << f[0] << "\n";
}
