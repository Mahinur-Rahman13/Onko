#include <iostream>
#include <cmath>
#include <iomanip>
#include <algorithm>

#include "Fields/VectorField3D.hpp"
#include "Operators/Divergence.hpp"

double testDivergence(double h)
{
    const std::size_t N = 25;

    Onko::VectorField3D u(N, N, N);

    // Fill vector field
    for (std::size_t k = 0; k < N; ++k)
    {
        for (std::size_t j = 0; j < N; ++j)
        {
            for (std::size_t i = 0; i < N; ++i)
            {
                double x = i * h;
                double y = j * h;
                double z = k * h;

                u.x()(i, j, k) = std::sin(x);
                u.y()(i, j, k) = std::sin(y);
                u.z()(i, j, k) = std::sin(z);
            }
        }
    }

    double maxError = 0.0;

    // Need two points on each side because
    // your derivative stencil is fourth-order central.
    for (std::size_t k = 2; k < N - 2; ++k)
    {
        for (std::size_t j = 2; j < N - 2; ++j)
        {
            for (std::size_t i = 2; i < N - 2; ++i)
            {
                double x = i * h;
                double y = j * h;
                double z = k * h;

                double numerical =
                    Onko::Divergence(u, i, j, k, h);

                double exact =
                    std::cos(x)
                    + std::cos(y)
                    + std::cos(z);

                double error =
                    std::abs(numerical - exact);

                maxError = std::max(maxError, error);
            }
        }
    }

    return maxError;
}

int main()
{
    std::cout << std::setprecision(10);

    double h1 = 0.2;
    double h2 = 0.1;
    double h3 = 0.05;
    double h4 = 0.025;

    double e1 = testDivergence(h1);
    double e2 = testDivergence(h2);
    double e3 = testDivergence(h3);
    double e4 = testDivergence(h4);

    std::cout << "h = " << h1 << "  Error = " << e1 << '\n';
    std::cout << "h = " << h2 << "  Error = " << e2 << '\n';
    std::cout << "h = " << h3 << "  Error = " << e3 << '\n';
    std::cout << "h = " << h4 << "  Error = " << e4 << '\n';

    std::cout << "\nError ratios:\n";

    std::cout << e1 / e2 << '\n';
    std::cout << e2 / e3 << '\n';
    std::cout << e3 / e4 << '\n';

    return 0;
}