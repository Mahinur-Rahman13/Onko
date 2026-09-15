#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstddef>

#include "Fields/ScalarField3D.hpp"
#include "Operators/Laplacian.hpp"

int main()
{
    std::cout << std::fixed << std::setprecision(10);

    std::cout << "============================================\n";
    std::cout << " FOURTH-ORDER LAPLACIAN TEST\n";
    std::cout << "============================================\n\n";

    // Point at which we evaluate the Laplacian.
    const double x = 0.7;
    const double y = 0.8;
    const double z = 0.9;

    // Grid spacings for the convergence test.
    const double spacings[] = {0.2, 0.1, 0.05, 0.025};

    // Our fourth-order central stencil requires
    // two points on each side of the point being tested.
    const std::size_t Nx = 5;
    const std::size_t Ny = 5;
    const std::size_t Nz = 5;

    // Center point of the grid.
    const std::size_t i = 2;
    const std::size_t j = 2;
    const std::size_t k = 2;

    double previous_error = 0.0;

    for (double h : spacings)
    {
        Onko::ScalarField3D f(Nx, Ny, Nz);

        // Define
        //
        // f(x,y,z) = sin(x) + sin(y) + sin(z)
        //
        // throughout the grid.

        for (std::size_t ii = 0; ii < Nx; ++ii)
        {
            for (std::size_t jj = 0; jj < Ny; ++jj)
            {
                for (std::size_t kk = 0; kk < Nz; ++kk)
                {
                    double xi =
                        x + (static_cast<double>(ii) - 2.0) * h;

                    double yj =
                        y + (static_cast<double>(jj) - 2.0) * h;

                    double zk =
                        z + (static_cast<double>(kk) - 2.0) * h;

                    f(ii, jj, kk) =
                        std::sin(xi)
                        + std::sin(yj)
                        + std::sin(zk);
                }
            }
        }

        // Numerical Laplacian.
        double numerical =
            Onko::Laplacian(f, i, j, k, h);

        // Exact Laplacian:
        //
        // ∇²f = -sin(x) - sin(y) - sin(z)

        double exact =
            -std::sin(x)
            -std::sin(y)
            -std::sin(z);

        // Absolute error.
        double error =
            std::abs(numerical - exact);

        std::cout << "h = " << h
                  << "    numerical = " << numerical
                  << "    exact = " << exact
                  << "    error = " << error;

        // For a fourth-order method, when h is halved:
        //
        // error(h) / error(h/2) → 2^4 = 16

        if (previous_error != 0.0)
        {
            double ratio =
                previous_error / error;

            std::cout << "    ratio = " << ratio;
        }

        std::cout << '\n';

        previous_error = error;
    }

    std::cout << "\nTest completed.\n";

    return 0;
}