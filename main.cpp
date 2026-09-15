#include <iostream>
#include <cmath>
#include <iomanip>

#include "Fields/VectorField3D.hpp"
#include "Operators/Divergence.hpp"

int main()
{
    const double pi = std::acos(-1.0);

    // Grid sizes
    const std::size_t Nx = 25;
    const std::size_t Ny = 25;
    const std::size_t Nz = 25;

    // Grid spacing
    const double h = 0.1;

    // Create vector field
    Onko::VectorField3D u(Nx, Ny, Nz);

    // --------------------------------------------------
    // Define analytical vector field:
    //
    // u_x = sin(x)
    // u_y = sin(y)
    // u_z = sin(z)
    // --------------------------------------------------

    for (std::size_t k = 0; k < Nz; ++k)
    {
        for (std::size_t j = 0; j < Ny; ++j)
        {
            for (std::size_t i = 0; i < Nx; ++i)
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

    // --------------------------------------------------
    // Test at one interior point
    // --------------------------------------------------

    std::size_t i = 10;
    std::size_t j = 10;
    std::size_t k = 10;

    double numerical = Onko::Divergence(u, i, j, k, h);

    double x = i * h;
    double y = j * h;
    double z = k * h;

    // Exact divergence:
    //
    // d(sin x)/dx + d(sin y)/dy + d(sin z)/dz
    // = cos x + cos y + cos z

    double exact =
        std::cos(x)
        + std::cos(y)
        + std::cos(z);

    double error = std::abs(numerical - exact);

    std::cout << std::setprecision(15);

    std::cout << "Numerical divergence = "
              << numerical << '\n';

    std::cout << "Exact divergence     = "
              << exact << '\n';

    std::cout << "Absolute error       = "
              << error << '\n';

    return 0;
}