#include <iostream>
#include <cmath>

#include "Fields/ScalarField3D.hpp"
#include "Operators/Gradient.hpp"

int main()
{
    const std::size_t N = 20;

    double h = 0.1;

    Onko::ScalarField3D f(N, N, N);

    // Fill the field
    for (std::size_t k = 0; k < N; ++k)
    {
        for (std::size_t j = 0; j < N; ++j)
        {
            for (std::size_t i = 0; i < N; ++i)
            {
                double x = i * h;
                double y = j * h;
                double z = k * h;

                f(i, j, k) = std::sin(x)
                            + std::sin(y)
                            + std::sin(z);
            }
        }
    }

    // Interior point
    std::size_t i = 8;
    std::size_t j = 9;
    std::size_t k = 7;

    double x = i * h;
    double y = j * h;
    double z = k * h;

    Onko::Vector3D numerical = Onko::Gradient(f, i, j, k, h);

    // Exact gradient
    double exact_x = std::cos(x);
    double exact_y = std::cos(y);
    double exact_z = std::cos(z);

    std::cout << "Numerical gradient:\n";
    std::cout << "Dx = " << numerical.x << '\n';
    std::cout << "Dy = " << numerical.y << '\n';
    std::cout << "Dz = " << numerical.z << '\n';

    std::cout << "\nExact gradient:\n";
    std::cout << "Dx = " << exact_x << '\n';
    std::cout << "Dy = " << exact_y << '\n';
    std::cout << "Dz = " << exact_z << '\n';

    std::cout << "\nErrors:\n";
    std::cout << "Error x = "
              << std::abs(numerical.x - exact_x) << '\n';

    std::cout << "Error y = "
              << std::abs(numerical.y - exact_y) << '\n';

    std::cout << "Error z = "
              << std::abs(numerical.z - exact_z) << '\n';
}