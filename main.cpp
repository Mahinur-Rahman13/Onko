#include <iostream>
#include <cmath>
#include <iomanip>

#include "Fields/ScalarField3D.hpp"
#include "Operators/FirstDerivative.hpp"

int main()
{
    // Grid size
    const std::size_t N = 100;

    // Test function:
    //
    // f(x) = sin(x)
    //
    // Exact derivative:
    //
    // f'(x) = cos(x)

    const double x0 = 1.0;

    double h_values[] = {
        0.2,
        0.1,
        0.05,
        0.025
    };

    std::cout << std::fixed << std::setprecision(10);

    std::cout << "Fourth-order first derivative test\n\n";

    double previous_error = 0.0;

    for (double h : h_values)
    {
        // Create a 1D-like 3D scalar field.
        // We only vary the x-direction.
        Onko::ScalarField3D f(N, 1, 1);

        // Fill the field with f(x) = sin(x)
        for (std::size_t i = 0; i < N; ++i)
        {
            double x = i * h;

            f(i, 0, 0) = std::sin(x);
        }

        // Find the grid point closest to x0
        std::size_t i = static_cast<std::size_t>(x0 / h);

        // Numerical derivative
        double numerical = Onko::Dx(f, i, 0, 0, h);

        // Actual coordinate of the grid point
        double x = i * h;

        // Exact derivative
        double exact = std::cos(x);

        // Absolute error
        double error = std::abs(numerical - exact);

        std::cout << "h = " << h
                  << "    numerical = " << numerical
                  << "    exact = " << exact
                  << "    error = " << error;

        // Calculate convergence ratio
        if (previous_error != 0.0)
        {
            double ratio = previous_error / error;

            std::cout << "    ratio = " << ratio;
        }

        std::cout << '\n';

        previous_error = error;
    }

    return 0;
}