#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstddef>

#include "Fields/ScalarField3D.hpp"
#include "Operators/FirstDerivative.hpp"
#include "Operators/SecondDerivative.hpp"

int main()
{
    std::cout << std::fixed << std::setprecision(10);

    // Point at which we evaluate the derivatives.
    const double x = 0.7;
    const double y = 0.8;
    const double z = 0.9;

    // Grid spacings used for the convergence test.
    const double spacings[] = {0.2, 0.1, 0.05, 0.025};

    // Five points are required because our fourth-order
    // central stencil uses i-2 through i+2.
    const std::size_t Nx = 5;
    const std::size_t Ny = 5;
    const std::size_t Nz = 5;

    // The point we are testing is always the center point.
    const std::size_t i = 2;
    const std::size_t j = 2;
    const std::size_t k = 2;


    // ============================================================
    // FIRST DERIVATIVES
    // ============================================================

    std::cout << "============================================\n";
    std::cout << " FOURTH-ORDER FIRST DERIVATIVE TESTS\n";
    std::cout << "============================================\n\n";


    // ------------------------------------------------------------
    // Dx
    // ------------------------------------------------------------

    std::cout << "--- Dx test ---\n";

    double previous_error = 0.0;

    for (double h : spacings)
    {
        Onko::ScalarField3D f(Nx, Ny, Nz);

        // Fill the 3D field:
        //
        // f(x,y,z) = sin(x) + sin(y) + sin(z)
        //
        // The center point is (x,y,z).

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

        double numerical = Onko::Dx(f, i, j, k, h);

        // Exact derivative:
        // d/dx [sin(x) + sin(y) + sin(z)] = cos(x)
        double exact = std::cos(x);

        double error = std::abs(numerical - exact);

        std::cout << "h = " << h
                  << "    numerical = " << numerical
                  << "    exact = " << exact
                  << "    error = " << error;

        if (previous_error != 0.0)
        {
            double ratio = previous_error / error;
            std::cout << "    ratio = " << ratio;
        }

        std::cout << '\n';

        previous_error = error;
    }


    // ------------------------------------------------------------
    // Dy
    // ------------------------------------------------------------

    std::cout << "\n--- Dy test ---\n";

    previous_error = 0.0;

    for (double h : spacings)
    {
        Onko::ScalarField3D f(Nx, Ny, Nz);

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

        double numerical = Onko::Dy(f, i, j, k, h);

        // Exact derivative with respect to y.
        double exact = std::cos(y);

        double error = std::abs(numerical - exact);

        std::cout << "h = " << h
                  << "    numerical = " << numerical
                  << "    exact = " << exact
                  << "    error = " << error;

        if (previous_error != 0.0)
        {
            double ratio = previous_error / error;
            std::cout << "    ratio = " << ratio;
        }

        std::cout << '\n';

        previous_error = error;
    }


    // ------------------------------------------------------------
    // Dz
    // ------------------------------------------------------------

    std::cout << "\n--- Dz test ---\n";

    previous_error = 0.0;

    for (double h : spacings)
    {
        Onko::ScalarField3D f(Nx, Ny, Nz);

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

        double numerical = Onko::Dz(f, i, j, k, h);

        // Exact derivative with respect to z.
        double exact = std::cos(z);

        double error = std::abs(numerical - exact);

        std::cout << "h = " << h
                  << "    numerical = " << numerical
                  << "    exact = " << exact
                  << "    error = " << error;

        if (previous_error != 0.0)
        {
            double ratio = previous_error / error;
            std::cout << "    ratio = " << ratio;
        }

        std::cout << '\n';

        previous_error = error;
    }


    // ============================================================
    // SECOND DERIVATIVES
    // ============================================================

    std::cout << "\n============================================\n";
    std::cout << " FOURTH-ORDER SECOND DERIVATIVE TESTS\n";
    std::cout << "============================================\n\n";


    // ------------------------------------------------------------
    // Dxx
    // ------------------------------------------------------------

    std::cout << "--- Dxx test ---\n";

    previous_error = 0.0;

    for (double h : spacings)
    {
        Onko::ScalarField3D f(Nx, Ny, Nz);

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

        double numerical = Onko::Dxx(f, i, j, k, h);

        // Exact second derivative with respect to x.
        double exact = -std::sin(x);

        double error = std::abs(numerical - exact);

        std::cout << "h = " << h
                  << "    numerical = " << numerical
                  << "    exact = " << exact
                  << "    error = " << error;

        if (previous_error != 0.0)
        {
            double ratio = previous_error / error;
            std::cout << "    ratio = " << ratio;
        }

        std::cout << '\n';

        previous_error = error;
    }


    // ------------------------------------------------------------
    // Dyy
    // ------------------------------------------------------------

    std::cout << "\n--- Dyy test ---\n";

    previous_error = 0.0;

    for (double h : spacings)
    {
        Onko::ScalarField3D f(Nx, Ny, Nz);

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

        double numerical = Onko::Dyy(f, i, j, k, h);

        // Exact second derivative with respect to y.
        double exact = -std::sin(y);

        double error = std::abs(numerical - exact);

        std::cout << "h = " << h
                  << "    numerical = " << numerical
                  << "    exact = " << exact
                  << "    error = " << error;

        if (previous_error != 0.0)
        {
            double ratio = previous_error / error;
            std::cout << "    ratio = " << ratio;
        }

        std::cout << '\n';

        previous_error = error;
    }


    // ------------------------------------------------------------
    // Dzz
    // ------------------------------------------------------------

    std::cout << "\n--- Dzz test ---\n";

    previous_error = 0.0;

    for (double h : spacings)
    {
        Onko::ScalarField3D f(Nx, Ny, Nz);

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

        double numerical = Onko::Dzz(f, i, j, k, h);

        // Exact second derivative with respect to z.
        double exact = -std::sin(z);

        double error = std::abs(numerical - exact);

        std::cout << "h = " << h
                  << "    numerical = " << numerical
                  << "    exact = " << exact
                  << "    error = " << error;

        if (previous_error != 0.0)
        {
            double ratio = previous_error / error;
            std::cout << "    ratio = " << ratio;
        }

        std::cout << '\n';

        previous_error = error;
    }


    std::cout << "\nAll tests completed.\n";

    return 0;
}