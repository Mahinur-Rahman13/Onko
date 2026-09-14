#include <iostream>
#include "VectorField3D.hpp"

int main()
{
    std::size_t Nx = 10;
    std::size_t Ny = 20;
    std::size_t Nz = 30;

    VectorField3D displacement(Nx, Ny, Nz);

    displacement(2, 3, 4).x = 1.0;
    displacement(2, 3, 4).y = 2.0;
    displacement(2, 3, 4).z = 3.0;

    std::cout << "ux = "
              << displacement(2, 3, 4).x << '\n';

    std::cout << "uy = "
              << displacement(2, 3, 4).y << '\n';

    std::cout << "uz = "
              << displacement(2, 3, 4).z << '\n';

    return 0;
}
