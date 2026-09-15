//This is still an interior point operator. Works on Ghost nodes or boundary reworks arrives soon. 
#pragma once
#include<cstddef>
#include "../Fields/ScalarField3D.hpp"
namespace Onko
{
    inline double Dx(
        const ScalarField3D& f, 
        std::size_t i,
        std::size_t j,
        std::size_t k,
        double dx
    )
    {
        return(-f(i + 2, j, k) + 8.0 * f(i + 1, j, k) - 8.0 * f(i - 1, j, k) + f(i - 2, j, k)) / (12.0 * dx);
    }
    inline double Dy(
        const ScalarField3D& f, 
        std::size_t i,
        std::size_t j,
        std::size_t k,
        double dy
    )
    {
        return(-f(i, j + 2, k) + 8.0 * f(i, j + 1, k) - 8.0 * f(i, j - 1, k) + f(i, j - 2, k)) / (12.0 * dy);
    }
    inline double Dz(
        const ScalarField3D& f, 
        std::size_t i,
        std::size_t j,
        std::size_t k,
        double dz
    )
    {
        return(-f(i, j, k + 2) + 8.0 * f(i, j, k + 1) - 8.0 * f(i, j, k - 1) + f(i, j, k - 2)) / (12.0 * dz);
    }
}