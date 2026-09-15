//This is still an interior point operator
#pragma once
#include<cstddef>
#include<cmath>
#include "../Fields/ScalarField3D.hpp"
namespace Onko
{
    inline double Dxx(
        const ScalarField3D& f,
        std::size_t i,
        std::size_t j, 
        std::size_t k,
        double dx
    )
    {
        return (-f(i + 2, j, k) + 16.0 * f(i + 1, j, k) - 30.0 * f(i, j, k) + 16.0 * f(i - 1, j, k) - f(i - 2, j, k)) / (12.0 * std::pow(dx, 2.0));
    }
    inline double Dyy(
        const ScalarField3D& f,
        std::size_t i, 
        std::size_t j,
        std::size_t k,
        double dy
    )
    {
        return (-f(i, j + 2, k) + 16.0 * f(i, j + 1, k) - 30.0 * f(i, j, k) + 16.0 * f(i, j - 1, k) - f(i, j - 2, k)) / (12.0 * std::pow(dy, 2.0));
    }
    inline double Dzz(
        const ScalarField3D& f,
        std::size_t i, 
        std::size_t j, 
        std::size_t k,
        double dz
    )
    {
        return (-f(i, j, k + 2) + 16.0 * f(i, j, k + 1) - 30.0 * f(i, j, k) + 16.0 * f(i, j, k - 1) - f(i, j, k - 2)) / (12.0 * std::pow(dz, 2.0));
    }
}