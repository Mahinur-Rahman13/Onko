//This is still and interior point operator
#pragma once
#include "../Fields/ScalarField3D.hpp"
#include "SecondDerivative.hpp"
#include<cstddef>

namespace Onko
{
    inline double Laplacian
    (
        const ScalarField3D& f, 
        std::size_t i,
        std::size_t j, 
        std::size_t k,
        double h
    ) 
    {
        return Dxx(f, i , j, k, h) + Dyy(f, i, j, k, h) + Dzz(f, i, j, k, h);
    }
}