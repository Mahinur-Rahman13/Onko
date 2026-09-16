#pragma once
#include "../Fields/VectorField3D.hpp"
#include "Divergence.hpp"

namespace Onko
{
    inline double VolStr(
        const Onko::VectorField3D& u, 
        std::size_t i, 
        std::size_t j, 
        std::size_t k, 
        double h
    )
    {
        return Divergence(u, i, j, k, h);
    }
}