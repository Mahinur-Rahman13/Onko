#pragma once
#include <cstddef>
#include "../Fields/VectorField3D.hpp"
#include "FirstDerivative.hpp"

namespace Onko{
    inline double Divergence(
        const Onko::VectorField3D& u, 
        std::size_t i, 
        std::size_t j, 
        std::size_t k, 
        double h
    )
    {
        return Dx(u.x(), i, j, k, h) + Dy(u.y(), i, j, k, h) + Dz(u.z(), i, j, k, h);
    }
}