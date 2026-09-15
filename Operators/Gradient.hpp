#pragma once
#include "../Fields/Vector3D.hpp"
#include "../Fields/ScalarField3D.hpp"
#include "FirstDerivative.hpp"

namespace Onko
{
    inline Vector3D Gradient(
        const ScalarField3D& f,
        std::size_t i,
        std::size_t j,
        std::size_t k,
        double h
    )
    {
        Vector3D result;
        result.x = Dx(f, i, j, k, h);
        result.y = Dy(f, i, j, k ,h);
        result.z = Dz(f, i, j, k, h);
        return result;

    }
}