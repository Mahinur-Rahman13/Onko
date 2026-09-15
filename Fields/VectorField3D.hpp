#pragma once
#include<vector>
#include<cstddef>
#include"Vector3D.hpp"
namespace Onko
{
        class VectorField3D
    {
    private:
        std::size_t Nx, Ny, Nz;
        std::vector<Onko::Vector3D> data;
    public:
        VectorField3D(std::size_t nx, std::size_t ny, std::size_t nz):
            Nx(nx), Ny(ny), Nz(nz), data(Nx * Ny * Nz, {0.0, 0.0, 0.0})
        {
        }
        Onko::Vector3D& operator()(std::size_t i, std::size_t j, std::size_t k)
        {
            return data[i + Nx * (j + Ny * k)];
        }

    };
}
