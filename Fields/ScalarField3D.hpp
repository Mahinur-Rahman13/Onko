#pragma once
#include<cstddef>
#include<vector>
namespace Onko
{
    class ScalarField3D
    {
    private:
        std::size_t Nx, Ny, Nz;
        std::vector<double> data;
    public:
        ScalarField3D(std::size_t nx, std::size_t ny, std::size_t nz):
            Nx(nx), Ny(ny), Nz(nz), data(Nx * Ny * Nz, 0.0)
        {
        }
        double& operator()(std::size_t i, std::size_t j, std::size_t k)
        {
            return data[i + Nx * (j + Ny * k)];
        }
        const double& operator()(
        std::size_t i, std::size_t j, std::size_t k) const
        {
            return data[i + Nx * (j + Ny * k)];
        }
    };

}