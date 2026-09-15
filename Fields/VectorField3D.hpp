#pragma once
#include<vector>
#include<cstddef>
#include"Vector3D.hpp"
#include"ScalarField3D.hpp"

namespace Onko
{
        class VectorField3D
    {
    private:
        Onko::ScalarField3D xField, yField, zField; 
    public:
        VectorField3D(std::size_t nx, std::size_t ny, std::size_t nz):
            xField(nx, ny, nz), 
            yField(nx, ny, nz), 
            zField(nx, ny, nz)
        {
        }
        Onko::ScalarField3D& x(){
            return xField;
        }
        Onko::ScalarField3D& y(){
            return yField;
        }
        Onko::ScalarField3D& z(){
            return zField; 
        }
        const Onko::ScalarField3D& x() const{
            return xField;
        }
        const Onko::ScalarField3D& y() const{
            return yField;
        }
        const Onko::ScalarField3D& z() const{
            return zField;
        }
    };
}
