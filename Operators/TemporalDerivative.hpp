//Backward-Euler method
#pragma once 
namespace Onko
{
    inline double Dt(
        double current, 
        double previous, 
        double dt
    )
    {
        return (current - previous) / dt;
    }
}