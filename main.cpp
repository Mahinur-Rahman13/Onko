#include <iostream>
#include <cmath>
#include <cstddef>
#include <algorithm>

#include "Fields/ScalarField3D.hpp"
#include "Fields/VectorField3D.hpp"

#include "Operators/FirstDerivative.hpp"
#include "Operators/SecondDerivative.hpp"
#include "Operators/Laplacian.hpp"
#include "Operators/Gradient.hpp"
#include "Operators/Divergence.hpp"
#include "Operators/TemporalDerivative.hpp"
#include "Operators/VolumetricStrain.hpp"


int main()
{
    // ============================================================
    // 1. MATERIAL PARAMETERS
    // ============================================================

    // SI units

    const double mu = 4.0e8;            // Pa
    const double lambda = 2.6666667e8;  // Pa

    const double alpha = 0.952381;

    const double Q = 1.47e9;            // Pa

    const double chi = 1.53061;


    // ============================================================
    // 2. COMPUTATIONAL DOMAIN
    // ============================================================

    const double Lx = 200.0;    // m
    const double Ly = 200.0;    // m
    const double Lz = 200.0;    // m

    // ------------------------------------------------------------
    // Start with a coarse grid.
    //
    // 41 points -> 40 intervals -> 5 m spacing
    //
    // Later:
    //
    // Nx = Ny = Nz = 201
    // dx = dy = dz = 1 m
    // ------------------------------------------------------------

    const std::size_t Nx = 41;
    const std::size_t Ny = 41;
    const std::size_t Nz = 41;

    const double dx = Lx / static_cast<double>(Nx - 1);
    const double dy = Ly / static_cast<double>(Ny - 1);
    const double dz = Lz / static_cast<double>(Nz - 1);

    // For now we assume a uniform grid.
    const double h = dx;


    // ============================================================
    // 3. TIME PARAMETERS
    // ============================================================

    double t = 0.0;

    const double dt = 1.0;       // s

    // Start with 1 hour of simulation time.
    const double t_end = 3600.0; // s


    // ============================================================
    // 4. SOURCE
    // ============================================================

    // Physical source location:
    //
    // x = 100 m
    // y = 100 m
    // z = 140 m
    //
    // For a 41^3 grid with h = 5 m:
    //
    // i = 20
    // j = 20
    // k = 28
    //

    const std::size_t source_i = 20;
    const std::size_t source_j = 20;
    const std::size_t source_k = 28;

    // Injection rate:
    //
    // 32 m^3/hour
    //
    // converted to m^3/s

    const double injection_rate = 32.0 / 3600.0;


    // ============================================================
    // 5. FIELDS
    // ============================================================

    Onko::ScalarField3D p(Nx, Ny, Nz);
    Onko::ScalarField3D p_old(Nx, Ny, Nz);

    Onko::VectorField3D u(Nx, Ny, Nz);
    Onko::VectorField3D u_old(Nx, Ny, Nz);


    // ============================================================
    // 6. SOLVER PARAMETERS
    // ============================================================

    const int max_iterations = 500;

    const double tolerance = 1.0e-6;


    // ============================================================
    // 7. INITIAL CONDITIONS
    // ============================================================

    // Constructors initialize fields to zero:
    //
    // p = 0
    // u = 0
    //
    // Therefore no additional initialization is required.


    // ============================================================
    // 8. TIME LOOP
    // ============================================================

    while (t < t_end)
    {
        std::cout << "\n========================================\n";
        std::cout << "Time = " << t << " s\n";
        std::cout << "========================================\n";


        // --------------------------------------------------------
        // Store previous time-step solution
        // --------------------------------------------------------

        p_old = p;
        u_old = u;


        // ========================================================
        // COUPLED ITERATION
        // ========================================================

        for (int iteration = 0;
             iteration < max_iterations;
             ++iteration)
        {
            // ----------------------------------------------------
            // Save current iteration
            // ----------------------------------------------------

            Onko::ScalarField3D p_previous = p;
            Onko::VectorField3D u_previous = u;


            // ====================================================
            // A. PRESSURE EQUATION
            // ====================================================

            //
            // Q^(-1) p_t
            //
            // + alpha d(epsilon_v)/dt
            //
            // - chi Laplacian(p)
            //
            // = q
            //
            //
            // where
            //
            // epsilon_v = div(u)
            //
            // and
            //
            // p_t = Dt(p, p_old, dt)
            //
            // d(epsilon_v)/dt =
            // Dt(epsilon_v, epsilon_v_old, dt)
            //
            // ====================================================


            for (std::size_t k = 2; k < Nz - 2; ++k)
            {
                for (std::size_t j = 2; j < Ny - 2; ++j)
                {
                    for (std::size_t i = 2; i < Nx - 2; ++i)
                    {
                        // ----------------------------------------
                        // Volumetric strain
                        // ----------------------------------------

                        double epsilon_v =
                            VolStr(u, i, j, k, h);

                        double epsilon_v_old =
                            VolStr(u_old, i, j, k, h);


                        // ----------------------------------------
                        // Time derivatives
                        // ----------------------------------------

                        double dp_dt =
                            Onko::Dt(
                                p(i,j,k),
                                p_old(i,j,k),
                                dt
                            );

                        double d_epsilon_v_dt =
                            Onko::Dt(
                                epsilon_v,
                                epsilon_v_old,
                                dt
                            );


                        // ----------------------------------------
                        // Source
                        // ----------------------------------------

                        double q = 0.0;

                        if (i == source_i &&
                            j == source_j &&
                            k == source_k)
                        {
                            q =
                                injection_rate
                                /
                                (dx * dy * dz);
                        }


                        // ----------------------------------------
                        // Laplacian of pressure
                        // ----------------------------------------

                        double lap_p =
                            Laplacian(
                                p,
                                i,
                                j,
                                k,
                                h
                            );


                        // ----------------------------------------
                        // Pressure PDE residual
                        // ----------------------------------------

                        double pressure_residual =
                            (1.0 / Q) * dp_dt
                            +
                            alpha * d_epsilon_v_dt
                            -
                            chi * lap_p
                            -
                            q;


                        // ----------------------------------------
                        // Temporary pressure correction
                        //
                        // IMPORTANT:
                        //
                        // This is a simple relaxation prototype.
                        // We will replace this with the properly
                        // derived discrete linear solve.
                        // ----------------------------------------

                        const double pressure_relaxation =
                            0.1 * Q * dt;

                        p(i,j,k) -=
                            pressure_relaxation
                            * pressure_residual;
                    }
                }
            }


            // ====================================================
            // B. PRESSURE BOUNDARY CONDITIONS
            // ====================================================

            //
            // p = 0 on ALL six boundaries.
            //


            // x = 0 and x = Lx

            for (std::size_t j = 0; j < Ny; ++j)
            {
                for (std::size_t k = 0; k < Nz; ++k)
                {
                    p(0,j,k) = 0.0;
                    p(Nx-1,j,k) = 0.0;
                }
            }


            // y = 0 and y = Ly

            for (std::size_t i = 0; i < Nx; ++i)
            {
                for (std::size_t k = 0; k < Nz; ++k)
                {
                    p(i,0,k) = 0.0;
                    p(i,Ny-1,k) = 0.0;
                }
            }


            // z = 0 and z = Lz

            for (std::size_t i = 0; i < Nx; ++i)
            {
                for (std::size_t j = 0; j < Ny; ++j)
                {
                    p(i,j,0) = 0.0;
                    p(i,j,Nz-1) = 0.0;
                }
            }


            // ====================================================
            // C. MECHANICAL EQUILIBRIUM
            // ====================================================

            //
            // Mechanical equation:
            //
            // (lambda + mu) grad(div u)
            //
            // + mu Laplacian(u)
            //
            // - alpha grad(p)
            //
            // = 0
            //
            // ====================================================

            //
            // At this stage we evaluate the mechanical residual.
            //
            // We are deliberately NOT hiding an arbitrary
            // relaxation coefficient here.
            //
            // The proper discrete mechanical system will be
            // derived and solved next.
            //


            double maximum_mechanical_residual = 0.0;


            for (std::size_t k = 2; k < Nz - 2; ++k)
            {
                for (std::size_t j = 2; j < Ny - 2; ++j)
                {
                    for (std::size_t i = 2; i < Nx - 2; ++i)
                    {
                        // ----------------------------------------
                        // Pressure gradient
                        // ----------------------------------------

                        double dp_dx =
                            Dx(p, i, j, k, dx);

                        double dp_dy =
                            Dy(p, i, j, k, dy);

                        double dp_dz =
                            Dz(p, i, j, k, dz);


                        // ----------------------------------------
                        // Displacement Laplacians
                        // ----------------------------------------

                        double lap_ux =
                            Dxx(u.x(), i,j,k,dx)
                            +
                            Dyy(u.x(), i,j,k,dy)
                            +
                            Dzz(u.x(), i,j,k,dz);

                        double lap_uy =
                            Dxx(u.y(), i,j,k,dx)
                            +
                            Dyy(u.y(), i,j,k,dy)
                            +
                            Dzz(u.y(), i,j,k,dz);

                        double lap_uz =
                            Dxx(u.z(), i,j,k,dx)
                            +
                            Dyy(u.z(), i,j,k,dy)
                            +
                            Dzz(u.z(), i,j,k,dz);


                        // ----------------------------------------
                        // Divergence
                        // ----------------------------------------

                        double div_u =
                            Divergence(
                                u,
                                i,j,k,
                                h
                            );


                        // ----------------------------------------
                        // grad(div u)
                        //
                        // For now use second-order central
                        // differences for these mixed terms.
                        //
                        // This will be upgraded to a consistent
                        // fourth-order formulation.
                        // ----------------------------------------

                        double d2ux_dx2 =
                            Dxx(u.x(), i,j,k,dx);

                        double d2uy_dy2 =
                            Dyy(u.y(), i,j,k,dy);

                        double d2uz_dz2 =
                            Dzz(u.z(), i,j,k,dz);


                        double d2uy_dxdy =
                            (
                                Dy(u.y(), i+1,j,k,dy)
                                -
                                Dy(u.y(), i-1,j,k,dy)
                            )
                            /
                            (2.0 * dx);


                        double d2uz_dxdz =
                            (
                                Dz(u.z(), i+1,j,k,dz)
                                -
                                Dz(u.z(), i-1,j,k,dz)
                            )
                            /
                            (2.0 * dx);


                        double d2ux_dxdy =
                            (
                                Dx(u.x(), i,j+1,k,dx)
                                -
                                Dx(u.x(), i,j-1,k,dx)
                            )
                            /
                            (2.0 * dy);


                        double d2uz_dydz =
                            (
                                Dz(u.z(), i,j+1,k,dz)
                                -
                                Dz(u.z(), i,j-1,k,dz)
                            )
                            /
                            (2.0 * dy);


                        double d2ux_dxdz =
                            (
                                Dx(u.x(), i,j,k+1,dx)
                                -
                                Dx(u.x(), i,j,k-1,dx)
                            )
                            /
                            (2.0 * dz);


                        double d2uy_dydz =
                            (
                                Dy(u.y(), i,j,k+1,dy)
                                -
                                Dy(u.y(), i,j,k-1,dy)
                            )
                            /
                            (2.0 * dz);


                        // ----------------------------------------
                        // grad(div u)
                        // ----------------------------------------

                        double grad_div_x =
                            d2ux_dx2
                            +
                            d2uy_dxdy
                            +
                            d2uz_dxdz;


                        double grad_div_y =
                            d2ux_dxdy
                            +
                            d2uy_dy2
                            +
                            d2uz_dydz;


                        double grad_div_z =
                            d2ux_dxdz
                            +
                            d2uy_dydz
                            +
                            d2uz_dz2;


                        // ----------------------------------------
                        // Mechanical residual
                        // ----------------------------------------

                        double Rx =
                            (lambda + mu)
                            * grad_div_x
                            +
                            mu * lap_ux
                            -
                            alpha * dp_dx;


                        double Ry =
                            (lambda + mu)
                            * grad_div_y
                            +
                            mu * lap_uy
                            -
                            alpha * dp_dy;


                        double Rz =
                            (lambda + mu)
                            * grad_div_z
                            +
                            mu * lap_uz
                            -
                            alpha * dp_dz;


                        double magnitude =
                            std::sqrt(
                                Rx*Rx +
                                Ry*Ry +
                                Rz*Rz
                            );


                        maximum_mechanical_residual =
                            std::max(
                                maximum_mechanical_residual,
                                magnitude
                            );
                    }
                }
            }


            // ====================================================
            // D. MECHANICAL BOUNDARY CONDITIONS
            // ====================================================

            //
            // u = 0 on:
            //
            // x = 0
            // x = Lx
            // y = 0
            // y = Ly
            // z = 0
            //
            // Top surface:
            //
            // traction-free
            //
            // These are currently applied only to the fixed
            // displacement boundaries.
            //
            // The traction-free top boundary will be implemented
            // explicitly after the interior mechanical equation
            // is verified.
            //


            // x boundaries

            for (std::size_t j = 0; j < Ny; ++j)
            {
                for (std::size_t k = 0; k < Nz; ++k)
                {
                    u.x()(0,j,k) = 0.0;
                    u.y()(0,j,k) = 0.0;
                    u.z()(0,j,k) = 0.0;

                    u.x()(Nx-1,j,k) = 0.0;
                    u.y()(Nx-1,j,k) = 0.0;
                    u.z()(Nx-1,j,k) = 0.0;
                }
            }


            // y boundaries

            for (std::size_t i = 0; i < Nx; ++i)
            {
                for (std::size_t k = 0; k < Nz; ++k)
                {
                    u.x()(i,0,k) = 0.0;
                    u.y()(i,0,k) = 0.0;
                    u.z()(i,0,k) = 0.0;

                    u.x()(i,Ny-1,k) = 0.0;
                    u.y()(i,Ny-1,k) = 0.0;
                    u.z()(i,Ny-1,k) = 0.0;
                }
            }


            // bottom boundary

            for (std::size_t i = 0; i < Nx; ++i)
            {
                for (std::size_t j = 0; j < Ny; ++j)
                {
                    u.x()(i,j,0) = 0.0;
                    u.y()(i,j,0) = 0.0;
                    u.z()(i,j,0) = 0.0;
                }
            }


            // ====================================================
            // E. ITERATION CONVERGENCE
            // ====================================================

            double p_difference = 0.0;
            double p_norm = 0.0;

            double u_difference = 0.0;
            double u_norm = 0.0;


            for (std::size_t k = 0; k < Nz; ++k)
            {
                for (std::size_t j = 0; j < Ny; ++j)
                {
                    for (std::size_t i = 0; i < Nx; ++i)
                    {
                        // Pressure

                        double dp =
                            p(i,j,k)
                            -
                            p_previous(i,j,k);

                        p_difference += dp * dp;

                        p_norm +=
                            p(i,j,k)
                            *
                            p(i,j,k);


                        // Displacement

                        double dux =
                            u.x()(i,j,k)
                            -
                            u_previous.x()(i,j,k);

                        double duy =
                            u.y()(i,j,k)
                            -
                            u_previous.y()(i,j,k);

                        double duz =
                            u.z()(i,j,k)
                            -
                            u_previous.z()(i,j,k);


                        u_difference +=
                            dux*dux +
                            duy*duy +
                            duz*duz;


                        u_norm +=
                            u.x()(i,j,k)
                            *
                            u.x()(i,j,k)
                            +
                            u.y()(i,j,k)
                            *
                            u.y()(i,j,k)
                            +
                            u.z()(i,j,k)
                            *
                            u.z()(i,j,k);
                    }
                }
            }


            p_difference =
                std::sqrt(p_difference);

            p_norm =
                std::sqrt(p_norm);

            u_difference =
                std::sqrt(u_difference);

            u_norm =
                std::sqrt(u_norm);


            double pressure_error =
                p_difference /
                std::max(p_norm, 1.0);


            double displacement_error =
                u_difference /
                std::max(u_norm, 1.0);


            // ----------------------------------------------------
            // Print progress
            // ----------------------------------------------------

            if (iteration % 10 == 0)
            {
                std::cout
                    << "Iteration "
                    << iteration
                    << " | dp = "
                    << pressure_error
                    << " | du = "
                    << displacement_error
                    << " | mechanical residual = "
                    << maximum_mechanical_residual
                    << '\n';
            }


            // ----------------------------------------------------
            // Coupled convergence
            // ----------------------------------------------------

            if (pressure_error < tolerance &&
                displacement_error < tolerance)
            {
                std::cout
                    << "Coupled iteration converged after "
                    << iteration
                    << " iterations.\n";

                break;
            }
        }


        // ========================================================
        // ADVANCE TIME
        // ========================================================

        t += dt;
    }


    // ============================================================
    // 9. FINAL RESULTS
    // ============================================================

    std::cout
        << "\n========================================\n";

    std::cout
        << "Simulation finished.\n";

    std::cout
        << "Final time = "
        << t
        << " s\n";


    std::cout
        << "Pressure at source = "
        << p(source_i, source_j, source_k)
        << " Pa\n";


    std::cout
        << "Displacement at source = ("
        << u.x()(source_i, source_j, source_k)
        << ", "
        << u.y()(source_i, source_j, source_k)
        << ", "
        << u.z()(source_i, source_j, source_k)
        << ") m\n";


    return 0;
}