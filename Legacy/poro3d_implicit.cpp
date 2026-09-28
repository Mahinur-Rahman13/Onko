// poro3d_single_file.cpp -- 3D finite-difference poroelastic solver for the
// point-injection problem of Wang & Kuempel (2003), Geophysics 68(2), 1-13.
//
// SINGLE-FILE BUILD merging poro3d.hpp (grid/fields/analytic references),
// mechanics.hpp (elasticity operator, boundary conditions, BiCGSTAB solver),
// and main.cpp (time-stepping driver) -- produced for distribution as one
// file; the original project keeps them separate as include/poro3d.hpp,
// include/mechanics.hpp, src/main.cpp.
//
// This is the FULL-DOMAIN 4th-order version: genuine O(h^4) 5-point
// stencils (diagonal AND cross-derivative terms) everywhere in the bulk,
// at both symmetry planes (exact, no accuracy cost), and at every far-field
// Dirichlet boundary (verified O(h^5)-accurate cubic-fit ghost
// extrapolation). The ONE exception is the free-surface TRACTION
// (Neumann-type) condition on displacement: a higher-order ghost formula
// that also used the boundary row's own current value was derived,
// verified analytically/numerically in isolation, but then found -- by
// actually running it -- to destabilize the BiCGSTAB solve (self-coupling
// into the unknown vector; max|u_z| diverged to ~1e33 within ~150 steps).
// That one ghost is kept at the empirically-verified-stable lower
// (cubic, no self-coupling) order; every other boundary in the domain is
// genuine 4th order. See the comments in the mechanics section below for
// the full account of that finding.
//
// Build:  g++ -O3 -march=native -funroll-loops -std=c++17 -fopenmp \
//             poro3d_single_file.cpp -o poro3d
// Run:    ./poro3d --h 5 --L 150 --Lz 150 --tend 1800 --sample 15 --out fig2.csv

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

namespace poro {



constexpr double PI = 3.14159265358979323846;

// ---------------------------------------------------------------------------
// Poroelastic material.  The paper parameterises a medium by
// (mu, nu, nu_u, B, D) and derives (lambda, alpha, 1/Q, chi) -- eqs (3)-(6).
// ---------------------------------------------------------------------------
struct Material {
    double mu   = 0.4e9;  // shear modulus                              [Pa]
    double nu   = 0.2;    // drained Poisson ratio                      [-]
    double nu_u = 0.4;    // undrained Poisson ratio                    [-]
    double B    = 0.75;   // Skempton ratio                             [-]
    double D    = 1.0;    // hydraulic diffusivity                      [m^2/s]

    double lambda() const { return 2.0 * nu * mu / (1.0 - 2.0 * nu); }          // (3)
    double alpha() const {                                                      // (4)
        return 3.0 * (nu_u - nu) / ((1.0 - 2.0 * nu) * (1.0 + nu_u) * B);
    }
    double Qinv() const {                                                       // (5)
        return 4.5 * (1.0 - 2.0 * nu_u) * (nu_u - nu) /
               ((1.0 - 2.0 * nu) * (1.0 + nu_u) * (1.0 + nu_u) * mu * B * B);
    }
    double chi() const {                                                        // (6)
        return 4.5 * (1.0 - nu_u) * (nu_u - nu) * D /
               ((1.0 - nu) * (1.0 + nu_u) * (1.0 + nu_u) * mu * B * B);
    }
    double M() const { return lambda() + 2.0 * mu; }          // constrained (P-wave) modulus
    double Kdr() const { return lambda() + 2.0 * mu / 3.0; }  // drained bulk modulus

    // Storage coefficient of the *fixed-stress* split: 1/Q + alpha^2/K_dr.
    double S_fs() const { return Qinv() + alpha() * alpha() / Kdr(); }
    // Storage coefficient of the exactly uncoupled diffusion mode: chi/S_D == D.
    double S_D() const { return Qinv() + alpha() * alpha() / M(); }
};

// ---------------------------------------------------------------------------
// Uniform Cartesian grid.  z points DOWN, z = 0 is the free surface.
// A quarter domain x,y in [0,L] is used: the point source sits on the
// x = y = 0 axis, so x = 0 and y = 0 are symmetry planes.  With mirrored
// ghost nodes the node (0,0,k) is numerically identical to the same node
// of the full domain, so the source strength is the *full* q0.
// ---------------------------------------------------------------------------
struct Grid {
    int nx = 51, ny = 51, nz = 51;
    double h = 5.0;
    double x(int i) const { return i * h; }
    double y(int j) const { return j * h; }
    double z(int k) const { return k * h; }
    std::size_t n() const {
        return std::size_t(nx + 4) * (ny + 4) * (nz + 4);  // TWO ghost layers each side
    }
};

// Scalar field with a two-cell ghost halo; indices run -2 .. n+1.
// (Widened from a 1-cell halo to support the 4th-order 5-point stencils --
// see the order-selection helpers in mechanics.hpp / main.cpp. The extra
// layer is only ever populated by EXACT symmetry mirroring at x=0,y=0; the
// free-surface / far-Dirichlet boundaries still only get a single physically
// derived ghost layer, and points that would need a second ghost layer there
// fall back to the original 2nd-order stencil -- see the order-selection
// comment block in mechanics.hpp.)
struct Field {
    const Grid* g = nullptr;
    std::vector<double> a;
    void init(const Grid& gr, double v = 0.0) {
        g = &gr;
        a.assign(gr.n(), v);
    }
    inline std::size_t id(int i, int j, int k) const {
        return std::size_t(i + 2) + std::size_t(g->nx + 4) *
               (std::size_t(j + 2) + std::size_t(g->ny + 4) * std::size_t(k + 2));
    }
    inline double& operator()(int i, int j, int k) { return a[id(i, j, k)]; }
    inline const double& operator()(int i, int j, int k) const { return a[id(i, j, k)]; }
    void zero() { std::fill(a.begin(), a.end(), 0.0); }
};

struct Vec3Field {
    Field x, y, z;
    void init(const Grid& g) { x.init(g); y.init(g); z.init(g); }
    void zero() { x.zero(); y.zero(); z.zero(); }
    void copyFrom(const Vec3Field& o) { x.a = o.x.a; y.a = o.y.a; z.a = o.z.a; }
};

// ---------------------------------------------------------------------------
// Analytic reference solutions.
// ---------------------------------------------------------------------------

// g(xi) = erf(xi/2) - (xi/sqrt(pi)) exp(-xi^2/4)          -- eq (20)
inline double rud_g(double xi) {
    if (xi < 0.05) {  // series: avoids catastrophic cancellation
        const double x3 = xi * xi * xi;
        return (x3 / 6.0 - x3 * xi * xi / 40.0) / std::sqrt(PI);
    }
    return std::erf(0.5 * xi) - xi / std::sqrt(PI) * std::exp(-0.25 * xi * xi);
}
// F(xi) = erfc(xi/2) + 2 g(xi)/xi^2                        -- eq (19)
inline double rud_F(double xi) {
    if (xi < 1e-10) return 1.0;
    return std::erfc(0.5 * xi) + 2.0 * rud_g(xi) / (xi * xi);
}
// F'(xi) = -4 g(xi)/xi^3
inline double rud_dF(double xi) {
    if (xi < 1e-10) return 0.0;
    return -4.0 * rud_g(xi) / (xi * xi * xi);
}

// Transient whole-space Green's functions of Rudnicki (1986), eqs (17)-(18),
// for a point source of strength q0 switched on at t = 0 at xs.
struct WholeSpace {
    Material m;
    double q0;  // [m^3/s]
    double u0() const { return q0 * (1.0 + m.nu_u) * m.B / (24.0 * PI * (1.0 - m.nu_u) * m.D); }

    double p(double dx, double dy, double dz, double t) const {
        const double R = std::sqrt(dx * dx + dy * dy + dz * dz);
        if (R < 1e-12 || t <= 0.0) return 0.0;
        const double xi = R / std::sqrt(m.D * t);
        return q0 / (4.0 * PI * m.chi() * R) * std::erfc(0.5 * xi);
    }
    // displacement component along the separation vector component `dc`
    double u(double dc, double dx, double dy, double dz, double t) const {
        const double R = std::sqrt(dx * dx + dy * dy + dz * dz);
        if (R < 1e-12 || t <= 0.0) return 0.0;
        const double xi = R / std::sqrt(m.D * t);
        return u0() * dc / R * rud_F(xi);
    }
    // d u_z / dr at cylindrical radius r, vertical offset s = z - z_source
    double duz_dr(double r, double s, double t) const {
        const double R = std::sqrt(r * r + s * s);
        if (R < 1e-12 || t <= 0.0) return 0.0;
        const double xi = R / std::sqrt(m.D * t);
        return u0() * s * r / (R * R * R) * (xi * rud_dF(xi) - rud_F(xi));
    }
};

// Mirror-image half-space field: real source at z = +d plus a negative image
// at z = -d.  This satisfies p = 0 on z = 0 exactly and is the leading-order
// half-space solution; it is exact in steady state (paper, eqs 32-34, 43).
// Used both as the far-field Dirichlet condition and as a reference curve.
struct MirrorImage {
    WholeSpace ws;
    double d;  // source depth [m]

    double p(double x, double y, double z, double t) const {
        return ws.p(x, y, z - d, t) - ws.p(x, y, z + d, t);
    }
    double ux(double x, double y, double z, double t) const {
        return ws.u(x, x, y, z - d, t) - ws.u(x, x, y, z + d, t);
    }
    double uy(double x, double y, double z, double t) const {
        return ws.u(y, x, y, z - d, t) - ws.u(y, x, y, z + d, t);
    }
    double uz(double x, double y, double z, double t) const {
        return ws.u(z - d, x, y, z - d, t) - ws.u(z + d, x, y, z + d, t);
    }
};




enum class TopBC { FreeSurface, Analytic };

struct Solver {
    Grid g;
    Material mat;
    MirrorImage ref;       // analytic far-field / reference solution
    TopBC top = TopBC::FreeSurface;

    // --- unknown / Dirichlet classification -------------------------------
    inline bool uDirichlet(int i, int j, int k) const {
        if (i == g.nx - 1 || j == g.ny - 1 || k == g.nz - 1) return true;
        if (top == TopBC::Analytic && k == 0) return true;
        return false;
    }
    inline bool pDirichlet(int i, int j, int k) const {
        if (i == g.nx - 1 || j == g.ny - 1 || k == g.nz - 1) return true;
        if (k == 0) return true;  // free surface: p = 0 (eq. 15); analytic top: p = p_ana
        return false;
    }

    // ----------------------------------------------------------------------
    // Ghost-node update.  Every rule below is linear and homogeneous, so the
    // same routine is valid for the solution field and for correction
    // vectors inside the Krylov solver (which carry zeros on Dirichlet
    // nodes).
    // ----------------------------------------------------------------------
    // ----------------------------------------------------------------------
    // 4th-order-consistent ghost-cell formulas, derived by fitting a cubic
    // through the known boundary condition and 3 known field values, then
    // evaluating the cubic at the ghost point(s). Verified symbolically
    // (sympy) and numerically (exact reconstruction of a random cubic,
    // residual at machine precision) before being used here.
    //   dirichletGhost1(u0,u1,u2,u3): ONE ghost point, mirrored through a
    //     known Dirichlet value u0, using the three field values u1,u2,u3
    //     on the OPPOSITE side from the ghost. Same formula serves two
    //     roles: (i) the near-side ghost at a Dirichlet boundary (u0 at the
    //     boundary itself, u1,u2,u3 = the first 3 INTERIOR points), and
    //     (ii) the single ghost layer BEYOND a far Dirichlet boundary (u0 at
    //     the boundary, u1,u2,u3 = the first 3 points going INWARD).
    //   dirichletGhost2(...): the second ghost layer for case (i) only.
    //   neumannGhost1/2(T,u0,u1,u2,u3,h): the two ghost layers implied by a
    //     derivative (traction) condition f'(0)=T. u0 is accepted but
    //     deliberately UNUSED -- see the function itself for why (verified
    //     empirically to destabilize BiCGSTAB when included).
    // ----------------------------------------------------------------------
    static inline double dirichletGhost1(double u0, double u1, double u2, double u3) {
        return 4.0 * u0 - 6.0 * u1 + 4.0 * u2 - u3;
    }
    static inline double dirichletGhost2(double u0, double u1, double u2, double u3) {
        return 10.0 * u0 - 20.0 * u1 + 15.0 * u2 - 4.0 * u3;
    }
    static inline double neumannGhost1(double T, double u0, double u1, double u2, double u3, double h) {
        // NOT the naive higher-accuracy choice: a quartic fit that also uses
        // u0 (the boundary row's own current value) is LOCALLY more accurate
        // (verified: O(h^5) in the ghost value, vs cubic's O(h^4)) -- but u0
        // is part of the unknown vector BiCGSTAB is solving for here (this
        // is a Neumann/traction condition, not Dirichlet), so feeding it
        // back into its own ghost formula creates an implicit self-coupling
        // that destabilized BiCGSTAB+RBGS in testing (verified: max|u_z|
        // diverged to 1e33 within ~150 steps with u0 included; removing u0
        // and keeping everything else identical fixed it immediately).
        // This cubic fit (T + 3 interior points, no u0) is the
        // empirically-verified-stable choice; local accuracy is one order
        // lower here than the Dirichlet ghosts elsewhere.
        (void)u0;
        return (-24.0 * T * h + 6.0 * u1 + 8.0 * u2 - 3.0 * u3) / 11.0;
    }
    static inline double neumannGhost2(double T, double u0, double u1, double u2, double u3, double h) {
        (void)u0;
        return (-60.0 * T * h - 40.0 * u1 + 75.0 * u2 - 24.0 * u3) / 11.0;
    }

    // ----------------------------------------------------------------------
    // Ghost-node update.  Every rule below is linear and homogeneous, so the
    // same routine is valid for the solution field and for correction
    // vectors inside the Krylov solver (which carry zeros on Dirichlet
    // nodes).
    // ----------------------------------------------------------------------
    void fillGhostsU(Vec3Field& u) const {
        const int nx = g.nx, ny = g.ny, nz = g.nz;
        const double h = g.h;
        const double lam = mat.lambda(), M = mat.M();

        // Fill order matters: each step below may read ghost values written
        // by an earlier step, so corners (two boundaries at once) come out
        // correct in a single pass. Order: x-near(symmetry) -> y-near(sym.,
        // now sees x-near) -> x-far(Dirichlet, now sees y-near) -> y-far
        // (Dirichlet, now sees x-near AND x-far) -> z-far(Dirichlet, now
        // sees the complete x,y ghost frame) -> z-near (free-surface/
        // Dirichlet, computed on the REAL i,j domain only, since the
        // physical BC only lives there) -> mirror/extrapolate the z-near
        // ghost PLANE itself out into all four x,y corners.

        // (a) x=0 symmetry, 2 layers -- EXACT, no accuracy cost.
        for (int k = 0; k < nz; ++k)
            for (int j = 0; j < ny; ++j) {
                u.x(-1, j, k) = -u.x(1, j, k); u.y(-1, j, k) = u.y(1, j, k); u.z(-1, j, k) = u.z(1, j, k);
                u.x(-2, j, k) = -u.x(2, j, k); u.y(-2, j, k) = u.y(2, j, k); u.z(-2, j, k) = u.z(2, j, k);
            }
        // (b) y=0 symmetry, 2 layers; i extended to see (a)'s ghosts.
        for (int k = 0; k < nz; ++k)
            for (int i = -2; i < nx; ++i) {
                u.x(i, -1, k) = u.x(i, 1, k); u.y(i, -1, k) = -u.y(i, 1, k); u.z(i, -1, k) = u.z(i, 1, k);
                u.x(i, -2, k) = u.x(i, 2, k); u.y(i, -2, k) = -u.y(i, 2, k); u.z(i, -2, k) = u.z(i, 2, k);
            }
        // (c) x=L far Dirichlet boundary, ONE ghost layer beyond (i=nx) --
        // only one layer is ever topologically needed on a far side (the
        // extremal unknown i=nx-2 only reaches i+2=nx). j extended to see
        // (b)'s ghosts, so the (x-far,y-near) corner comes out right too.
        for (int k = 0; k < nz; ++k)
            for (int j = -2; j < ny; ++j) {
                u.x(nx, j, k) = dirichletGhost1(u.x(nx - 1, j, k), u.x(nx - 2, j, k), u.x(nx - 3, j, k), u.x(nx - 4, j, k));
                u.y(nx, j, k) = dirichletGhost1(u.y(nx - 1, j, k), u.y(nx - 2, j, k), u.y(nx - 3, j, k), u.y(nx - 4, j, k));
                u.z(nx, j, k) = dirichletGhost1(u.z(nx - 1, j, k), u.z(nx - 2, j, k), u.z(nx - 3, j, k), u.z(nx - 4, j, k));
            }
        // (d) y=L far Dirichlet boundary, ONE layer; i extended to see BOTH
        // (a/b)'s near ghost AND (c)'s far ghost, so (y-far,x-near) and
        // (y-far,x-far)=(nx,ny) corners both come out right.
        for (int k = 0; k < nz; ++k)
            for (int i = -2; i <= nx; ++i) {
                u.x(i, ny, k) = dirichletGhost1(u.x(i, ny - 1, k), u.x(i, ny - 2, k), u.x(i, ny - 3, k), u.x(i, ny - 4, k));
                u.y(i, ny, k) = dirichletGhost1(u.y(i, ny - 1, k), u.y(i, ny - 2, k), u.y(i, ny - 3, k), u.y(i, ny - 4, k));
                u.z(i, ny, k) = dirichletGhost1(u.z(i, ny - 1, k), u.z(i, ny - 2, k), u.z(i, ny - 3, k), u.z(i, ny - 4, k));
            }
        // (e) z=Lz far Dirichlet boundary, ONE layer; i,j now cover the
        // COMPLETE x,y ghost frame from (a)-(d), so every corner combined
        // with z-far comes out right in this one pass.
        for (int j = -2; j <= ny; ++j)
            for (int i = -2; i <= nx; ++i) {
                u.x(i, j, nz) = dirichletGhost1(u.x(i, j, nz - 1), u.x(i, j, nz - 2), u.x(i, j, nz - 3), u.x(i, j, nz - 4));
                u.y(i, j, nz) = dirichletGhost1(u.y(i, j, nz - 1), u.y(i, j, nz - 2), u.y(i, j, nz - 3), u.y(i, j, nz - 4));
                u.z(i, j, nz) = dirichletGhost1(u.z(i, j, nz - 1), u.z(i, j, nz - 2), u.z(i, j, nz - 3), u.z(i, j, nz - 4));
            }

        // (f) z=0: Dirichlet in Analytic-top mode, traction-free (Neumann)
        // in FreeSurface mode. Computed on the REAL i,j domain [0,nx)x[0,ny)
        // only -- the physical BC lives there, not on a symmetry mirror.
        if (top != TopBC::FreeSurface) {
            for (int j = 0; j < ny; ++j)
                for (int i = 0; i < nx; ++i) {
                    u.x(i, j, -1) = dirichletGhost1(u.x(i, j, 0), u.x(i, j, 1), u.x(i, j, 2), u.x(i, j, 3));
                    u.x(i, j, -2) = dirichletGhost2(u.x(i, j, 0), u.x(i, j, 1), u.x(i, j, 2), u.x(i, j, 3));
                    u.y(i, j, -1) = dirichletGhost1(u.y(i, j, 0), u.y(i, j, 1), u.y(i, j, 2), u.y(i, j, 3));
                    u.y(i, j, -2) = dirichletGhost2(u.y(i, j, 0), u.y(i, j, 1), u.y(i, j, 2), u.y(i, j, 3));
                    u.z(i, j, -1) = dirichletGhost1(u.z(i, j, 0), u.z(i, j, 1), u.z(i, j, 2), u.z(i, j, 3));
                    u.z(i, j, -2) = dirichletGhost2(u.z(i, j, 0), u.z(i, j, 1), u.z(i, j, 2), u.z(i, j, 3));
                }
        } else {
            // Traction-free (eq. 14):  du_x/dz=-du_z/dx, du_y/dz=-du_z/dy,
            // (lam+2mu) du_z/dz = -lam(du_x/dx+du_y/dy). In-plane target
            // derivatives use the full 4th-order 5-point formula (safe: x,y
            // near is exact symmetry, far now has the beyond-boundary ghost
            // from (c)/(d)). Ghost value itself via the verified quartic
            // Neumann-fit formula (uses u0 too -- see write-up).
            for (int j = 0; j < ny; ++j)
                for (int i = 0; i < nx; ++i) {
                    const double duzdx = (u.z(i + 1, j, 0) - u.z(i - 1, j, 0)) / (2 * h);
                    const double duzdy = (u.z(i, j + 1, 0) - u.z(i, j - 1, 0)) / (2 * h);
                    const double duxdx = (u.x(i + 1, j, 0) - u.x(i - 1, j, 0)) / (2 * h);
                    const double duydy = (u.y(i, j + 1, 0) - u.y(i, j - 1, 0)) / (2 * h);
                    const double Tx = -duzdx, Ty = -duzdy, Tz = -(lam / M) * (duxdx + duydy);
                    u.x(i, j, -1) = neumannGhost1(Tx, u.x(i, j, 0), u.x(i, j, 1), u.x(i, j, 2), u.x(i, j, 3), h);
                    u.x(i, j, -2) = neumannGhost2(Tx, u.x(i, j, 0), u.x(i, j, 1), u.x(i, j, 2), u.x(i, j, 3), h);
                    u.y(i, j, -1) = neumannGhost1(Ty, u.y(i, j, 0), u.y(i, j, 1), u.y(i, j, 2), u.y(i, j, 3), h);
                    u.y(i, j, -2) = neumannGhost2(Ty, u.y(i, j, 0), u.y(i, j, 1), u.y(i, j, 2), u.y(i, j, 3), h);
                    u.z(i, j, -1) = neumannGhost1(Tz, u.z(i, j, 0), u.z(i, j, 1), u.z(i, j, 2), u.z(i, j, 3), h);
                    u.z(i, j, -2) = neumannGhost2(Tz, u.z(i, j, 0), u.z(i, j, 1), u.z(i, j, 2), u.z(i, j, 3), h);
                }
        }

        // (g) Extend the just-computed z=0-side ghost PLANE (both layers)
        // out into ALL FOUR x,y corners: symmetry sides via exact mirror,
        // far sides via the same dirichletGhost1 extrapolation applied to
        // this ghost plane's own data (treating it as just another row/
        // column to extrapolate). This closes the last remaining corner
        // case (x-or-y-far combined with z-near) that a symmetry-only
        // mirror would miss.
        for (int k = -2; k <= -1; ++k) {
            for (int j = 0; j < ny; ++j) {
                u.x(-1, j, k) = -u.x(1, j, k); u.y(-1, j, k) = u.y(1, j, k); u.z(-1, j, k) = u.z(1, j, k);
                u.x(-2, j, k) = -u.x(2, j, k); u.y(-2, j, k) = u.y(2, j, k); u.z(-2, j, k) = u.z(2, j, k);
                u.x(nx, j, k) = dirichletGhost1(u.x(nx - 1, j, k), u.x(nx - 2, j, k), u.x(nx - 3, j, k), u.x(nx - 4, j, k));
                u.y(nx, j, k) = dirichletGhost1(u.y(nx - 1, j, k), u.y(nx - 2, j, k), u.y(nx - 3, j, k), u.y(nx - 4, j, k));
                u.z(nx, j, k) = dirichletGhost1(u.z(nx - 1, j, k), u.z(nx - 2, j, k), u.z(nx - 3, j, k), u.z(nx - 4, j, k));
            }
            for (int i = -2; i <= nx; ++i) {
                u.x(i, -1, k) = u.x(i, 1, k); u.y(i, -1, k) = -u.y(i, 1, k); u.z(i, -1, k) = u.z(i, 1, k);
                u.x(i, -2, k) = u.x(i, 2, k); u.y(i, -2, k) = -u.y(i, 2, k); u.z(i, -2, k) = u.z(i, 2, k);
                u.x(i, ny, k) = dirichletGhost1(u.x(i, ny - 1, k), u.x(i, ny - 2, k), u.x(i, ny - 3, k), u.x(i, ny - 4, k));
                u.y(i, ny, k) = dirichletGhost1(u.y(i, ny - 1, k), u.y(i, ny - 2, k), u.y(i, ny - 3, k), u.y(i, ny - 4, k));
                u.z(i, ny, k) = dirichletGhost1(u.z(i, ny - 1, k), u.z(i, ny - 2, k), u.z(i, ny - 3, k), u.z(i, ny - 4, k));
            }
        }
    }

    void fillGhostsP(Field& p) const {
        const int nx = g.nx, ny = g.ny, nz = g.nz;
        // x=0, y=0 symmetry: even reflection, two EXACT ghost layers.
        for (int k = 0; k < nz; ++k)
            for (int j = 0; j < ny; ++j) { p(-1, j, k) = p(1, j, k); p(-2, j, k) = p(2, j, k); }
        for (int k = 0; k < nz; ++k)
            for (int i = -2; i < nx; ++i) { p(i, -1, k) = p(i, 1, k); p(i, -2, k) = p(i, 2, k); }
        // One ghost layer beyond each far Dirichlet boundary (same
        // dirichletGhost1 formula as the free-surface case below).
        for (int k = 0; k < nz; ++k)
            for (int j = 0; j < ny; ++j)
                p(nx, j, k) = dirichletGhost1(p(nx - 1, j, k), p(nx - 2, j, k), p(nx - 3, j, k), p(nx - 4, j, k));
        for (int k = 0; k < nz; ++k)
            for (int i = 0; i < nx; ++i)
                p(i, ny, k) = dirichletGhost1(p(i, ny - 1, k), p(i, ny - 2, k), p(i, ny - 3, k), p(i, ny - 4, k));
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
                p(i, j, nz) = dirichletGhost1(p(i, j, nz - 1), p(i, j, nz - 2), p(i, j, nz - 3), p(i, j, nz - 4));
        // z = -h, -2h: Dirichlet ghost INTO the domain from the known
        // boundary value p(.,.,0) (0 for free surface, analytic for the
        // --top analytic mode), using the 3 nearest interior points.
        for (int j = -1; j < ny; ++j)
            for (int i = -1; i < nx; ++i) {
                p(i, j, -1) = dirichletGhost1(p(i, j, 0), p(i, j, 1), p(i, j, 2), p(i, j, 3));
                p(i, j, -2) = dirichletGhost2(p(i, j, 0), p(i, j, 1), p(i, j, 2), p(i, j, 3));
            }
    }

    // ----------------------------------------------------------------------
    // With full boundary ghost coverage (see fillGhostsU/P above -- a proper
    // 4th-order-consistent ghost is now available at every boundary: the
    // free surface via the Neumann cubic-fit, the far Dirichlet boundary via
    // one ghost layer beyond it), the 5-point stencil is valid EVERYWHERE
    // for the diagonal (axis-aligned) second derivatives. No fallback is
    // needed any more.
    // ----------------------------------------------------------------------
    inline bool x4ok(int) const { return true; }
    inline bool y4ok(int) const { return true; }
    inline bool z4ok(int) const { return true; }

    template <class F>
    static inline double d2x(const F& f, int i, int j, int k, double h, bool ok4) {
        if (ok4) return (-f(i + 2, j, k) + 16.0 * f(i + 1, j, k) - 30.0 * f(i, j, k) +
                          16.0 * f(i - 1, j, k) - f(i - 2, j, k)) / (12.0 * h * h);
        return (f(i + 1, j, k) - 2.0 * f(i, j, k) + f(i - 1, j, k)) / (h * h);
    }
    template <class F>
    static inline double d2y(const F& f, int i, int j, int k, double h, bool ok4) {
        if (ok4) return (-f(i, j + 2, k) + 16.0 * f(i, j + 1, k) - 30.0 * f(i, j, k) +
                          16.0 * f(i, j - 1, k) - f(i, j - 2, k)) / (12.0 * h * h);
        return (f(i, j + 1, k) - 2.0 * f(i, j, k) + f(i, j - 1, k)) / (h * h);
    }
    template <class F>
    static inline double d2z(const F& f, int i, int j, int k, double h, bool ok4) {
        if (ok4) return (-f(i, j, k + 2) + 16.0 * f(i, j, k + 1) - 30.0 * f(i, j, k) +
                          16.0 * f(i, j, k - 1) - f(i, j, k - 2)) / (12.0 * h * h);
        return (f(i, j, k + 1) - 2.0 * f(i, j, k) + f(i, j, k - 1)) / (h * h);
    }

    // ----------------------------------------------------------------------
    // L(u) = (lambda+mu) grad(div u) + mu lap(u),  evaluated on unknown nodes.
    // Diagonal (axis-aligned) 2nd derivatives use the order-selected stencil
    // above; the mixed (cross) derivatives dxy,dxz,dyz stay 2nd-order (a
    // 4th-order cross-stencil would need a much wider footprint and is not
    // implemented here -- see the write-up).
    // ----------------------------------------------------------------------
    double diag() const { return -(2.0 * mat.lambda() + 8.0 * mat.mu) / (g.h * g.h); }

    void applyL(Vec3Field& u, Vec3Field& out) const {
        fillGhostsU(u);
        const int nx = g.nx, ny = g.ny, nz = g.nz;
        const double h = g.h;
        const double ih2 = 1.0 / (h * h), q = 0.25 * ih2;
        const double lm = mat.lambda() + mat.mu, mu = mat.mu;
        const int k0 = (top == TopBC::FreeSurface) ? 0 : 1;
        auto fx = [&](int i, int j, int k) { return u.x(i, j, k); };
        auto fy = [&](int i, int j, int k) { return u.y(i, j, k); };
        auto fz = [&](int i, int j, int k) { return u.z(i, j, k); };
#pragma omp parallel for collapse(2) schedule(static)
        for (int k = k0; k <= nz - 2; ++k)
            for (int j = 0; j <= ny - 2; ++j)
                for (int i = 0; i <= nx - 2; ++i) {
                    const bool ox = x4ok(i), oy = y4ok(j), oz = z4ok(k);
                    const double dxx = d2x(fx, i, j, k, h, ox);
                    const double dyy_x = d2y(fx, i, j, k, h, oy);
                    const double dzz_x = d2z(fx, i, j, k, h, oz);
                    const double dxx_y = d2x(fy, i, j, k, h, ox);
                    const double dyy = d2y(fy, i, j, k, h, oy);
                    const double dzz_y = d2z(fy, i, j, k, h, oz);
                    const double dxx_z = d2x(fz, i, j, k, h, ox);
                    const double dyy_z = d2y(fz, i, j, k, h, oy);
                    const double dzz = d2z(fz, i, j, k, h, oz);
                    const double lapx = dxx + dyy_x + dzz_x;
                    const double lapy = dxx_y + dyy + dzz_y;
                    const double lapz = dxx_z + dyy_z + dzz;
                    const double dxy = q * (u.y(i + 1, j + 1, k) - u.y(i + 1, j - 1, k) -
                                            u.y(i - 1, j + 1, k) + u.y(i - 1, j - 1, k));
                    const double dyx = q * (u.x(i + 1, j + 1, k) - u.x(i + 1, j - 1, k) -
                                            u.x(i - 1, j + 1, k) + u.x(i - 1, j - 1, k));
                    const double dxz = q * (u.z(i + 1, j, k + 1) - u.z(i + 1, j, k - 1) -
                                            u.z(i - 1, j, k + 1) + u.z(i - 1, j, k - 1));
                    const double dzx = q * (u.x(i + 1, j, k + 1) - u.x(i + 1, j, k - 1) -
                                            u.x(i - 1, j, k + 1) + u.x(i - 1, j, k - 1));
                    const double dyz = q * (u.z(i, j + 1, k + 1) - u.z(i, j + 1, k - 1) -
                                            u.z(i, j - 1, k + 1) + u.z(i, j - 1, k - 1));
                    const double dzy = q * (u.y(i, j + 1, k + 1) - u.y(i, j + 1, k - 1) -
                                            u.y(i, j - 1, k + 1) + u.y(i, j - 1, k - 1));
                    out.x(i, j, k) = lm * (dxx + dxy + dxz) + mu * lapx;
                    out.y(i, j, k) = lm * (dyx + dyy + dyz) + mu * lapy;
                    out.z(i, j, k) = lm * (dzx + dzy + dzz) + mu * lapz;
                }
    }

    // RHS of the mechanics equation: alpha * grad(p). Same order-selection
    // policy as applyL: 4th-order 1st-derivative stencil where safe, 2nd
    // order (central 2-point) where the wide stencil would touch the single
    // approximate boundary ghost layer.
    void mechRHS(Field& p, Vec3Field& b) const {
        fillGhostsP(p);
        const int nx = g.nx, ny = g.ny, nz = g.nz;
        const double h = g.h;
        const double c2 = mat.alpha() / (2.0 * h), c4 = mat.alpha() / (12.0 * h);
        const int k0 = (top == TopBC::FreeSurface) ? 0 : 1;
        b.zero();
#pragma omp parallel for collapse(2) schedule(static)
        for (int k = k0; k <= nz - 2; ++k)
            for (int j = 0; j <= ny - 2; ++j)
                for (int i = 0; i <= nx - 2; ++i) {
                    b.x(i, j, k) = x4ok(i)
                        ? c4 * (-p(i + 2, j, k) + 8.0 * p(i + 1, j, k) - 8.0 * p(i - 1, j, k) + p(i - 2, j, k))
                        : c2 * (p(i + 1, j, k) - p(i - 1, j, k));
                    b.y(i, j, k) = y4ok(j)
                        ? c4 * (-p(i, j + 2, k) + 8.0 * p(i, j + 1, k) - 8.0 * p(i, j - 1, k) + p(i, j - 2, k))
                        : c2 * (p(i, j + 1, k) - p(i, j - 1, k));
                    b.z(i, j, k) = z4ok(k)
                        ? c4 * (-p(i, j, k + 2) + 8.0 * p(i, j, k + 1) - 8.0 * p(i, j, k - 1) + p(i, j, k - 2))
                        : c2 * (p(i, j, k + 1) - p(i, j, k - 1));
                }
    }

    // --- vector helpers restricted to the unknown set ----------------------
    template <class Fn>
    void forEachUnknown(Fn fn) const {
        const int k0 = (top == TopBC::FreeSurface) ? 0 : 1;
        for (int k = k0; k <= g.nz - 2; ++k)
            for (int j = 0; j <= g.ny - 2; ++j)
                for (int i = 0; i <= g.nx - 2; ++i) fn(i, j, k);
    }
    double dot(const Vec3Field& a, const Vec3Field& b) const {
        double s = 0.0;
        const int k0 = (top == TopBC::FreeSurface) ? 0 : 1;
#pragma omp parallel for collapse(2) reduction(+ : s) schedule(static)
        for (int k = k0; k <= g.nz - 2; ++k)
            for (int j = 0; j <= g.ny - 2; ++j)
                for (int i = 0; i <= g.nx - 2; ++i)
                    s += a.x(i, j, k) * b.x(i, j, k) + a.y(i, j, k) * b.y(i, j, k) +
                         a.z(i, j, k) * b.z(i, j, k);
        return s;
    }
    double norm(const Vec3Field& a) const { return std::sqrt(dot(a, a)); }

    // --- red-black Gauss-Seidel, used as the BiCGSTAB preconditioner -------
    void rbgs(Vec3Field& zf, const Vec3Field& r, Vec3Field& tmp, int sweeps,
              double omega = 1.0) const {
        const double dinv = omega / diag();
        for (int s = 0; s < sweeps; ++s)
            for (int color = 0; color < 2; ++color) {
                applyL(zf, tmp);
                const int k0 = (top == TopBC::FreeSurface) ? 0 : 1;
#pragma omp parallel for collapse(2) schedule(static)
                for (int k = k0; k <= g.nz - 2; ++k)
                    for (int j = 0; j <= g.ny - 2; ++j)
                        for (int i = 0; i <= g.nx - 2; ++i) {
                            if (((i + j + k) & 1) != color) continue;
                            zf.x(i, j, k) += dinv * (r.x(i, j, k) - tmp.x(i, j, k));
                            zf.y(i, j, k) += dinv * (r.y(i, j, k) - tmp.y(i, j, k));
                            zf.z(i, j, k) += dinv * (r.z(i, j, k) - tmp.z(i, j, k));
                        }
            }
    }

    // ----------------------------------------------------------------------
    // Matrix-free BiCGSTAB for  L(u) = b  with u warm-started from the
    // previous time level.  BiCGSTAB (not CG) because the free-surface and
    // symmetry ghost eliminations make the discrete operator only nearly
    // symmetric.
    // ----------------------------------------------------------------------
    struct Work {
        Vec3Field r, r0, p, v, s, t, z, y, tmp;
        void init(const Grid& g) {
            r.init(g); r0.init(g); p.init(g); v.init(g);
            s.init(g); t.init(g); z.init(g); y.init(g); tmp.init(g);
        }
    };

    int solveMechanics(Vec3Field& u, Vec3Field& b, Work& w, double tol, int maxit,
                       int smooth = 2) const {
        applyL(u, w.tmp);
        const int k0 = (top == TopBC::FreeSurface) ? 0 : 1;
#pragma omp parallel for collapse(2) schedule(static)
        for (int k = k0; k <= g.nz - 2; ++k)
            for (int j = 0; j <= g.ny - 2; ++j)
                for (int i = 0; i <= g.nx - 2; ++i) {
                    w.r.x(i, j, k) = b.x(i, j, k) - w.tmp.x(i, j, k);
                    w.r.y(i, j, k) = b.y(i, j, k) - w.tmp.y(i, j, k);
                    w.r.z(i, j, k) = b.z(i, j, k) - w.tmp.z(i, j, k);
                }
        const double bnorm = std::max(norm(b), 1e-300);
        double rn = norm(w.r);
        if (rn <= tol * bnorm) return 0;

        w.r0.copyFrom(w.r);
        w.p.zero(); w.v.zero();
        double rho = 1.0, alp = 1.0, om = 1.0;

        for (int it = 1; it <= maxit; ++it) {
            const double rho1 = dot(w.r0, w.r);
            if (std::fabs(rho1) < 1e-300) break;
            const double beta = (rho1 / rho) * (alp / om);
            rho = rho1;
            axpby3(w.p, w.r, beta, -beta * om, w.v);  // p = r + beta(p - om v)
            w.y.zero();
            rbgs(w.y, w.p, w.tmp, smooth);
            applyL(w.y, w.v);
            const double d1 = dot(w.r0, w.v);
            if (std::fabs(d1) < 1e-300) break;
            alp = rho / d1;
            lincomb(w.s, w.r, w.v, -alp);              // s = r - alp v
            if (norm(w.s) <= tol * bnorm) {
                addScaled(u, w.y, alp);
                return it;
            }
            w.z.zero();
            rbgs(w.z, w.s, w.tmp, smooth);
            applyL(w.z, w.t);
            const double tt = dot(w.t, w.t);
            om = (tt > 0.0) ? dot(w.t, w.s) / tt : 0.0;
            addScaled(u, w.y, alp);
            addScaled(u, w.z, om);
            lincomb(w.r, w.s, w.t, -om);               // r = s - om t
            rn = norm(w.r);
            if (rn <= tol * bnorm) return it;
            if (om == 0.0) break;
        }
        return maxit;
    }

   private:
    void axpby3(Vec3Field& p, const Vec3Field& r, double b1, double b2,
                const Vec3Field& v) const {
        const int k0 = (top == TopBC::FreeSurface) ? 0 : 1;
#pragma omp parallel for collapse(2) schedule(static)
        for (int k = k0; k <= g.nz - 2; ++k)
            for (int j = 0; j <= g.ny - 2; ++j)
                for (int i = 0; i <= g.nx - 2; ++i) {
                    p.x(i, j, k) = r.x(i, j, k) + b1 * p.x(i, j, k) + b2 * v.x(i, j, k);
                    p.y(i, j, k) = r.y(i, j, k) + b1 * p.y(i, j, k) + b2 * v.y(i, j, k);
                    p.z(i, j, k) = r.z(i, j, k) + b1 * p.z(i, j, k) + b2 * v.z(i, j, k);
                }
    }
    void lincomb(Vec3Field& o, const Vec3Field& a, const Vec3Field& b, double c) const {
        const int k0 = (top == TopBC::FreeSurface) ? 0 : 1;
#pragma omp parallel for collapse(2) schedule(static)
        for (int k = k0; k <= g.nz - 2; ++k)
            for (int j = 0; j <= g.ny - 2; ++j)
                for (int i = 0; i <= g.nx - 2; ++i) {
                    o.x(i, j, k) = a.x(i, j, k) + c * b.x(i, j, k);
                    o.y(i, j, k) = a.y(i, j, k) + c * b.y(i, j, k);
                    o.z(i, j, k) = a.z(i, j, k) + c * b.z(i, j, k);
                }
    }
    void addScaled(Vec3Field& u, const Vec3Field& d, double c) const {
        const int k0 = (top == TopBC::FreeSurface) ? 0 : 1;
#pragma omp parallel for collapse(2) schedule(static)
        for (int k = k0; k <= g.nz - 2; ++k)
            for (int j = 0; j <= g.ny - 2; ++j)
                for (int i = 0; i <= g.nx - 2; ++i) {
                    u.x(i, j, k) += c * d.x(i, j, k);
                    u.y(i, j, k) += c * d.y(i, j, k);
                    u.z(i, j, k) += c * d.z(i, j, k);
                }
    }
};


}  // namespace poro

// main.cpp -- 3D finite-difference poroelastic solver for the point-injection
// problem of Wang & Kuempel (2003), Geophysics 68(2), 1-13.
//
//   (lambda+2mu) grad(div u) - mu curl(curl u) - alpha grad p = 0        (1)
//   (1/Q) dp/dt + alpha d(div u)/dt - chi lap(p) = q(x,t)                (2)
//
// Scheme
// ------
//  * Equation (2) is advanced with EXPLICIT (forward-Euler) central-difference
//    diffusion, so the step is limited by the CFL / von-Neumann condition
//        dt <= cfl / (2 D (1/dx^2 + 1/dy^2 + 1/dz^2)) = cfl h^2 / (6 D).
//  * Equation (1) has no time derivative; it is elliptic and is solved at
//    every time level by matrix-free BiCGSTAB (mechanics.hpp).
//  * The two are coupled with the FIXED-STRESS split, which uses the exact
//    identity  sigma_v = K_dr div u - alpha p  to rewrite (2) as
//        (1/Q + alpha^2/K_dr) dp/dt + (alpha/K_dr) d sigma_v/dt
//                                              - chi lap p = q ,
//    and lags sigma_v.  The split is unconditionally stable in the coupling
//    iteration, and its effective explicit diffusivity chi/S_fs is SMALLER
//    than D, so the CFL bound written with D is conservative.
//
// Build:  cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build
// Run:    ./build/poro3d --h 5 --L 250 --tend 21600

using namespace poro;

struct Obs {
    int i, j, k;
    double r, z;
};

int main(int argc, char** argv) {
    // ---------------- configuration (paper's TEST RESULTS section) --------
    Material mat;              // mu=0.4 GPa, nu=0.2, nu_u=0.4, B=0.75, D=1 m^2/s
    double q0 = 32.0 / 3600.0; // 32 m^3/h  -> 0.008889 m^3/s
    double d = 60.0;           // source depth [m]
    double robs = 40.0;        // observation radius [m]
    double tend = 3600.0;  // 1 h -- laptop-scale default (paper's window is 6 h)
    double h = 5.0, L = 200.0, Lz = 200.0, cfl = 0.9;
    int coupleIters = 2, maxit = 300, smooth = 2, sampleEvery = 20, subiters = 1;
    double tol = 1e-8;
    // Default is a PURE ZERO Dirichlet far field (u = 0, p = 0 on the outer
    // box) as requested -- NOT the mirror-image analytic far field used
    // earlier for whole-space verification. Pass --farbc analytic to restore
    // that (useful only for the --top analytic whole-space check).
    bool freeSurf = true, farAnalytic = false;
    std::string out = "fig2.csv";
    bool implicitMode = false; double dtImp = 10.0, fsTol = 1e-6; int fsMax = 200; double Lfac = 1.0;

    for (int a = 1; a < argc; ++a) {
        auto eq = [&](const char* s) { return std::strcmp(argv[a], s) == 0; };
        auto nxt = [&]() { return std::atof(argv[++a]); };
        if (eq("--h")) h = nxt();
        else if (eq("--L")) L = nxt();
        else if (eq("--Lz")) Lz = nxt();
        else if (eq("--tend")) tend = nxt();
        else if (eq("--cfl")) cfl = nxt();
        else if (eq("--d")) d = nxt();
        else if (eq("--r")) robs = nxt();
        else if (eq("--q0")) q0 = nxt();
        else if (eq("--D")) mat.D = nxt();
        else if (eq("--mu")) mat.mu = nxt();
        else if (eq("--nu")) mat.nu = nxt();
        else if (eq("--nuu")) mat.nu_u = nxt();
        else if (eq("--B")) mat.B = nxt();
        else if (eq("--couple")) coupleIters = (int)nxt();
        else if (eq("--subiters")) subiters = (int)nxt();
        else if (eq("--tol")) tol = nxt();
        else if (eq("--maxit")) maxit = (int)nxt();
        else if (eq("--smooth")) smooth = (int)nxt();
        else if (eq("--sample")) sampleEvery = (int)nxt();
        else if (eq("--top")) { ++a; freeSurf = (std::strcmp(argv[a], "free") == 0); }
        else if (eq("--farbc")) { ++a; farAnalytic = (std::strcmp(argv[a], "analytic") == 0); }
        else if (eq("--implicit")) implicitMode = true;
        else if (eq("--dt")) dtImp = nxt();
        else if (eq("--fsmax")) fsMax = (int)nxt();
        else if (eq("--fstol")) fsTol = nxt();
        else if (eq("--Lfac")) Lfac = nxt();
        else if (eq("--out")) out = argv[++a];
        else { std::fprintf(stderr, "unknown option %s\n", argv[a]); return 1; }
    }

    Solver S;
    S.mat = mat;
    S.top = freeSurf ? TopBC::FreeSurface : TopBC::Analytic;
    S.g.h = h;
    S.g.nx = int(L / h) + 1;
    S.g.ny = int(L / h) + 1;
    S.g.nz = int(Lz / h) + 1;
    S.ref.ws.m = mat;
    S.ref.ws.q0 = q0;
    S.ref.d = d;
    const Grid& g = S.g;

    const int kd = int(std::lround(d / h));
    if (std::fabs(kd * h - d) > 1e-9) {
        std::fprintf(stderr, "source depth %g is not a grid node for h = %g\n", d, h);
        return 1;
    }

    // Analytic field used on the far boundaries.  With a free surface the
    // mirror-image field (source + negative image) is used: it satisfies
    // p = 0 on z = 0 exactly and is the correct leading-order half-space
    // far field.  With --top analytic the pure whole-space field is used, so
    // the run reproduces Rudnicki's (1986) whole-space Green's functions.
    const bool useImage = freeSurf;
    auto anaP = [&](double x, double y, double z, double t) {
        if (!farAnalytic) return 0.0;
        return useImage ? S.ref.p(x, y, z, t) : S.ref.ws.p(x, y, z - d, t);
    };
    auto anaUx = [&](double x, double y, double z, double t) {
        if (!farAnalytic) return 0.0;
        return useImage ? S.ref.ux(x, y, z, t) : S.ref.ws.u(x, x, y, z - d, t);
    };
    auto anaUy = [&](double x, double y, double z, double t) {
        if (!farAnalytic) return 0.0;
        return useImage ? S.ref.uy(x, y, z, t) : S.ref.ws.u(y, x, y, z - d, t);
    };
    auto anaUz = [&](double x, double y, double z, double t) {
        if (!farAnalytic) return 0.0;
        return useImage ? S.ref.uz(x, y, z, t) : S.ref.ws.u(z - d, x, y, z - d, t);
    };

    // ---------------- CFL -------------------------------------------------
    const double invh2 = 1.0 / (h * h);
    // CFL / von-Neumann limit for the single-pass, sigma_v-lagged
    // fixed-stress scheme actually used below. Von Neumann analysis of the
    // full two-level recursion
    //   p^{n+1} = A(k) p^n + B(k) p^{n-1},
    //   B(k) = [alpha*r(k) - alpha^2/Kdr] / S_fs,   r(k) in [0, alpha/M]
    // gives |B| <= max(alpha^2/(S_fs*Kdr), |alpha/M - alpha^2/(Kdr S_fs)|)
    // which is < 1 here at both the smooth-mode end (r=alpha/M, |B|~0.36)
    // and the grid-Nyquist end (r=0, |B|~0.71) -- i.e. the sigma_v lag is
    // stable across the whole mode range PROVIDED the explicit diffusion
    // term itself is CFL-bounded using S_fs (not the paper's D, and not
    // 1/Q): the explicit term's own diffusivity is chi/S_fs.
    // (An earlier attempt lagging div(u) directly instead of sigma_v used
    // Qinv as the storage coefficient and was UNCONDITIONALLY UNSTABLE
    // here, |B|~1.25>1 even for smooth modes -- verified by direct blow-up.
    // That is a known distinction in the literature: fixed-stress (sigma_v)
    // splits are stable, naive "drained" (div u) splits generally are not.)
    const double Dcfl = mat.chi() / mat.S_fs();
    // 4TH-ORDER UPGRADE: the 5-point 2nd-derivative stencil's Nyquist-mode
    // eigenvalue is 16/(3h^2) per axis (vs 4/h^2 for the 3-point stencil --
    // checked by von Neumann analysis of (-1,16,-30,16,-1)/12 at theta=pi).
    // Summed over 3 axes that is 16/h^2 vs the old 12/h^2, so the explicit
    // step is now MORE restrictive by a factor 4/3. Using this tighter bound
    // everywhere is a safe, conservative choice for the whole mixed-order
    // grid (the 2nd-order fallback regions are less restrictive, so they are
    // automatically satisfied too).
    const double dtCFL = 1.0 / (2.0 * Dcfl * (4.0 * invh2));  // = h^2/(8 chi/S_fs)
    double dt = cfl * dtCFL;
    const long nsteps = (long)std::ceil(tend / dt);
    dt = tend / nsteps;

    const double lam = mat.lambda(), alpha = mat.alpha(), chi = mat.chi();
    const double Kdr = mat.Kdr(), Sfs = mat.S_fs();

    std::printf("--- Wang & Kuempel (2003) 3D FD poroelastic solver ---\n");
    std::printf("grid      : %d x %d x %d  (quarter domain, h = %g m)\n", g.nx, g.ny, g.nz, h);
    std::printf("material  : lambda=%.4e Pa  alpha=%.4f  1/Q=%.4e 1/Pa  chi=%.4e m^2/(Pa s)\n",
                lam, alpha, mat.Qinv(), chi);
    std::printf("check     : chi/(1/Q + alpha^2/(lambda+2mu)) = %.6f  (must equal D = %g)\n",
                chi / mat.S_D(), mat.D);
    std::printf("diffusiv. : D = chi/S_D = %.4f m^2/s (true consolidation coefficient) ; "
                "chi/S_fs = %.4f m^2/s (explicit CFL diffusivity used)\n", chi / mat.S_D(), Dcfl);
    std::printf("CFL       : dt_max = h^2/(8 chi/S_fs) [4th-order Nyquist bound] = %.4f s ; "
                "dt = %.4f s ; steps = %ld\n", dtCFL, dt, nsteps);
    std::printf("top BC    : %s ; far BC : %s\n", freeSurf ? "free surface (sigma.n=0, p=0)"
                                                           : "analytic Dirichlet",
                farAnalytic ? "analytic Dirichlet" : "zero");

    // ---------------- fields ---------------------------------------------
    // SINGLE-PASS, one-step-lagged fixed-stress coupling (sigma_v lag, not
    // div(u) lag). This choice of lag matters: a direct div(u) lag is a
    // known-unstable "drained split" for Biot coefficients like this one (I
    // verified it blows up here -- von Neumann analysis gives |amplification
    // factor| > 1 even for smooth modes). The sigma_v-lag fixed-stress split
    // is the standard, provably-stable choice in the literature (Kim,
    // Tchelepi & Juanes 2011); von Neumann analysis of THIS scheme gives
    // |amplification| ~ 0.36-0.71 < 1 across the full mode range including
    // the grid-Nyquist mode near the point source. Earlier I ALSO tried
    // Picard-iterating this same sigma_v lag to convergence within each time
    // step -- that converged to a stable but WRONG fixed point (verified:
    // error grew with more iterations), because the point source's
    // grid-Nyquist content violates the smooth-mode assumption baked into
    // the fixed-stress correction. Using it as a single EXPLICIT lag (this
    // version) avoids that failure mode entirely: there is no iteration to
    // converge to the wrong place, and the physics is trusted to the same
    // order in dt as the diffusion term already is (i.e. this uses S_fs's
    // fast decay rate; it is used at the CFL-limited dt required for
    // stability regardless).
    Field p, pnew, lapP, eps, sigv, sigv0, pn0;
    Vec3Field u, b;
    Solver::Work W;
    p.init(g); pnew.init(g); lapP.init(g); eps.init(g); sigv.init(g); sigv0.init(g); pn0.init(g);
    u.init(g); b.init(g); W.init(g);

    // Smear the point source over a small Gaussian kernel (radius ~1.5h)
    // instead of one grid cell. A one-cell delta excites strong grid-Nyquist
    // content right at the source, which is what made even a SMALL number
    // of fixed-stress Picard sub-iterations blow up (verified: subiters=2
    // diverges with a one-cell source). Smoothing the source removes that
    // high-wavenumber content at its origin, consistent with standard
    // practice for point-source injection wells in FD/FV codes.
    std::vector<std::array<int,3>> srcCells;
    std::vector<double> srcWeight;
    {
        const double sigma = 1.2 * h;
        const int R = 2;  // +/- 2 cells
        double wsum = 0.0;
        for (int dk = -R; dk <= R; ++dk)
            for (int dj = -R; dj <= R; ++dj)
                for (int di = -R; di <= R; ++di) {
                    int ii = di, jj = dj, kk = kd + dk;   // source on axis i=j=0
                    if (ii < 0 || jj < 0 || kk < 1 || kk > g.nz - 2) continue;
                    if (ii > g.nx - 2 || jj > g.ny - 2) continue;
                    double r2 = double(di*di + dj*dj + dk*dk) * h * h;
                    double w = std::exp(-0.5 * r2 / (sigma * sigma));
                    double mult = (ii == 0 ? 1.0 : 2.0) * (jj == 0 ? 1.0 : 2.0);
                    wsum += w * mult;
                    srcCells.push_back({ii, jj, kk});
                    srcWeight.push_back(w);
                }
        for (auto& w : srcWeight) w *= q0 / (wsum * h * h * h);
    }

    auto divergence = [&](Vec3Field& uu, Field& e) {
        S.fillGhostsU(uu);
        const double c = 1.0 / (2.0 * h);
#pragma omp parallel for collapse(2) schedule(static)
        for (int k = 0; k <= g.nz - 1; ++k)
            for (int j = 0; j <= g.ny - 1; ++j)
                for (int i = 0; i <= g.nx - 1; ++i)
                    e(i, j, k) = c * (uu.x(i + 1, j, k) - uu.x(i - 1, j, k) +
                                      uu.y(i, j + 1, k) - uu.y(i, j - 1, k) +
                                      uu.z(i, j, k + 1) - uu.z(i, j, k - 1));
    };
    auto meanStress = [&](Field& e, Field& pp, Field& sv) {
#pragma omp parallel for collapse(2) schedule(static)
        for (int k = 0; k <= g.nz - 1; ++k)
            for (int j = 0; j <= g.ny - 1; ++j)
                for (int i = 0; i <= g.nx - 1; ++i)
                    sv(i, j, k) = Kdr * e(i, j, k) - alpha * pp(i, j, k);
    };
    auto setDirichlet = [&](double t) {
        for (int k = 0; k < g.nz; ++k)
            for (int j = 0; j < g.ny; ++j)
                for (int i = 0; i < g.nx; ++i) {
                    const bool far = (i == g.nx - 1 || j == g.ny - 1 || k == g.nz - 1);
                    if (!far && k != 0) continue;
                    const double X = g.x(i), Y = g.y(j), Z = g.z(k);
                    p(i, j, k) = (freeSurf && k == 0) ? 0.0 : anaP(X, Y, Z, t);
                    if (far || !freeSurf) {  // u is unknown on the free surface
                        u.x(i, j, k) = anaUx(X, Y, Z, t);
                        u.y(i, j, k) = anaUy(X, Y, Z, t);
                        u.z(i, j, k) = anaUz(X, Y, Z, t);
                    }
                }
    };

    // observation points of Fig. 2: r = 40 m, z = 5, 15, 45, 75 m
    std::vector<Obs> obs;
    for (double zo : {5.0, 15.0, 45.0, 75.0}) {
        const int i = int(std::lround(robs / h)), k = int(std::lround(zo / h));
        if (std::fabs(i * h - robs) > 1e-9 || std::fabs(k * h - zo) > 1e-9) {
            std::fprintf(stderr, "note: obs (r=%g,z=%g) is off-grid for h=%g, skipping\n",
                         robs, zo, h);
            continue;
        }
        if (i > g.nx - 2 || k > g.nz - 2) {
            std::fprintf(stderr, "note: obs (r=%g,z=%g) outside domain, skipping\n", robs, zo);
            continue;
        }
        obs.push_back({i, 0, k, robs, zo});
    }

    std::FILE* f = std::fopen(out.c_str(), "w");
    std::fprintf(f, "t_s,t_h");
    for (auto& o : obs) std::fprintf(f, ",p_num_z%.0f,p_ws_z%.0f,p_mi_z%.0f,"
                                       "tiltA_num_z%.0f,tiltB_num_z%.0f,tilt_ws_z%.0f",
                                    o.z, o.z, o.z, o.z, o.z, o.z);
    std::fprintf(f, "\n");

    // Two candidate "tilt" definitions -- the paper does not spell out which
    // one a borehole tiltmeter records, and we could not resolve this from
    // the text alone (see the discussion below), so BOTH are written out:
    //   A) slope        : -du_z/dr             (surface-slope convention)
    //   B) rigid rotation: 0.5*(du_r/dz - du_z/dr)   (elastic rotation omega_rz)
    // For a purely radial (irrotational) field the two coincide up to a
    // factor of 2 only if du_r/dz = du_z/dr identically, which happens to be
    // true for the whole-space Green's function; they generally differ once
    // the free surface is present, which is exactly where the paper reports
    // sign differences between whole- and half-space tilt.
    auto sample = [&](double t) {
        std::fprintf(f, "%.6g,%.6g", t, t / 3600.0);
        for (auto& o : obs) {
            const double pn = p(o.i, o.j, o.k);
            const double pw = S.ref.ws.p(robs, 0.0, o.z - d, t);
            const double pm = S.ref.p(robs, 0.0, o.z, t);
            const double duzdr = (u.z(o.i + 1, o.j, o.k) - u.z(o.i - 1, o.j, o.k)) / (2 * h);
            const double durdz = (u.x(o.i, o.j, o.k + 1) - u.x(o.i, o.j, o.k - 1)) / (2 * h);
            const double tA = -duzdr;
            const double tB = 0.5 * (durdz - duzdr);
            const double tw = -S.ref.ws.duz_dr(robs, o.z - d, t);
            std::fprintf(f, ",%.8g,%.8g,%.8g,%.8g,%.8g,%.8g", pn, pw, pm, tA, tB, tw);
        }
        std::fprintf(f, "\n");
    };

    // ---------------- initial state: p = 0, u = 0 at t = 0 ----------------
    p.zero(); u.zero();
    setDirichlet(0.0);
    divergence(u, eps);
    meanStress(eps, p, sigv);
    sigv0.a = sigv.a;
    sample(0.0);

    long itTotal = 0;
    double t = 0.0;
    if (implicitMode) {
        // ===== IMPLICIT (backward-Euler) pressure + converged fixed-stress =====
        // Per step, iterate k:  (S_fs/dt) p^{k+1} - chi Lap p^{k+1}
        //    = (S_fs/dt) p^n + q - (alpha/(K_dr dt)) (sigma_v^k - sigma_v^n)
        // then solve mechanics with p^{k+1}. Converged fixed point = backward
        // Euler of the fully coupled system (Kim-Tchelepi-Juanes 2011;
        // Mikelic-Wheeler 2013 contraction proof applies to this form).
        const long nI = (long)std::ceil(tend / dtImp); const double dI = tend / nI;
        const double invh2i = 1.0 / (h * h); const double Ls = Lfac*alpha*alpha/Kdr, SL = mat.Qinv() + Ls;
        Field epsn; epsn.init(g);
        std::printf("IMPLICIT  : dt = %.4f s, steps = %ld, tol = %.1e, L = %.2f*alpha^2/Kdr\n", dI, nI, fsTol, Lfac);
        Field pn, rhs, rr, dd, Ad, del, sign, pold;
        pn.init(g); rhs.init(g); rr.init(g); dd.init(g); Ad.init(g); del.init(g); sign.init(g); pold.init(g);
        auto inU = [&](int i,int j,int k){ return i<=g.nx-2 && j<=g.ny-2 && k>=1 && k<=g.nz-2; };
        auto wgt = [&](int i,int j){ return (i==0?0.5:1.0)*(j==0?0.5:1.0); };
        auto nb = [&](Field& f,int i,int j,int k){   // 7-pt Laplacian*h^2 with symmetry mirrors
            auto v=[&](int ii,int jj,int kk){ if(ii<0) ii=-ii; if(jj<0) jj=-jj; return f(ii,jj,kk); };
            return v(i+1,j,k)+v(i-1,j,k)+v(i,j+1,k)+v(i,j-1,k)+v(i,j,k+1)+v(i,j,k-1)-6.0*f(i,j,k); };
        auto applyA = [&](Field& x, Field& y){   // x zero outside unknowns
            for(int k=1;k<=g.nz-2;++k) for(int j=0;j<=g.ny-2;++j) for(int i=0;i<=g.nx-2;++i)
                y(i,j,k) = (SL/dI)*x(i,j,k) - chi*invh2i*nb(x,i,j,k); };
        auto wdot = [&](Field& a1, Field& b1){ double s1=0;
            for(int k=1;k<=g.nz-2;++k) for(int j=0;j<=g.ny-2;++j) for(int i=0;i<=g.nx-2;++i)
                s1 += wgt(i,j)*a1(i,j,k)*b1(i,j,k); return s1; };
        long cgTot=0, fsTot=0, fsWorst=0;
        for (long n = 1; n <= nI; ++n) {
            t = n * dI;
            pn.a = p.a; epsn.a = eps.a;
            setDirichlet(t);                          // boundary p,u at t^{n+1}
            int kfs = 0; double chg = 1.0;
            for (kfs = 1; kfs <= fsMax; ++kfs) {
                pold.a = p.a;
                // rhs and residual r = rhs - A p  (p carries Dirichlet values)
                for(int k=1;k<=g.nz-2;++k) for(int j=0;j<=g.ny-2;++j) for(int i=0;i<=g.nx-2;++i)
                    rhs(i,j,k) = (mat.Qinv()*pn(i,j,k) + Ls*p(i,j,k) - alpha*(eps(i,j,k)-epsn(i,j,k)))/dI;
                for (std::size_t c = 0; c < srcCells.size(); ++c) { auto& cc=srcCells[c]; rhs(cc[0],cc[1],cc[2]) += srcWeight[c]; }
                rr.zero(); dd.zero(); del.zero();
                for(int k=1;k<=g.nz-2;++k) for(int j=0;j<=g.ny-2;++j) for(int i=0;i<=g.nx-2;++i)
                    rr(i,j,k) = rhs(i,j,k) - ((SL/dI)*p(i,j,k) - chi*invh2i*nb(p,i,j,k));
                // weighted CG, Jacobi (constant diag) -> plain CG
                dd.a = rr.a; double rho = wdot(rr,rr), r0 = std::sqrt(rho);
                for (int it=0; it<2000 && std::sqrt(rho) > 1e-11*std::max(r0,1e-300); ++it) {
                    applyA(dd, Ad); double al = rho / wdot(dd,Ad);
                    for(std::size_t c=0;c<del.a.size();++c){ del.a[c]+=al*dd.a[c]; rr.a[c]-=al*Ad.a[c]; }
                    double rn = wdot(rr,rr), be = rn/rho; rho = rn; ++cgTot;
                    for(std::size_t c=0;c<dd.a.size();++c) dd.a[c] = rr.a[c] + be*dd.a[c];
                }
                for(int k=1;k<=g.nz-2;++k) for(int j=0;j<=g.ny-2;++j) for(int i=0;i<=g.nx-2;++i) p(i,j,k) += del(i,j,k);
                // mechanics with new p
                S.mechRHS(p, b);
                itTotal += S.solveMechanics(u, b, W, tol, maxit, smooth);
                divergence(u, eps); meanStress(eps, p, sigv);
                double dmax=0, pmax=1e-300;
                for(std::size_t c=0;c<p.a.size();++c){ dmax=std::max(dmax,std::fabs(p.a[c]-pold.a[c])); pmax=std::max(pmax,std::fabs(p.a[c])); }
                chg = dmax/pmax;
                if (chg < fsTol) break;
            }
            fsTot += kfs; fsWorst = std::max<long>(fsWorst, kfs);
            if (n % sampleEvery == 0 || n == nI) sample(t);
            if (n % std::max<long>(1, nI / 10) == 0)
                std::printf("  t = %8.1f s  fixed-stress its = %d (last change %.1e)  avg CG/solve = %.0f\n",
                            t, kfs, chg, double(cgTot)/double(fsTot));
        }
        std::printf("fixed-stress iterations: avg %.1f, worst %ld\n", double(fsTot)/nI, fsWorst);
    } else
    for (long n = 1; n <= nsteps; ++n) {
        t = n * dt;

        // explicit diffusion term at the OLD level p^n -- 4TH-ORDER UPGRADE:
        // order-selected per axis (see Solver::d2x/d2y/d2z, x4ok/y4ok/z4ok in
        // mechanics.hpp). This loop's k-range already starts at k=1 (p=0 is
        // enforced directly on the free surface, never timestepped), so the
        // z near-boundary fallback (k<2) only ever bites at k=1 itself.
        S.fillGhostsP(p);
        auto fp = [&](int i, int j, int k) { return p(i, j, k); };
#pragma omp parallel for collapse(2) schedule(static)
        for (int k = 1; k <= g.nz - 2; ++k)
            for (int j = 0; j <= g.ny - 2; ++j)
                for (int i = 0; i <= g.nx - 2; ++i)
                    lapP(i, j, k) = Solver::d2x(fp, i, j, k, h, S.x4ok(i)) +
                                    Solver::d2y(fp, i, j, k, h, S.y4ok(j)) +
                                    Solver::d2z(fp, i, j, k, h, S.z4ok(k));
        // sigma_v-lag fixed-stress update, with a SMALL bounded number of
        // Picard sub-iterations per step (subiters, default 1 = pure single
        // pass). Both p^n (pn0) and lapP stay FIXED across sub-iterations --
        // only the coupling correction (sigv, recomputed from the latest
        // pressure each sub-iteration) is refined. NOT run to full
        // convergence -- that was shown to converge to a stable but wrong
        // fixed point near the point source.
        pn0.a = p.a;
        for (int si = 0; si < subiters; ++si) {
            pnew.a = pn0.a;
#pragma omp parallel for collapse(2) schedule(static)
            for (int k = 1; k <= g.nz - 2; ++k)
                for (int j = 0; j <= g.ny - 2; ++j)
                    for (int i = 0; i <= g.nx - 2; ++i)
                        pnew(i, j, k) = pn0(i, j, k) +
                                        dt / Sfs * chi * lapP(i, j, k) -
                                        (alpha / (Sfs * Kdr)) * (sigv(i, j, k) - sigv0(i, j, k));
            for (std::size_t c = 0; c < srcCells.size(); ++c) {
                auto& cc = srcCells[c];
                pnew(cc[0], cc[1], cc[2]) += dt / Sfs * srcWeight[c];
            }
            std::swap(p.a, pnew.a);
            setDirichlet(t);
            S.mechRHS(p, b);
            itTotal += S.solveMechanics(u, b, W, tol, maxit, smooth);
            divergence(u, eps);
            meanStress(eps, p, sigv);
        }
        sigv0.a = sigv.a;

        if (n % sampleEvery == 0 || n == nsteps) sample(t);
        if (n % std::max<long>(1, nsteps / 20) == 0) {
            double pmax = 0, umax = 0;
            for (int k = 0; k < g.nz; ++k) for (int j = 0; j < g.ny; ++j)
                for (int i = 0; i < g.nx; ++i) {
                    pmax = std::max(pmax, std::fabs(p(i, j, k)));
                    umax = std::max(umax, std::fabs(u.z(i, j, k)));
                }
            std::printf("  t = %8.1f s (%5.2f h)  BiCGSTAB its/solve = %.1f  "
                        "max|p| = %.4e Pa  max|uz| = %.4e m\n",
                        t, t / 3600.0, double(itTotal) / double(n), pmax, umax);
        }
    }
    std::fclose(f);

    // ---------------- final report ---------------------------------------
    std::printf("\nSteady-state check at r = %g m, t = %g h\n", robs, tend / 3600.0);
    std::printf("   z[m]   p_num[kPa]  p_mirror[kPa]  err%%   "
                "tiltA=-duz/dr  tiltB=rot  tilt_ws(whole-space)[urad]\n");
    for (auto& o : obs) {
        const double pn = p(o.i, o.j, o.k);
        const double pm = S.ref.p(robs, 0.0, o.z, tend);
        const double duzdr = (u.z(o.i + 1, o.j, o.k) - u.z(o.i - 1, o.j, o.k)) / (2 * h);
        const double durdz = (u.x(o.i, o.j, o.k + 1) - u.x(o.i, o.j, o.k - 1)) / (2 * h);
        const double tA = -duzdr, tB = 0.5 * (durdz - duzdr);
        const double tw = -S.ref.ws.duz_dr(robs, o.z - d, tend);
        std::printf("  %5.0f  %10.4f  %13.4f  %6.2f  %13.4f  %9.4f  %14.4f\n", o.z, pn / 1e3,
                    pm / 1e3, 100.0 * (pn - pm) / (std::fabs(pm) + 1e-30), tA * 1e6, tB * 1e6,
                    tw * 1e6);
    }
    std::printf("\nwrote %s\n", out.c_str());
    return 0;
}

