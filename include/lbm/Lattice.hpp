#pragma once

// D2Q9 lattice constants and the moment basis used by every collision model.
// Lattice units throughout: node spacing = time step = 1, sound speed c_s^2 = 1/3.

#include <cstdint>

namespace lbm {

// Population storage precision. The website build defines LBM_SINGLE_PRECISION;
// validation runs use double. Collision arithmetic is always done in double.
#ifdef LBM_SINGLE_PRECISION
using real = float;
#else
using real = double;
#endif

constexpr int Q = 9;

// Lattice velocities c_i = (cx[i], cy[i]); opposite[i] has c = -c_i.
//   6 2 5
//   3 0 1
//   7 4 8      (y grows downward on the page, matching the river's reference frame)
constexpr int cx[Q] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
constexpr int cy[Q] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
constexpr int opposite[Q] = {0, 3, 4, 1, 2, 7, 8, 5, 6};
constexpr double w[Q] = {4.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9,
                         1.0 / 36, 1.0 / 36, 1.0 / 36, 1.0 / 36};
constexpr double cs2 = 1.0 / 3.0;
constexpr double pi = 3.14159265358979323846;

// Moment basis e_k(c_i) of Dünweg, Schiller & Ladd (2007). The rows are
// orthogonal under the lattice weights, sum_i w_i e_k(c_i) e_l(c_i) = norm[k] delta_kl,
// so the inverse transform is f_i = w_i sum_k e_k(c_i) m_k / norm[k].
//   k = 0  1                     density
//   k = 1  cx,  k = 2  cy        momentum
//   k = 3  3c^2 - 2              bulk stress (trace)
//   k = 4  cx^2 - cy^2           shear stress, normal difference
//   k = 5  cx cy                 shear stress, off-diagonal
//   k = 6  (3c^2 - 4) cx         ghost (odd)
//   k = 7  (3c^2 - 4) cy         ghost (odd)
//   k = 8  9c^4 - 15c^2 + 2      ghost (even)
constexpr double basis[Q][Q] = {
    //  0    1    2    3    4    5    6    7    8   (direction i)
    {   1,   1,   1,   1,   1,   1,   1,   1,   1},
    {   0,   1,   0,  -1,   0,   1,  -1,  -1,   1},
    {   0,   0,   1,   0,  -1,   1,   1,  -1,  -1},
    {  -2,   1,   1,   1,   1,   4,   4,   4,   4},
    {   0,   1,  -1,   1,  -1,   0,   0,   0,   0},
    {   0,   0,   0,   0,   0,   1,  -1,   1,  -1},
    {   0,  -1,   0,   1,   0,   2,  -2,  -2,   2},
    {   0,   0,  -1,   0,   1,   2,   2,  -2,  -2},
    {   2,  -4,  -4,  -4,  -4,   8,   8,   8,   8},
};
constexpr double norm[Q] = {1.0, 1.0 / 3, 1.0 / 3, 4.0, 4.0 / 9, 1.0 / 9,
                            2.0 / 3, 2.0 / 3, 16.0};

// m_k = sum_i e_k(c_i) f_i, written out from the basis table above (the kernel test
// checks it against the table).
inline void toMoments(const double f[Q], double m[Q]) {
    const double axes = f[1] + f[2] + f[3] + f[4];
    const double diagonals = f[5] + f[6] + f[7] + f[8];
    m[0] = f[0] + axes + diagonals;
    m[1] = f[1] - f[3] + f[5] - f[6] - f[7] + f[8];
    m[2] = f[2] - f[4] + f[5] + f[6] - f[7] - f[8];
    m[3] = -2.0 * f[0] + axes + 4.0 * diagonals;
    m[4] = f[1] - f[2] + f[3] - f[4];
    m[5] = f[5] - f[6] + f[7] - f[8];
    m[6] = -f[1] + f[3] + 2.0 * (f[5] - f[6] - f[7] + f[8]);
    m[7] = -f[2] + f[4] + 2.0 * (f[5] + f[6] - f[7] - f[8]);
    m[8] = 2.0 * f[0] - 4.0 * axes + 8.0 * diagonals;
}

// f_i = w_i sum_k e_k(c_i) m_k / norm[k], written out likewise.
inline void fromMoments(const double m[Q], double f[Q]) {
    const double n0 = m[0], n1 = 3.0 * m[1], n2 = 3.0 * m[2], n3 = 0.25 * m[3];
    const double n4 = 2.25 * m[4], n5 = 9.0 * m[5], n6 = 1.5 * m[6], n7 = 1.5 * m[7];
    const double n8 = m[8] / 16.0;
    const double axisCommon = n0 + n3 - 4.0 * n8;
    const double diagCommon = n0 + 4.0 * n3 + 8.0 * n8;
    f[0] = (4.0 / 9.0) * (n0 - 2.0 * n3 + 2.0 * n8);
    f[1] = (1.0 / 9.0) * (axisCommon + n1 + n4 - n6);
    f[2] = (1.0 / 9.0) * (axisCommon + n2 - n4 - n7);
    f[3] = (1.0 / 9.0) * (axisCommon - n1 + n4 + n6);
    f[4] = (1.0 / 9.0) * (axisCommon - n2 - n4 + n7);
    f[5] = (1.0 / 36.0) * (diagCommon + n1 + n2 + n5 + 2.0 * n6 + 2.0 * n7);
    f[6] = (1.0 / 36.0) * (diagCommon - n1 + n2 - n5 - 2.0 * n6 + 2.0 * n7);
    f[7] = (1.0 / 36.0) * (diagCommon - n1 - n2 + n5 - 2.0 * n6 - 2.0 * n7);
    f[8] = (1.0 / 36.0) * (diagCommon + n1 - n2 - n5 + 2.0 * n6 - 2.0 * n7);
}

// Second-order equilibrium in population space:
// f_i^eq = w_i rho [1 + 3 c.u + 4.5 (c.u)^2 - 1.5 u^2]
inline double equilibrium(int i, double rho, double ux, double uy) {
    const double cu = cx[i] * ux + cy[i] * uy;
    return w[i] * rho * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * (ux * ux + uy * uy));
}

// Equilibrium moments of the non-conserved modes (k = 3..8); the ghosts vanish.
struct EquilibriumMoments {
    double bulk, normalDiff, shear;   // k = 3, 4, 5
};
inline EquilibriumMoments equilibriumMoments(double rho, double ux, double uy) {
    return {3.0 * rho * (ux * ux + uy * uy), rho * (ux * ux - uy * uy), rho * ux * uy};
}

}  // namespace lbm
