#include "ramses/RiemannSolver.hpp"
#include <cmath>
#include <algorithm>

namespace ramses {

void RiemannSolver::solve_llf(const real_t ql[], const real_t qr[], real_t flux[], real_t gamma) {
    // q = [rho, u, v, w, p]
    // Cons = [rho, rho*u, rho*v, rho*w, E]
    
    real_t fl[5], fr[5];
    compute_flux(ql, fl, gamma);
    compute_flux(qr, fr, gamma);

    real_t ul[5], ur[5];
    prim_to_cons(ql, ul, gamma);
    prim_to_cons(qr, ur, gamma);

    // Max wave speed
    real_t al = std::sqrt(gamma * ql[4] / ql[0]);
    real_t ar = std::sqrt(gamma * qr[4] / qr[0]);
    real_t a_max = std::max(std::abs(ql[1]) + al, std::abs(qr[1]) + ar);

    // LLF Flux: F = 0.5*(Fl + Fr) - 0.5*a_max*(Ur - Ul)
    for (int i = 0; i < 5; ++i) {
        flux[i] = 0.5 * (fl[i] + fr[i]) - 0.5 * a_max * (ur[i] - ul[i]);
    }
}

void RiemannSolver::prim_to_cons(const real_t q[], real_t u[], real_t gamma) {
    u[0] = q[0]; // rho
    u[1] = q[0] * q[1]; // rho*u
    u[2] = q[0] * q[2]; // rho*v
    u[3] = q[0] * q[3]; // rho*w
    real_t e_kin = 0.5 * q[0] * (q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
    real_t e_int = q[4] / (gamma - 1.0);
    u[4] = e_kin + e_int; // E
}

void RiemannSolver::compute_flux(const real_t q[], real_t f[], real_t gamma) {
    real_t rho = q[0];
    real_t u = q[1];
    real_t v = q[2];
    real_t w = q[3];
    real_t p = q[4];
    
    real_t e_kin = 0.5 * rho * (u*u + v*v + w*w);
    real_t e_int = p / (gamma - 1.0);
    real_t E = e_kin + e_int;
    
    f[0] = rho * u;
    f[1] = rho * u * u + p;
    f[2] = rho * u * v;
    f[3] = rho * u * w;
    f[4] = (E + p) * u;
}

} // namespace ramses
