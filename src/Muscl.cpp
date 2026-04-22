#include "ramses/Muscl.hpp"
#include <algorithm>

namespace ramses {

void Muscl::predict(const real_t q[], const real_t dq[], const real_t s0[], 
                   real_t dt_dx, real_t qm[], real_t qp[], int nvar) {
    const real_t smallr = 1e-10;
    
    // MUSCL Prediction (Hancock):
    // q_left  = q + 0.5 * dq + 0.5 * s0 * dt/dx
    // q_right = q - 0.5 * dq + 0.5 * s0 * dt/dx
    
    // In RAMSES:
    // qp (right state at left interface) = q - 0.5*dq + 0.5*s0*dtdx
    // qm (left state at right interface) = q + 0.5*dq + 0.5*s0*dtdx

    for (int iv = 0; iv < nvar; ++iv) {
        real_t src_term = s0[iv] * dt_dx * 0.5;
        qp[iv] = q[iv] - 0.5 * dq[iv] + src_term;
        qm[iv] = q[iv] + 0.5 * dq[iv] + src_term;
    }

    // Density floors
    if (qp[0] < smallr) qp[0] = q[0];
    if (qm[0] < smallr) qm[0] = q[0];
}

} // namespace ramses
