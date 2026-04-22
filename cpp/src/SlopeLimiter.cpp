#include "ramses/SlopeLimiter.hpp"
#include <cmath>
#include <algorithm>

namespace ramses {

static inline real_t sign(real_t a, real_t b) {
    return b >= 0 ? std::abs(a) : -std::abs(a);
}

real_t SlopeLimiter::compute_slope(real_t ql, real_t qc, real_t qr, int slope_type, real_t theta) {
    if (slope_type == 0) return 0.0;

    real_t dlft = qc - ql;
    real_t drgt = qr - qc;

    if (slope_type >= 1 && slope_type <= 3) {
        real_t mult = static_cast<real_t>(std::min(slope_type, 2));
        real_t dl = mult * dlft;
        real_t dr = mult * drgt;
        real_t dcen = 0.5 * (dl + dr) / mult;
        real_t dsgn = sign(1.0, dcen);
        real_t dlim = std::min(std::abs(dl), std::abs(dr));
        if (dlft * drgt <= 0.0) dlim = 0.0;
        return dsgn * std::min(dlim, std::abs(dcen));
    } 
    else if (slope_type == 7) { // van Leer
        if (dlft * drgt <= 0.0) return 0.0;
        return (2.0 * dlft * drgt) / (dlft + drgt);
    }
    else if (slope_type == 8) { // Generalized MonCen
        real_t dcen = 0.5 * (dlft + drgt);
        real_t dsgn = sign(1.0, dcen);
        real_t dlim = std::min(theta * std::abs(dlft), theta * std::abs(drgt));
        if (dlft * drgt <= 0.0) dlim = 0.0;
        return dsgn * std::min(dlim, std::abs(dcen));
    }

    return 0.0; // Default or unknown
}

} // namespace ramses
