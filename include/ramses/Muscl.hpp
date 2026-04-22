#ifndef RAMSES_MUSCL_HPP
#define RAMSES_MUSCL_HPP

#include "Types.hpp"
#include <vector>

namespace ramses {

/**
 * @brief Implements MUSCL reconstruction and tracing logic.
 */
class Muscl {
public:
    /**
     * @brief Performs MUSCL prediction (tracing) for a single cell and dimension.
     * @param q Cell-centered primitive state
     * @param dq Limited slope for this cell/dimension
     * @param s0 Source term (predicted change over dt/2)
     * @param dt_dx dt/dx factor
     * @param qm Output left-interface state (predicted)
     * @param qp Output right-interface state (predicted)
     */
    static void predict(const real_t q[], const real_t dq[], const real_t s0[], 
                       real_t dt_dx, real_t qm[], real_t qp[], int nvar);
};

} // namespace ramses

#endif // RAMSES_MUSCL_HPP
