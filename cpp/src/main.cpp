#include <iostream>
#include <vector>
#include "ramses/SlopeLimiter.hpp"

int main() {
    std::cout << "--- Testing Slope Limiter ---" << std::endl;
    
    using ramses::real_t;
    real_t ql = 1.0, qc = 2.0, qr = 4.0;
    
    // MinMod (slope_type 1)
    real_t slope_minmod = ramses::SlopeLimiter::compute_slope(ql, qc, qr, 1);
    // MonCen (slope_type 2)
    real_t slope_moncen = ramses::SlopeLimiter::compute_slope(ql, qc, qr, 2);
    
    std::cout << "Values: ql=" << ql << " qc=" << qc << " qr=" << qr << std::endl;
    std::cout << "  MinMod: " << slope_minmod << " (Expected 1.0)" << std::endl;
    std::cout << "  MonCen: " << slope_moncen << " (Expected 1.5)" << std::endl;

    return 0;
}
