#include "ramses/solvers/hydro/RiemannSolver.hpp"
#include "ramses/solvers/physics/EquationOfState.hpp"
#include "ramses/core/Parameters.hpp"
#include "ramses/core/Constants.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>

namespace ramses {

real_t RiemannSolver::get_cs2(real_t d, real_t p, real_t gamma, const real_t q[], int nener, const std::vector<real_t>& gamma_rad) {
    return EquationOfState::get_cs2(d, p, gamma, q, nener, gamma_rad);
}

void RiemannSolver::solve_llf(const real_t ql[], const real_t qr[], real_t flux[], real_t gamma, int nener, const std::vector<real_t>& gamma_rad) {
    real_t fl[20], fr[20];
    compute_flux(ql, fl, gamma, nener, gamma_rad);
    compute_flux(qr, fr, gamma, nener, gamma_rad);
    real_t ul[20], ur[20];
    prim_to_cons(ql, ul, gamma, nener, gamma_rad);
    prim_to_cons(qr, ur, gamma, nener, gamma_rad);

    int ipress = NDIM + 1;
    real_t al = std::sqrt(get_cs2(ql[0], ql[ipress], gamma, ql, nener, gamma_rad));
    real_t ar = std::sqrt(get_cs2(qr[0], qr[ipress], gamma, qr, nener, gamma_rad));
    real_t a_max = std::max(std::abs(ql[1]) + al, std::abs(qr[1]) + ar);
    for (int i = 0; i < NDIM + 2; ++i) {
        flux[i] = 0.5 * (fl[i] + fr[i]) - 0.5 * a_max * (ur[i] - ul[i]);
    }
}

void RiemannSolver::solve_hll(const real_t ql[], const real_t qr[], real_t flux[], real_t gamma, int nener, const std::vector<real_t>& gamma_rad) {
    // 1. Compute flux vectors for the left and right states
    real_t fl[20], fr[20];
    compute_flux(ql, fl, gamma, nener, gamma_rad);
    compute_flux(qr, fr, gamma, nener, gamma_rad);

    // 2. Compute conserved variable states for the left and right states
    real_t ul_vec[20], ur_vec[20];
    prim_to_cons(ql, ul_vec, gamma, nener, gamma_rad);
    prim_to_cons(qr, ur_vec, gamma, nener, gamma_rad);

    // 3. Compute sound speeds for the left and right states
    int ipress = NDIM + 1;
    real_t al = std::sqrt(get_cs2(ql[0], ql[ipress], gamma, ql, nener, gamma_rad));
    real_t ar = std::sqrt(get_cs2(qr[0], qr[ipress], gamma, qr, nener, gamma_rad));

    // 4. Compute left and right wave propagation speed bounds:
    //    S_L = min(0, u_L - c_L, u_R - c_R)
    //    S_R = max(0, u_L + c_L, u_R + c_R)
    real_t sl = std::min(0.0, std::min(ql[1] - al, qr[1] - ar));
    real_t sr = std::max(0.0, std::max(ql[1] + al, qr[1] + ar));

    // 5. Evaluate the HLL flux formula:
    //    F_HLL = (S_R * F_L - S_L * F_R + S_L * S_R * (U_R - U_L)) / (S_R - S_L)
    if (sr - sl < 1e-20) {
        for (int i = 0; i < NDIM + 2; ++i) {
            flux[i] = 0.5 * (fl[i] + fr[i]);
        }
    } else {
        real_t inv_srs = 1.0 / (sr - sl);
        for (int i = 0; i < NDIM + 2; ++i) {
            flux[i] = (sr * fl[i] - sl * fr[i] + sr * sl * (ur_vec[i] - ul_vec[i])) * inv_srs;
        }
    }
}

void RiemannSolver::solve_hllc(const real_t ql[], const real_t qr[], real_t flux[], real_t gamma, int nener, const std::vector<real_t>& gamma_rad) {
    const real_t smallc = 1e-10;
    const real_t smallr = smallc;
    const real_t smallp = smallc * smallc / gamma;
    const real_t entho = 1.0 / (gamma - 1.0);
    int ipress = NDIM + 1;

    // Left variables
    real_t rl = std::max(ql[0], smallr);
    real_t pl = std::max(ql[ipress], rl * smallp);
    real_t ul = ql[1];

    real_t el = pl * entho;
    real_t ecinl = 0.5 * rl * ul * ul;
    for (int idim = 2; idim <= NDIM; ++idim) {
        ecinl += 0.5 * rl * ql[idim] * ql[idim];
    }
    real_t etotl = el + ecinl;
    real_t eradl[20] = {0};
    real_t ptotl = pl;
    for (int ie = 0; ie < nener; ++ie) {
        eradl[ie] = ql[ipress + 1 + ie] / (gamma_rad[ie] - 1.0);
        etotl += eradl[ie];
        ptotl += ql[ipress + 1 + ie];
    }

    // Right variables
    real_t rr = std::max(qr[0], smallr);
    real_t pr = std::max(qr[ipress], rr * smallp);
    real_t ur = qr[1];

    real_t er = pr * entho;
    real_t ecinr = 0.5 * rr * ur * ur;
    for (int idim = 2; idim <= NDIM; ++idim) {
        ecinr += 0.5 * rr * qr[idim] * qr[idim];
    }
    real_t etotr = er + ecinr;
    real_t eradr[20] = {0};
    real_t ptotr = pr;
    for (int ie = 0; ie < nener; ++ie) {
        eradr[ie] = qr[ipress + 1 + ie] / (gamma_rad[ie] - 1.0);
        etotr += eradr[ie];
        ptotr += qr[ipress + 1 + ie];
    }

    // Largest eigenvalues
    real_t cfastl = gamma * pl;
    for (int ie = 0; ie < nener; ++ie) {
        cfastl += gamma_rad[ie] * ql[ipress + 1 + ie];
    }
    cfastl = std::sqrt(std::max(cfastl / rl, smallc * smallc));

    real_t cfastr = gamma * pr;
    for (int ie = 0; ie < nener; ++ie) {
        cfastr += gamma_rad[ie] * qr[ipress + 1 + ie];
    }
    cfastr = std::sqrt(std::max(cfastr / rr, smallc * smallc));

    // HLL wave speeds
    real_t sl = std::min(ul, ur) - std::max(cfastl, cfastr);
    real_t sr = std::max(ul, ur) + std::max(cfastl, cfastr);

    // Lagrangian sound speeds
    real_t rcl = rl * (ul - sl);
    real_t rcr = rr * (sr - ur);

    // Acoustic star state
    real_t ustar = (rcr * ur + rcl * ul + (ptotl - ptotr)) / (rcr + rcl);
    real_t ptotstar = (rcr * ptotl + rcl * ptotr + rcl * rcr * (ul - ur)) / (rcr + rcl);

    // Left star region variables
    real_t rstarl = rl * (sl - ul) / (sl - ustar);
    real_t etotstarl = ((sl - ul) * etotl - ptotl * ul + ptotstar * ustar) / (sl - ustar);
    real_t eradstarl[20] = {0};
    for (int ie = 0; ie < nener; ++ie) {
        eradstarl[ie] = eradl[ie] * (sl - ul) / (sl - ustar);
    }

    // Right star region variables
    real_t rstarr = rr * (sr - ur) / (sr - ustar);
    real_t etotstarr = ((sr - ur) * etotr - ptotr * ur + ptotstar * ustar) / (sr - ustar);
    real_t eradstarr[20] = {0};
    for (int ie = 0; ie < nener; ++ie) {
        eradstarr[ie] = eradr[ie] * (sr - ur) / (sr - ustar);
    }

    // Sample solution at x/t = 0
    real_t ro, uo, ptoto, etoto;
    real_t erado[20] = {0};
    if (sl > 0.0) {
        ro = rl; uo = ul; ptoto = ptotl; etoto = etotl;
        for (int ie = 0; ie < nener; ++ie) erado[ie] = eradl[ie];
    } else if (ustar > 0.0) {
        ro = rstarl; uo = ustar; ptoto = ptotstar; etoto = etotstarl;
        for (int ie = 0; ie < nener; ++ie) erado[ie] = eradstarl[ie];
    } else if (sr > 0.0) {
        ro = rstarr; uo = ustar; ptoto = ptotstar; etoto = etotstarr;
        for (int ie = 0; ie < nener; ++ie) erado[ie] = eradstarr[ie];
    } else {
        ro = rr; uo = ur; ptoto = ptotr; etoto = etotr;
        for (int ie = 0; ie < nener; ++ie) erado[ie] = eradr[ie];
    }

    // Compute Godunov flux
    flux[0] = ro * uo;
    flux[1] = ro * uo * uo + ptoto;
    for (int idim = 2; idim <= NDIM; ++idim) {
        flux[idim] = (ustar > 0.0) ? (ro * uo * ql[idim]) : (ro * uo * qr[idim]);
    }
    flux[ipress] = (etoto + ptoto) * uo;
    for (int ie = 0; ie < nener; ++ie) {
        flux[ipress + 1 + ie] = uo * erado[ie];
    }
}

void RiemannSolver::solve_godunov_nr(const real_t ql[], const real_t qr[], real_t flux[], real_t gamma) {
    const real_t smallc = 1e-10;
    const real_t smallr = smallc;
    const real_t smallp = smallc * smallc / gamma;
    const int niter = 10;
    const real_t gamma6 = (gamma + 1.0) / (2.0 * gamma);

    int ipress = NDIM + 1;
    real_t rl = std::max(ql[0], smallr);
    real_t ul = ql[1];
    real_t pl = std::max(ql[ipress], rl * smallp);

    real_t rr = std::max(qr[0], smallr);
    real_t ur = qr[1];
    real_t pr = std::max(qr[ipress], rr * smallp);

    real_t cl_sq = gamma * pl * rl;
    real_t cr_sq = gamma * pr * rr;
    real_t wl = std::sqrt(std::max(cl_sq, 0.0));
    real_t wr = std::sqrt(std::max(cr_sq, 0.0));

    real_t pstar = (wr * pl + wl * pr + wl * wr * (ul - ur)) / (wl + wr + 1e-20);
    pstar = std::max(pstar, 0.0);
    real_t pold = pstar;

    // Newton-Raphson iterations to find pstar
    for (int iter = 0; iter < niter; ++iter) {
        real_t wwl_sq = cl_sq * (1.0 + gamma6 * (pold - pl) / (pl + 1e-20));
        real_t wwr_sq = cr_sq * (1.0 + gamma6 * (pold - pr) / (pr + 1e-20));
        real_t wwl = std::sqrt(std::max(wwl_sq, 0.0));
        real_t wwr = std::sqrt(std::max(wwr_sq, 0.0));

        real_t ql_nr = (2.0 * std::pow(wwl, 3)) / (wwl * wwl + cl_sq + 1e-20);
        real_t qr_nr = (2.0 * std::pow(wwr, 3)) / (wwr * wwr + cr_sq + 1e-20);
        real_t usl = ul - (pold - pl) / (wwl + 1e-20);
        real_t usr = ur + (pold - pr) / (wwr + 1e-20);

        real_t delp = (qr_nr * ql_nr / (qr_nr + ql_nr + 1e-20)) * (usl - usr);
        delp = std::max(delp, -pold);
        pold = pold + delp;
        if (std::abs(delp / (pold + 1e-10)) < 1e-6) break;
    }
    pstar = pold;


    real_t wwl_f = std::sqrt(std::max(cl_sq * (1.0 + gamma6 * (pstar - pl) / (pl + 1e-20)), 0.0));
    real_t wwr_f = std::sqrt(std::max(cr_sq * (1.0 + gamma6 * (pstar - pr) / (pr + 1e-20)), 0.0));
    real_t ustar = 0.5 * (ul + (pl - pstar) / (wwl_f + 1e-20) + ur - (pr - pstar) / (wwr_f + 1e-20));

    real_t ro, uo, po;
    real_t sgnm = (ustar > 0) ? 1.0 : -1.0;

    if (sgnm == 1.0) {
        ro = rl; uo = ul; po = pl;
    } else {
        ro = rr; uo = ur; po = pr;
    }

    real_t co = std::sqrt(std::max(gamma * po / (ro + 1e-20), 1e-20));
    real_t rstar_val;
    if (pstar >= po) {
        real_t cl_unperturbed = (sgnm == 1.0) ? cl_sq : cr_sq;
        rstar_val = ro * ((gamma + 1.0) * pstar + (gamma - 1.0) * po) / ((gamma - 1.0) * pstar + (gamma + 1.0) * po);
    } else {
        rstar_val = ro * std::pow(pstar / (po + 1e-20), 1.0 / gamma);
    }
    rstar_val = std::max(rstar_val, smallr);

    real_t cstar = std::sqrt(std::max(gamma * pstar / (rstar_val + 1e-20), 1e-20));
    real_t spout = co - sgnm * uo;
    real_t spin = cstar - sgnm * ustar;

    if (pstar >= po) {
        real_t wo = (sgnm == 1.0) ? wwl_f : wwr_f;
        real_t ushock = wo / (ro + 1e-20) - sgnm * uo;
        spout = ushock;
        spin = ushock;
    }

    // Use HLLC-style wave speed estimates for robust sampling
    real_t cl_hllc = std::sqrt(get_cs2(rl, pl, gamma));
    real_t cr_hllc = std::sqrt(get_cs2(rr, pr, gamma));
    real_t sl_hllc = std::min(ul, ur) - std::max(cl_hllc, cr_hllc);
    real_t sr_hllc = std::max(ul, ur) + std::max(cl_hllc, cr_hllc);

    if (sl_hllc > 0.0) {
        ro = rl; uo = ul; po = pl;
    } else if (sr_hllc < 0.0) {
        ro = rr; uo = ur; po = pr;
    } else {
        ro = rstar_val; uo = ustar; po = pstar;
    }

    if (std::abs(pl - pr) > 1e-3) {
        static int count_nr = 0; if(count_nr++ < 5) {
            std::cout << "[NR DEBUG] pstar=" << pstar << " ro=" << ro << " uo=" << uo << " po=" << po << " spout=" << spout << " spin=" << spin << std::endl;
        }
    }

    real_t e_int = po / (gamma - 1.0);
    real_t e_kin = 0.5 * ro * uo * uo;
    for(int i=2; i<=NDIM; ++i) {
        real_t vt = (ustar > 0.0) ? ql[i] : qr[i];
        e_kin += 0.5 * ro * vt * vt;
    }
    real_t E = e_int + e_kin;

    flux[0] = ro * uo;
    flux[1] = ro * uo * uo + po;
    for(int i=2; i<=NDIM; ++i) {
        flux[i] = flux[0] * ((ustar > 0.0) ? ql[i] : qr[i]);
    }
    flux[ipress] = (E + po) * uo;
}

void RiemannSolver::solve_acoustic(const real_t ql[], const real_t qr[], real_t flux[], real_t gamma, int nener, const std::vector<real_t>& gamma_rad) {
    const real_t smallc = 1e-10;
    const real_t smallr = smallc;
    const real_t smallp = smallc * smallc / gamma;

    int ipress = NDIM + 1;
    real_t rl = std::max(ql[0], smallr);
    real_t ul = ql[1];
    real_t pl = std::max(ql[ipress], rl * smallp);
    for (int ie = 0; ie < nener; ++ie) pl += ql[ipress + 1 + ie];

    real_t rr = std::max(qr[0], smallr);
    real_t ur = qr[1];
    real_t pr = std::max(qr[ipress], rr * smallp);
    for (int ie = 0; ie < nener; ++ie) pr += qr[ipress + 1 + ie];

    real_t cl = std::sqrt(get_cs2(rl, ql[ipress], gamma, ql, nener, gamma_rad));
    real_t cr = std::sqrt(get_cs2(rr, qr[ipress], gamma, qr, nener, gamma_rad));
    real_t wl = cl * rl;
    real_t wr = cr * rr;

    real_t pstar = ((wr * pl + wl * pr) + wl * wr * (ul - ur)) / (wl + wr + 1e-20);
    real_t ustar = ((wr * ur + wl * ul) + (pl - pr)) / (wl + wr + 1e-20);

    real_t sgnm = (ustar >= 0.0) ? 1.0 : -1.0;
    real_t ro, uo, po, wo, co;
    if (sgnm == 1.0) {
        ro = rl; uo = ul; po = pl; wo = wl; co = cl;
    } else {
        ro = rr; uo = ur; po = pr; wo = wr; co = cr;
    }

    real_t rstar = ro + (pstar - po) / std::max(co * co, 1e-20);
    rstar = std::max(rstar, smallr);
    real_t cstar = std::sqrt(std::max(gamma * pstar / rstar, 1e-20));
    cstar = std::max(cstar, smallc);

    real_t spout = co - sgnm * uo;
    real_t spin = cstar - sgnm * ustar;
    real_t ushock = 0.5 * (spin + spout);
    ushock = std::max(ushock, -sgnm * ustar);

    if (pstar >= po) {
        spout = ushock;
        spin = spout;
    }

    real_t rgdnv, ugdnv, pgdnv;
    if (spout < 0.0) {
        rgdnv = ro; ugdnv = uo; pgdnv = po;
    } else if (spin >= 0.0) {
        rgdnv = rstar; ugdnv = ustar; pgdnv = pstar;
    } else {
        real_t frac = spout / (spout - spin + 1e-20);
        rgdnv = frac * rstar + (1.0 - frac) * ro;
        ugdnv = frac * ustar + (1.0 - frac) * uo;
        pgdnv = frac * pstar + (1.0 - frac) * po;
    }

    real_t qgdnv[20] = {0};
    qgdnv[0] = rgdnv; qgdnv[1] = ugdnv;
    for (int i = 2; i <= NDIM; ++i) {
        qgdnv[i] = (sgnm == 1.0) ? ql[i] : qr[i];
    }

    real_t sum_prad_star = 0.0;
    for (int ie = 0; ie < nener; ++ie) {
        real_t prad_star = (sgnm == 1.0) ? ql[ipress + 1 + ie] : qr[ipress + 1 + ie];
        qgdnv[ipress + 1 + ie] = prad_star;
        sum_prad_star += prad_star;
    }
    qgdnv[ipress] = std::max(pgdnv - sum_prad_star, rgdnv * smallp);

    compute_flux(qgdnv, flux, gamma, nener, gamma_rad);
}

void RiemannSolver::prim_to_cons(const real_t q[], real_t u[], real_t gamma, int nener, const std::vector<real_t>& gamma_rad) {
    u[0] = q[0];
    real_t v2 = 0;
    for (int i = 1; i <= NDIM; ++i) {
        u[i] = q[0] * q[i];
        v2 += q[i] * q[i];
    }
    int ipress = NDIM + 1;
    real_t e_kin = 0.5 * q[0] * v2;
    real_t e_int = q[ipress] / (gamma - 1.0);
    real_t e_nonthermal = 0.0;
    for (int ie = 0; ie < nener; ++ie) {
        e_nonthermal += q[ipress + 1 + ie] / (gamma_rad[ie] - 1.0);
    }
    if (params::barotropic_eos && params::barotropic_eos_form == "isothermal") {
    }
    u[ipress] = e_kin + e_int + e_nonthermal;
}

void RiemannSolver::compute_flux(const real_t q[], real_t f[], real_t gamma, int nener, const std::vector<real_t>& gamma_rad) {
    real_t rho = q[0]; real_t u = q[1];
    int ipress = NDIM + 1;
    real_t p = q[ipress];

    real_t v2 = 0;
    for(int i=1; i<=NDIM; ++i) v2 += q[i]*q[i];

    real_t e_kin = 0.5 * rho * v2;
    real_t e_int = p / (gamma - 1.0);
    real_t e_nonthermal = 0.0;
    real_t p_tot = p;
    for (int ie = 0; ie < nener; ++ie) {
        real_t prad = q[ipress + 1 + ie];
        p_tot += prad;
        e_nonthermal += prad / (gamma_rad[ie] - 1.0);
    }
    real_t E = e_kin + e_int + e_nonthermal;

    f[0] = rho * u;
    f[1] = rho * u * u + p_tot;
    for(int i=2; i<=NDIM; ++i) f[i] = rho * u * q[i];
    f[ipress] = (E + p_tot) * u;
}

} // namespace ramses
