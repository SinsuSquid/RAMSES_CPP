#include "ramses/RamsesReader.hpp"
#include <iostream>
#include <algorithm>
#include <cmath>

namespace ramses {

bool RamsesReader::load_amr(AmrGrid& grid) {
    if (!file_.is_open()) return false;

    // 1. Grid Variables (Each is a separate record)
    int ncpu = read_single<int>();
    int ndim = read_single<int>();
    
    std::vector<int> nxyz;
    read_record(nxyz);
    int nx = nxyz[0], ny = nxyz[1], nz = nxyz[2];
    
    int nlevelmax = read_single<int>();
    int ngridmax = read_single<int>();
    int nboundary = read_single<int>();
    int ngrid_current = read_single<int>();
    real_t boxlen = read_single<real_t>();

    std::cout << "  [Reader] ncpu=" << ncpu << " ngrid_current=" << ngrid_current << " ncell_calc=" << (nx*ny*nz + 8*ngridmax) << std::endl;

    // Allocate grid
    grid.allocate(nx, ny, nz, ngridmax, 5, ncpu, nlevelmax);
    grid.ndim = ndim;

    // 2. Time Variables
    for(int i=0; i<11; ++i) skip_record(); 

    // 3. Level/List Variables
    std::vector<int> buf;
    read_record(buf); // headl
    for(int i=0; i<ncpu*nlevelmax; ++i) grid.headl.data()[i] = buf[i];
    
    read_record(buf); // taill
    for(int i=0; i<ncpu*nlevelmax; ++i) grid.taill.data()[i] = buf[i];

    read_record(buf); // numbl
    for(int i=0; i<ncpu*nlevelmax; ++i) grid.numbl.data()[i] = buf[i];

    skip_record(); // numbtot
    if (nboundary > 0) { skip_record(); skip_record(); skip_record(); }

    // 4. Free Memory / Ordering
    skip_record(); skip_record(); skip_record();

    // 5. Coarse Level
    read_record(buf); // son
    for(int i=0; i<grid.ncoarse; ++i) grid.son[i+1] = buf[i];
    read_record(buf); // flag1
    for(int i=0; i<grid.ncoarse; ++i) grid.flag1[i+1] = buf[i];
    read_record(buf); // cpu_map
    for(int i=0; i<grid.ncoarse; ++i) grid.cpu_map[i+1] = buf[i];

    return true;
}

bool RamsesReader::load_hydro(AmrGrid& grid) {
    return true; 
}

} // namespace ramses
