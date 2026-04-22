#include "ramses/RamsesWriter.hpp"
#include "ramses/Parameters.hpp"
#include <iostream>

namespace ramses {

void RamsesWriter::write_amr(const AmrGrid& grid) {
    if (!file_.is_open()) return;

    std::cout << "  [Writer] Saving Snapshot: ncpu=" << grid.ncpu << " ncoarse=" << grid.ncoarse << std::endl;
    // 0. Rank record (Fortran writes ncpu first)
    write_single<int>(grid.ncpu);

    // 1. Grid Variables Header
    uint32_t header_tag = 9 * 4 + sizeof(real_t);
    write_tag(header_tag);
    file_.write(reinterpret_cast<const char*>(&grid.ncpu), 4);
    int ndim = NDIM; file_.write(reinterpret_cast<const char*>(&ndim), 4);
    
    // Access global params for nx, ny, nz
    file_.write(reinterpret_cast<const char*>(&params::nx), 4);
    file_.write(reinterpret_cast<const char*>(&params::ny), 4);
    file_.write(reinterpret_cast<const char*>(&params::nz), 4);
    
    file_.write(reinterpret_cast<const char*>(&grid.nlevelmax), 4);
    file_.write(reinterpret_cast<const char*>(&grid.ngridmax), 4);
    int nboundary = 0; file_.write(reinterpret_cast<const char*>(&nboundary), 4);
    int ngrid_current = 0;
    file_.write(reinterpret_cast<const char*>(&ngrid_current), 4);
    real_t boxlen = 1.0; file_.write(reinterpret_cast<const char*>(&boxlen), sizeof(real_t));
    write_tag(header_tag);

    // 2. Dummy Time Variables (11 records to match reader)
    for(int i=0; i<11; ++i) write_single<int>(0);

    // 3. Level/List Variables
    std::vector<int> buf(grid.ncpu * grid.nlevelmax);
    for(size_t i=0; i<buf.size(); ++i) buf[i] = grid.headl.data()[i];
    write_record(buf);
    for(size_t i=0; i<buf.size(); ++i) buf[i] = grid.taill.data()[i];
    write_record(buf);
    for(size_t i=0; i<buf.size(); ++i) buf[i] = grid.numbl.data()[i];
    write_record(buf);
    
    write_single<int>(0); // numbtot

    // 4. Free Memory / Ordering dummies
    write_single<int>(0); 
    write_single<int>(0); 
    write_single<int>(0);

    // 5. Coarse Level
    std::vector<int> c_buf(grid.ncoarse);
    for(int i=0; i<grid.ncoarse; ++i) c_buf[i] = grid.son[i+1];
    write_record(c_buf);
    for(int i=0; i<grid.ncoarse; ++i) c_buf[i] = grid.flag1[i+1];
    write_record(c_buf);
    for(int i=0; i<grid.ncoarse; ++i) c_buf[i] = grid.cpu_map[i+1];
    write_record(c_buf);
}

} // namespace ramses
