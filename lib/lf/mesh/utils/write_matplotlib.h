#ifndef WRITE_MATPLOTLIB_H
#define WRITE_MATPLOTLIB_H

#include <lf/mesh/mesh.h>

#include <string>

namespace lf::mesh::utils {
    /**
     * @brief Write affine triangulation data to file in matplotlib format
     * @param mesh the mesh to be stored to file
     * @param filename name of output file: .csv appended unless present
     *
     *
     */
    void writeMatplotlib(const lf::mesh::Mesh &mesh, std::string filename);

} // namespace lf::mesh::utils

#endif /* WRITE_MATPLOTLIB_H */
