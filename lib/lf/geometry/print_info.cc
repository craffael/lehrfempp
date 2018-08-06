#include "print_info.h"


namespace lf::geometry {

CONTROLDECLARECOMMENT(Geometry, output_ctrl_, "output_ctrl_", "Diagnostics control for Geometry");

void PrintInfo(const Geometry &geom, std::ostream &o){

    lf::base::RefEl geom_refel = geom.RefEl();
    o << "Reference element type for geometry: " << geom_refel << std::endl;

    if(Geometry::output_ctrl_ > 10){
        // Dimensions and the derived type
        int dim_glob = geom.DimGlobal();
        int dim_local = geom.DimLocal();
        o << "Global dimension: " << dim_glob << std::endl;
        o << "Local dimension: " << dim_local << std::endl;
        o << "Derived type: " << typeid(geom).name() << std::endl;

        if(Geometry::output_ctrl_ > 90){
            // Geometry (coordinates of vertices)
            const Eigen::MatrixXd &ref_el_corners(geom_refel.NodeCoords());
            o << geom.Global(ref_el_corners) << std::endl;
        }
    }


} // void PrintInfo

std::ostream& operator<<(std::ostream& stream, const Geometry& geom) {

   if (Geometry::output_ctrl_ == 0){
       return stream << geom.RefEl();
   } else {
       PrintInfo(geom, stream);
   }
}


} // namespace lf::geometry
