#ifndef EXAMPLES_LECTUREDEMOS_CONVERGENCESTUDIES_FESPACELAGRANGEON_H_
#define EXAMPLES_LECTUREDEMOS_CONVERGENCESTUDIES_FESPACELAGRANGEON_H_

#include "felagrangeonquad.h"
#include "felagrangeontria.h"
#include "felagrangeonsegment.h"
#include <lf/uscalfe/uscalfe.h>
#include <lf/mesh/mesh.h>
#include <memory>



template<typename SCALAR>
class FeSpaceLagrangeON : public lf::uscalfe::UniformScalarFESpace<SCALAR> {
public:
    using Scalar = SCALAR;

    FeSpaceLagrangeON() = delete;
    FeSpaceLagrangeON(const FeSpaceLagrangeON&) = delete;
    FeSpaceLagrangeON(FeSpaceLagrangeON&&) noexcept = default;
    FeSpaceLagrangeON& operator=(const FeSpaceLagrangeON&) = delete;
    FeSpaceLagrangeON& operator=(FeSpaceLagrangeON&&) noexcept = default;
    explicit FeSpaceLagrangeON(const std::shared_ptr<const lf::mesh::Mesh> &mesh_p, unsigned N)
	: lf::uscalfe::UniformScalarFESpace<SCALAR>(
		mesh_p,
		std::make_shared<FeLagrangeONTria<SCALAR>>(N),
		std::make_shared<FeLagrangeONQuad<SCALAR>>(N),
		std::make_shared<FeLagrangeONSegment<SCALAR>>(N),
		std::make_shared<lf::uscalfe::FeLagrangePoint<SCALAR>>(N)) { }
    ~FeSpaceLagrangeON() override = default;
};



#endif // EXAMPLES_LECTUREDEMOS_CONVERGENCESTUDIES_FESPACELAGRANGEON_H_
