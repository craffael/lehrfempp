#ifndef EXAMPLES_LECTUREDEMOS_CONVERGENCESTUDIES_FELAGRANGEONSEGMENT_H
#define EXAMPLES_LECTUREDEMOS_CONVERGENCESTUDIES_FELAGRANGEONSEGMENT_H


#include <lf/uscalfe/uscalfe.h>
#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <cmath>


template<typename SCALAR>
class FeLagrangeONSegment final : public lf::uscalfe::ScalarReferenceFiniteElement<SCALAR> {
public:
    FeLagrangeONSegment(const FeLagrangeONSegment&) = default;
    FeLagrangeONSegment(FeLagrangeONSegment&&) = default;
    FeLagrangeONSegment& operator=(const FeLagrangeONSegment&) = default;
    FeLagrangeONSegment& operator=(FeLagrangeONSegment&&) = default;
    ~FeLagrangeONSegment() = default;

    FeLagrangeONSegment(unsigned degree) : degree_(degree), ref_func_coeffs_() {
	eval_nodes_ = ComputeEvaluationNodes(degree);
	ref_func_coeffs_ = ComputePolyBasis(eval_nodes_).inverse().transpose();
    }

    [[nodiscard]] lf::base::RefEl RefEl() const override {
	return lf::base::RefEl::kSegment();
    }

    [[nodiscard]] unsigned Degree() const override {
	return degree_;
    }

    [[nodiscard]] lf::base::size_type NumRefShapeFunctions() const override {
	return degree_ + 1;
    }

    [[nodiscard]] lf::base::size_type NumRefShapeFunctions(lf::assemble::dim_t codim) const override {
	return codim == 0 ? degree_-1 : 1;
    }

    [[nodiscard]] lf::base::size_type NumRefShapeFunctions(lf::assemble::dim_t codim, lf::base::sub_idx_t /*subidx*/) const override {
	return NumRefShapeFunctions(codim);
    }

    [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> EvalReferenceShapeFunctions(const Eigen::MatrixXd &refcoords) const override {
	const auto poly_basis = ComputePolyBasis(refcoords);
	return ref_func_coeffs_ * poly_basis.transpose();
    }

    [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> GradientsReferenceShapeFunctions(const Eigen::MatrixXd &refcoords) const override {
	const auto dx = ComputePolyBasisDerivative(refcoords);
	return ref_func_coeffs_ * dx.transpose();
    }

    [[nodiscard]] Eigen::MatrixXd EvaluationNodes() const override {
	return eval_nodes_;
    }

    [[nodiscard]] lf::base::size_type NumEvaluationNodes() const override {
	return NumRefShapeFunctions();
    }


private:
    unsigned degree_;
    Eigen::MatrixXd eval_nodes_;
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> ref_func_coeffs_;

    Eigen::MatrixXd ComputeEvaluationNodes(unsigned p) const {
	Eigen::MatrixXd nodes(1, p+1);
	// Add the evaluation nodes corresponding to the vertices of the segment
	nodes(0, 0) = 0;
	nodes(0, 1) = 1;
	// Add the evaluation nodes corresponding to the interior of the segment
	for (int i = 0 ; i < p-1 ; ++i) {
	    nodes(0, i+2) = static_cast<double>(i+1) / p;
	}
	return nodes;
    }

    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> ComputePolyBasis(const Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> &refcoords) const {
	Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> poly_basis(refcoords.cols(), degree_+1);
	for (unsigned exponent = 0 ; exponent <= degree_ ; ++exponent) {
	    for (unsigned node_idx = 0 ; node_idx < refcoords.cols() ; ++node_idx) {
		poly_basis(node_idx, exponent) = std::pow(refcoords(0, node_idx), exponent);
	    }
	}
	return poly_basis;
    }

    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> ComputePolyBasisDerivative(const Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> &refcoords) const {
	Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> dx(refcoords.cols(), degree_+1);
	dx.col(0).setZero();
	for (unsigned exponent = 1 ; exponent <= degree_ ; ++exponent) {
	    for (unsigned node_idx = 0 ; node_idx < refcoords.cols() ; ++node_idx) {
		dx(node_idx, exponent) = exponent * std::pow(refcoords(0, node_idx), exponent-1);
	    }
	}
	return dx;
    }
};


#endif // EXAMPLES_LECTUREDEMOS_CONVERGENCESTUDIES_FELAGRANGEONSEGMENT_H
