#ifndef THESIS_POST_PROCESSING_NORMS_H
#define THESIS_POST_PROCESSING_NORMS_H

/**
 * @file norms.h
 * @brief Functions to compute norms of discontinuous vector valued functions
 * over a mesh
 */

#include <lf/mesh/entity.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>
#include <lf/quad/quad.h>
#include <lf/uscalfe/uscalfe.h>
#include <utils.h>

#include <functional>

namespace projects::ipdg_stokes {

namespace post_processing {

/**
 * @brief Compute the @f$L^2@f$-norm of a vector valued function over a mesh
 * @tparam MF The type of mesh function to take the norm from
 * @tparam QR_SELECTOR The type of the `qr_selector` parameter
 * @param mesh The mesh over which the function is defined
 * @param f The function to take the norm from
 * @param qr_selector A function giving the quadrature rule to use for a given
 * element
 * @returns The @f$L^2@f$-norm of the function over the mesh
 *
 * This function computes the value of the following formula
 * @f[
 *  ||f||_0^2 = \sum_{K\in\mathcal{M}} \int_K\! ||f(x)||^2 \,\mathrm{d}x
 * @f]
 */
template <typename MF, typename QR_SELECTOR>
double L2norm(const std::shared_ptr<const lf::mesh::Mesh> &mesh, const MF &f,
              const QR_SELECTOR &qr_selector) {
  return std::sqrt(lf::uscalfe::IntegrateMeshFunction(
      *mesh, lf::uscalfe::squaredNorm(f), qr_selector));
}

/**
 * @brief Compute the DG-norm of a vector valued function over a mesh
 * @tparam MF_F The type of mesh function to take the norm from
 * @tparam MF_GRAD The type of the gradient of the mesh function to take the
 * norm from
 * @tparam QR_SELECTOR The type of the `qr_selector` parameter
 * @param mesh The mesh over which the function is defined
 * @param f The function to take the norm from
 * @param f_grad The gradient of the function to take the norm from
 * @param qr_selector A function giving the quadrature rule to use for a given
 * element
 * @returns The DG-norm of the function over the mesh
 *
 * This functoin computes the following formula
 * @f[
 *  ||f||_{1,h}^2 = \sum_{K\in\mathcal{M}} \int_K\! ||\nabla f(x)||^2
 * \,\mathrm{d}x + \sum_{e\in\mathcal{E}(\mathcal{M})} \int_e\!
 * \frac{1}{|e|}||\left[\!\left[ f(x) \otimes n \right]\!\right]||^2
 * \,\mathrm{d}x
 * @f]
 */
template <typename MF_F, typename MF_GRAD, typename QR_SELECTOR>
double DGnorm(const std::shared_ptr<const lf::mesh::Mesh> &mesh, const MF_F &f,
              const MF_GRAD &f_grad, const QR_SELECTOR &qr_selector) {
  double norm2 = 0;
  // Compute the DG norm on the elements
  norm2 += lf::uscalfe::IntegrateMeshFunction(
      *mesh, lf::uscalfe::squaredNorm(f_grad), qr_selector);
  // Compute the norm on the edges
  const auto boundary = lf::mesh::utils::flagEntitiesOnBoundary(mesh, 1);
  using jump_t = std::remove_cv_t<decltype(
      (std::declval<lf::mesh::utils::MeshFunctionReturnType<MF_F>>() *
       std::declval<Eigen::Vector2d>().transpose())
          .eval())>;
  lf::mesh::utils::CodimMeshDataSet<std::vector<jump_t>> eval(mesh, 1);
  for (const auto entity : mesh->Entities(0)) {
    const auto edges = entity->SubEntities(1);
    const auto orientations = entity->RelativeOrientations();
    const auto normals =
        projects::ipdg_stokes::mesh::computeOutwardNormals(*entity);
    std::array<lf::quad::QuadRule, 3> quadrules;
    std::array<Eigen::MatrixXd, 3> quadpoints;
    for (int i = 0; i < 3; ++i) {
      quadrules[i] = qr_selector(*edges[i]);
      quadpoints[i] = Eigen::MatrixXd(2, quadrules[i].NumPoints());
    }
    quadpoints[0] << quadrules[0].Points(),
        Eigen::VectorXd::Zero(quadrules[0].NumPoints());
    quadpoints[1] << Eigen::VectorXd::Ones(quadrules[1].NumPoints()) -
                         quadrules[1].Points(),
        quadrules[1].Points();
    quadpoints[2] << Eigen::VectorXd::Zero(quadrules[2].NumPoints()),
        Eigen::VectorXd::Ones(quadrules[2].NumPoints()) - quadrules[2].Points();
    for (int i = 0; i < 3; ++i) {
      const auto edge = edges[i];
      const auto segment_eval = f(*entity, quadpoints[i]);
      if (orientations[i] == lf::mesh::Orientation::positive) {
        if (eval(*edge).size() == 0) {
          eval(*edge).resize(segment_eval.size());
          for (int j = 0; j < segment_eval.size(); ++j) {
            eval(*edge)[j] = segment_eval[j] * normals.col(i).transpose();
          }
        } else {
          for (int j = 0; j < segment_eval.size(); ++j) {
            eval(*edge)[j] += segment_eval[j] * normals.col(i).transpose();
          }
        }
      } else {
        if (eval(*edge).size() == 0) {
          eval(*edge).resize(segment_eval.size());
          for (int j = 0; j < segment_eval.size(); ++j) {
            eval(*edge)[j] = segment_eval[segment_eval.size() - j - 1] *
                             normals.col(i).transpose();
          }
        } else {
          for (int j = 0; j < segment_eval.size(); ++j) {
            eval(*edge)[j] += segment_eval[segment_eval.size() - j - 1] *
                              normals.col(i).transpose();
          }
        }
      }
    }
  }
  for (const auto edge : mesh->Entities(1)) {
    if (!boundary(*edge)) {
      const auto qr = qr_selector(*edge);
      const auto geom = edge->Geometry();
      const auto ie = geom->IntegrationElement(qr.Points());
      const Eigen::VectorXd weights =
          (qr.Weights().array() * ie.array()).matrix() /
          lf::geometry::Volume(*geom);
      const auto &p = eval(*edge);
      for (int i = 0; i < p.size(); ++i) {
        norm2 += weights[i] * p[i].squaredNorm();
      }
    }
  }
  return std::sqrt(norm2);
}

}  // end namespace post_processing

}  // end namespace projects::ipdg_stokes

#endif  // THESIS_POST_PROCESSING_NORMS_H
