/**
 * @file
 * @brief Define the LambdaMeshDataSet class.
 * @author Raffael Casagrande
 * @date   2018-08-04 03:24:39
 * @copyright MIT License
 */

#ifndef __72c5049b6331454e8fbba79135b5319b
#define __72c5049b6331454e8fbba79135b5319b

#include "mesh_data_set.h"

namespace lf::mesh::utils {

/**
 * @brief Calls a provided Lambda function to retrieve the data of an entity
 * @tparam ValueLambda The type of the lambda function that provides the values.
 * @tparam DefinedOnPredicate The type of the predicate that defines on which
 *                            entities this MeshDataSet is defined.
 *
 * #### Sample usage:
 * @snippet lambda_mesh_data_set.cc usage
 *
 * @note It is of course essential, that the two lambdas return the same value
 * if they are called twice with the same argument (deterministic)!
 */
template <class ValueLambda, class DefinedOnPredicate>
class LambdaMeshDataSet
    : public MeshDataSet<std::remove_reference_t<decltype(
          std::declval<ValueLambda>()(std::declval<Entity>()))>> {
  using lambda_return_t =
      decltype(std::declval<ValueLambda>()(std::declval<Entity>()));
  using T = std::remove_reference_t<lambda_return_t>;
  using base_t = MeshDataSet<T>;

 public:
  // NOLINTNEXTLINE(readability-const-return-type)
  [[nodiscard]] const T operator()(const Entity& e) const override {
    LF_ASSERT_MSG(DefinedOn(e), "MeshDataSet not defined on this entity.");
    return value_lambda_(e);
  }

  [[nodiscard]] bool DefinedOn(const Entity& e) const override {
    return defined_on_lambda_(e);
  }

 private:
  ValueLambda value_lambda_;
  DefinedOnPredicate defined_on_lambda_;

  explicit LambdaMeshDataSet(ValueLambda vl, DefinedOnPredicate dol)
      : value_lambda_(std::move(vl)), defined_on_lambda_(std::move(dol)) {}

  // Friends:
  template <class ValueLambda_, class DefinedOnPredicate_>
  friend std::shared_ptr<LambdaMeshDataSet<ValueLambda_, DefinedOnPredicate_>>
  make_LambdaMeshDataSet(ValueLambda_ vl, DefinedOnPredicate_ dol);  // NOLINT
};

/**
 * @brief Construct LambdaMeshDataSet
 * @param vl This lambda function should accept a `const Entity&` and provide
 *           the value for this entity.
 * @param dol (optional) This predicate should accept a `const Entity&` and
 * return a `bool` that specifies whether this mesh data set attaches data
 * to this entity. If this argument is not specified, the MeshDataSet is defined
 * for all entities.
 */
template <class ValueLambda, class DefinedOnPredicate = base::PredicateTrue>
std::shared_ptr<LambdaMeshDataSet<ValueLambda, DefinedOnPredicate>>
make_LambdaMeshDataSet(ValueLambda vl,
                       DefinedOnPredicate dol = base::PredicateTrue{}) {
  using impl_t = LambdaMeshDataSet<ValueLambda, DefinedOnPredicate>;
  return std::shared_ptr<impl_t>(new impl_t(vl, dol));
}

}  // namespace lf::mesh::utils

#endif  // __72c5049b6331454e8fbba79135b5319b
