#ifndef LBM_COLLISION_MODEL_HPP
#define LBM_COLLISION_MODEL_HPP

#include <lbm/multiple_relaxation_time_model.hpp>
#include <lbm/single_relaxation_time_model.hpp>
#include <variant>

namespace lbm {

using CollisionParameters = std::variant<SingleRelaxationTimeModelParameters,
                                         MultipleRelaxationTimeModelParameters>;
using CollisionModel =
    std::variant<SingleRelaxationTimeModel, MultipleRelaxationTimeModel>;

/**
 * @brief Create a CollisionModel object
 *
 * @param params Collision model parameters
 * @return CollisionModel
 * @throw std::runtime_error
 */
inline CollisionModel create_collision_model(
    const CollisionParameters& params) {
  if (std::holds_alternative<SingleRelaxationTimeModelParameters>(params)) {
    return SingleRelaxationTimeModel{
        std::get<SingleRelaxationTimeModelParameters>(params)};
  } else if (std::holds_alternative<MultipleRelaxationTimeModelParameters>(
                 params)) {
    return MultipleRelaxationTimeModel{
        std::get<MultipleRelaxationTimeModelParameters>(params)};
  } else {
    throw std::runtime_error(
        "Error: CollisionParameters holds neither of "
        "SingleRelaxationTimeModel or MultipleRelaxationTimeModel.");
  }
}

inline double get_relaxation_time(const CollisionParameters& params) {
  if (std::holds_alternative<SingleRelaxationTimeModelParameters>(params)) {
    return std::get<SingleRelaxationTimeModelParameters>(params).tau;
  } else if (std::holds_alternative<MultipleRelaxationTimeModelParameters>(
                 params)) {
    return std::get<MultipleRelaxationTimeModelParameters>(params).tau;
  } else {
    throw std::runtime_error(
        "Error: CollisionParameters holds neither of "
        "SingleRelaxationTimeModel or MultipleRelaxationTimeModel.");
  }
}

}  // namespace lbm

#endif  // LBM_COLLISION_MODEL_HPP