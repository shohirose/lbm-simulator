#include "lbm/collision_model.hpp"

namespace lbm {

CollisionModel create_collision_model(const CollisionParameters& params) {
  if (std::holds_alternative<SingleRelaxationTimeModelParameters>(params)) {
    return SingleRelaxationTimeModel{
        std::get<SingleRelaxationTimeModelParameters>(params)};
  } else if (std::holds_alternative<MultipleRelaxationTimeModelParameters>(
                 params)) {
    return MultipleRelaxationTimeModel{
        std::get<MultipleRelaxationTimeModelParameters>(params)};
  } else if (std::holds_alternative<CentralMomentModelParameters>(params)) {
    return CentralMomentModel(std::get<CentralMomentModelParameters>(params));
  } else {
    throw std::runtime_error(
        "Error: CollisionParameters holds neither of "
        "SingleRelaxationTimeModel, MultipleRelaxationTimeModel, or "
        "CentralMomentModel.");
  }
}

double get_relaxation_time(const CollisionParameters& params) {
  return std::visit(
      [](const auto& params) { return params.get_relaxation_time(); }, params);
}

}  // namespace lbm
