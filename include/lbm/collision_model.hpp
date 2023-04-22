#ifndef LBM_COLLISION_MODEL_HPP
#define LBM_COLLISION_MODEL_HPP

#include <lbm/central_moment_model.hpp>
#include <lbm/multiple_relaxation_time_model.hpp>
#include <lbm/single_relaxation_time_model.hpp>
#include <variant>

namespace lbm {

using CollisionParameters = std::variant<SingleRelaxationTimeModelParameters,
                                         MultipleRelaxationTimeModelParameters,
                                         CentralMomentModelParameters>;
using CollisionModel =
    std::variant<SingleRelaxationTimeModel, MultipleRelaxationTimeModel,
                 CentralMomentModel>;

/**
 * @brief Create a CollisionModel object
 *
 * @param params Collision model parameters
 * @return CollisionModel
 * @throw std::runtime_error
 */
CollisionModel create_collision_model(const CollisionParameters& params);

/**
 * @brief Get the relaxation time for dynamic viscosity
 *
 * @param params Collision model parameters
 * @return double Relaxation time
 */
double get_relaxation_time(const CollisionParameters& params);

}  // namespace lbm

#endif  // LBM_COLLISION_MODEL_HPP