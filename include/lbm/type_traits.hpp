#ifndef LBM_TYPE_TRAITS_HPP
#define LBM_TYPE_TRAITS_HPP

namespace lbm {

class SingleRelaxationTimeModel;
struct SingleRelaxationTimeModelParameters;
class MultipleRelaxationTimeModel;
struct MultipleRelaxationTimeModelParameters;

template <typename T>
struct get_collision_parameters {};

template <>
struct get_collision_parameters<SingleRelaxationTimeModel> {
  using type = SingleRelaxationTimeModelParameters;
};

template <>
struct get_collision_parameters<MultipleRelaxationTimeModel> {
  using type = MultipleRelaxationTimeModelParameters;
};

template <typename T>
using get_collision_parameters_t = typename get_collision_parameters<T>::type;

}  // namespace lbm

#endif  // LBM_TYPE_TRAITS_HPP