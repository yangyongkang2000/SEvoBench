#pragma once
#include <cmath>
#include <iterator>
#include <type_traits>
#include <utility>

template <typename InputIt>
auto mean_std(InputIt beg, InputIt end) -> std::pair<
    std::conditional_t<
        std::is_integral_v<typename std::iterator_traits<InputIt>::value_type>,
        double, typename std::iterator_traits<InputIt>::value_type>,
    std::conditional_t<
        std::is_integral_v<typename std::iterator_traits<InputIt>::value_type>,
        double, typename std::iterator_traits<InputIt>::value_type>> {
  using value_type = typename std::iterator_traits<InputIt>::value_type;
  using promoted_type =
      std::conditional_t<std::is_integral_v<value_type>, double, value_type>;

  promoted_type mean = 0;
  promoted_type M2 = 0;
  size_t count = 0;

  for (auto it = beg; it != end; ++it) {
    const auto x = static_cast<promoted_type>(*it);
    ++count;

    const auto delta = x - mean;
    mean += delta / count;
    M2 += delta * (x - mean);
  }

  return count < 2 ? std::pair{promoted_type{0}, promoted_type{0}}
                   : std::pair{mean, std::sqrt(M2 / (count - 1))};
}