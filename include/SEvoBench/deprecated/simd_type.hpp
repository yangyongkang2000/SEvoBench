#pragma once
#if __arm64
#include "sse2neon/sse2neon.h"

// limit to 128byte, since we want to use ARM-neon
#define MAX_VECTOR_SIZE 128

// limit to sse4.2, sse2neon does not have any AVX instructions ( so far )
#define INSTRSET 6

// define unknown function
#define _mm_getcsr() 1

// simulate header included
#define __X86INTRIN_H
#endif
// finally include vectorclass
#include "../version2/vectorclass.h"
#include <concepts>
#include <tuple>
#include <type_traits>

namespace sevobench {
constexpr auto simd_id() {
#ifdef INSTRSET
  if constexpr (INSTRSET < 6)
    return -1;
  if constexpr (INSTRSET == 6)
    return 0;
  if constexpr (INSTRSET < 9) {
    return 1;
  } else {
    return 2;
  }
#else
  return -1;
#endif
}

} // namespace sevobench

namespace sevobench::simd_type_detail {

template <std::floating_point T>
  requires(!std::same_as<T, long double>) && (simd_id() >= 0)
using simd_type =
    std::conditional_t<std::is_same_v<T, float>,
                       std::tuple_element_t<simd_id() < 0 ? 0 : simd_id(),
                                            std::tuple<Vec4f, Vec8f, Vec16f>>,
                       std::tuple_element_t<simd_id() < 0 ? 0 : simd_id(),
                                            std::tuple<Vec2d, Vec4d, Vec8d>>>;

template <std::floating_point T>
  requires(!std::same_as<T, long double>) && (simd_id() >= 0)
constexpr int simd_width = std::is_same_v<T, float> ? (1 << (simd_id() + 2))
                                                    : (1 << (simd_id() + 1));

template <typename T, int Dim>
concept simd_dim = (simd_id() < 0 && Dim >= 1) ||
                   (simd_id() > 0 && Dim >= 1 &&
                    Dim % (simd_type_detail::simd_width<T>) == 0 &&
                    (!std::is_same_v<T, long double>));

} // namespace sevobench::simd_type_detail
