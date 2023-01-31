#ifndef ENCODINGS_SRC_ENCODINGS_H_
#define ENCODINGS_SRC_ENCODINGS_H_

#include <cstdlib>
#include <stdexcept>

namespace enc::lib {

enum Encoding {
  NONE,DFE,EDFE
};

template<class T>
static inline consteval T upper_size_helper() {
  if constexpr (sizeof(T) == 8) {
    return (T) 6;
  } else if constexpr (sizeof(T) == 4) {
    return (T) 5;
  } else if constexpr (sizeof(T) == 2) {
    return (T) 4;
  } else if constexpr (sizeof(T) == 1) {
    return (T) 3;
  } else {
    return (T) 0;
  }
}

template<class T>
static inline T ToDFE(const T in) {
  constexpr T B = sizeof(T) * 8;
  constexpr T L_U = upper_size_helper<T>();
  constexpr T L_L = B - L_U;
  constexpr T mask_L = (((T) 1) << L_L) - 1;
  if (in == (T) 0) [[unlikely]] {
    return (T) 0;
  }
  T clz = ((T) __builtin_clzll((unsigned long long) in)) - (64 - B);
  if (clz < (L_U - 1)) [[unlikely]] {
    throw std::out_of_range("ToDFE: CLZ < (L_U - 1)");
  }
  const T upper = (B - clz) << L_L;
  const T lower = in << (clz - L_U + 1);
  return (upper | (lower & mask_L));
}

template<class T>
static T FromDFE(const T in) {
  constexpr T B = sizeof(T) * 8;
  constexpr T L_U = upper_size_helper<T>();
  constexpr T L_L = B - L_U;
  constexpr T mask_L = (((T) 1) << L_L) - 1;
  if (in == (T) 0) [[unlikely]] {
    return (T) 0;
  }
  const T val = in & mask_L | (mask_L + 1);
  const T shift = (L_L + 1) - (in >> L_L);
  return val >> shift;
}

template<class T>
static inline T ToEDFE(const T in) {
  constexpr T B = sizeof(T) * 8;
  constexpr T L_U = upper_size_helper<T>();
  constexpr T L_L = B - L_U;
  constexpr T mask_L = (((T) 1) << (L_L - 2)) - 1;
  constexpr T format = ((T) 1) << (B - 2);
  if (in == (T) 0) [[unlikely]] {
    return (T) 0;
  }
  const T in_abs = std::abs(in);
  T clz = ((T) __builtin_clzll((unsigned long long) in_abs)) - (64 - B);
  if (clz <= 1) [[unlikely]] {
    throw std::out_of_range("ToEDFE: CLZ <= 1");
  } else if (clz <= L_U) {
    return in ^ format;
  } else {
    const T upper = (B - clz) << (L_L - 2);
    const T lower = in_abs << (clz - L_U - 1);
    const T out = upper | (lower & mask_L);
    return (in < 0) ? (~out) : (out);
  }
}

template<class T>
static T FromEDFE(const T in) {
  constexpr T B = sizeof(T) * 8;
  constexpr T L_U = upper_size_helper<T>();
  constexpr T L_L = B - L_U;
  constexpr T mask_L = (((T) 1) << (L_L - 2)) - 1;
  constexpr T format = ((T) 1) << (B - 2);
  if (in == (T) 0) [[unlikely]] {
    return (T) 0;
  }
  const T in_flip = (in < 0) ? (~in) : (in);
  if (in_flip & format) {
    return in ^ format;
  } else {
    const T val = (in_flip & mask_L) | (mask_L + 1);
    const T shift = (L_L - 1) - (in_flip >> (L_L - 2));
    const T out = val >> shift;
    return (in < 0) ? (-out) : (out);
  }
}

}// namespace enc::lib

#endif//ENCODINGS_SRC_ENCODINGS_H_
