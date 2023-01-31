#ifndef ENCODINGS_ROOT_SRC_WRAPPERS_WRAPPER_UTILITIES_H_
#define ENCODINGS_ROOT_SRC_WRAPPERS_WRAPPER_UTILITIES_H_

#include <string>
#include <cstdlib>

namespace enc::lib {

int32_t cents_decimal_to_int(const std::string &value, bool as_cents);

uint32_t date_string_to_int(const std::string &value, bool add_jdc_offset);

} //

#endif//ENCODINGS_ROOT_SRC_WRAPPERS_WRAPPER_UTILITIES_H_
