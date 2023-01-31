#include "wrapper_utilities.h"
#include <date/date.h>
#include <iostream>

namespace enc::lib {

int32_t cents_decimal_to_int(const std::string &value, bool as_cents) {
  const auto find_decimal = value.find('.');
  if (find_decimal != std::string::npos) {
    int integer_half = std::stoi(value.substr(0, find_decimal));
    if (as_cents) {
      int decimal_half = 0;
      const auto decimal_substr = value.substr(find_decimal + 1);
      if (decimal_substr.length() == 2) {
        decimal_half = std::stoi(decimal_substr);
      }
      if (decimal_substr.length() == 1) {
        decimal_half = std::stoi(decimal_substr) * 10;
      }
      return ((integer_half * 100) + decimal_half);
    } else {
      return integer_half;
    }
  } else {
    if (as_cents) {
      return std::stoi(value) * 100;
    } else {
      return std::stoi(value);
    }
  }
}

uint32_t date_string_to_int(const std::string &value, bool add_jdc_offset) {
  tm time_struct;
  if (strptime(value.c_str(), "%Y-%m-%d", &time_struct)) {
    int year_conv = time_struct.tm_year + 1900;
    auto ymd =
        date::year(year_conv) / (time_struct.tm_mon + 1) / time_struct.tm_mday;
    auto out = (unsigned int) (date::sys_days{ymd}.time_since_epoch().count());
    if (add_jdc_offset) {
      return out + 2440587;
    } else {
      return out;
    }
  } else {
    return 0;
  }
}

}// namespace enc::lib