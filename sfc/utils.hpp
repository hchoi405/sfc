#ifndef HJ_UTILS
#define HJ_UTILS

#include <iostream>

namespace hj {
template <typename Arg>
void print(const std::string &name, const Arg &arg) {
    std::cout << name << ": " << arg << std::endl;
}

template <typename Arg>
void printArray(const std::string &name, const Arg &arg, const int n) {
    std::cout << name << "(" << n << "):[";
    for (int i = 0; i < n; ++i) {
        std::cout << arg[i];
        if (i < n-1) std::cout << ",";
    }
    std::cout << "]";
}

#define print(v) print(#v, v)
#define printArray(v, ...) printArray(#v, v, __VA_ARGS__)
}  // namespace hj

#endif