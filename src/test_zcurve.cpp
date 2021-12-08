#include <chrono>
#include <iostream>
#include <random>

#include "sfc/Header.h"
#include "sfc/Vector.h"
#include "sfc/zcurve.hpp"

using namespace std;

const float eps = 1e-6;
using hclock = chrono::high_resolution_clock;
using duration = chrono::duration<double>;

template <size_t n>
std::array<float, n> create_random_data() {
    std::random_device r;
    std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
    std::mt19937 eng(seed);  // a source of random data

    std::uniform_real_distribution<float> dist;
    std::array<float, n> v;

    generate(begin(v), end(v), bind(dist, eng));
    return v;
}

bool dimensionTimeComparison(bool print) {
    int numData = 10000;

    {
        hj::Zcurve<3, float, uint64_t> zcurve3;
        using vec3 = hj::Vector<float, 3>;
        vec3 pmin(0.f), pmax(1.f);
        std::vector<vec3> arr1;
        for (int i = 0; i < numData; ++i) arr1.push_back(create_random_data<3>());

        auto start = hclock::now();
        zcurve3.order(vec3(0.f) - eps, vec3(1.f) + eps, arr1);
        duration time = hclock::now() - start;
        std::cout << "64 bits - 3D took " << time.count() << std::endl;
    }

    {
        hj::Zcurve<5, float, uint64_t> zcurve5;
        using vec5 = hj::Vector<float, 5>;
        vec5 pmin(0.f), pmax(1.f);
        std::vector<vec5> arr5;
        for (int i = 0; i < numData; ++i) arr5.push_back(create_random_data<5>());

        auto start = hclock::now();
        zcurve5.order(vec5(0.f) - eps, vec5(1.f) + eps, arr5);
        duration time = hclock::now() - start;
        std::cout << "64 bits - 5D took " << time.count() << std::endl;
    }

    {
        hj::Zcurve<3, float, uint128_t> zcurve3;
        using vec3 = hj::Vector<float, 3>;
        std::vector<vec3> arr2;
        for (int i = 0; i < numData; ++i) arr2.push_back(create_random_data<3>());

        auto start = hclock::now();
        zcurve3.order(vec3(0.f) - eps, vec3(1.f) + eps, arr2);
        duration time = hclock::now() - start;
        std::cout << "128 bits - 3D took " << time.count() << std::endl;
    }

    {
        hj::Zcurve<5, float, uint128_t> zcurve5;
        using vec5 = hj::Vector<float, 5>;
        std::vector<vec5> arr2;
        for (int i = 0; i < numData; ++i) arr2.push_back(create_random_data<5>());

        auto start = hclock::now();
        zcurve5.order(vec5(0.f) - eps, vec5(1.f) + eps, arr2);
        duration time = hclock::now() - start;
        std::cout << "128 bits - 5D took " << time.count() << std::endl;
    }

    {
        hj::Zcurve<18, float, uint128_t> zcurve18;
        using vec18 = hj::Vector<float, 18>;
        std::vector<vec18> arr2;
        for (int i = 0; i < numData; ++i) arr2.push_back(create_random_data<18>());

        auto start = hclock::now();
        zcurve18.order(vec18(0.f) - eps, vec18(1.f) + eps, arr2);
        duration time = hclock::now() - start;
        std::cout << "128 bits - 18D took " << time.count() << std::endl;
    }

    return false;
}

bool higherPrecisionTest(bool print) {
    const int DIM = 3;
    const int NUM_RANDOM_DATA = 10;
    using item = hj::Vector<float, DIM>;

    // Generate Data
    vector<item> arr;

    for (int i = 0; i < NUM_RANDOM_DATA; ++i) arr.push_back(create_random_data<DIM>());

    /*
            After ordering:
            0, [0, 0, 0]
            1, [0, 0.5, 0]
            2, [0.5, 0.5, 0]
            3, [0, 0, 0.5]
            4, [0.5, 0, 0.5]
            5, [0, 0.5, 0.5]
            6, [0.5, 0.5, 0.5]
    */

    /*arr.push_back({ .5f, .5f, .5f });
    arr.push_back({ 0.f, .5f, 0.f });
    arr.push_back({ 0.f, .5f, .5f });
    arr.push_back({ .5f, .5f, 0.f });
    arr.push_back({ 0.f, 0.f, 0.f });
    arr.push_back({ .5f, 0.f, .5f });
    arr.push_back({ 0.f, 0.f, .5f });*/

    uint64_t maxCells = uint64_t(1) << (64 / 3);
    hj::Vector<uint64_t, 3> floored;
    for (int i = 0; i < 3; ++i) {
        floored[i] = uint64_t(boost::multiprecision::cpp_dec_float_50(arr[0][i]) *
                              boost::multiprecision::cpp_dec_float_50(maxCells));
        std::cout << floored[i] << " ";
    }
    std::cout << std::endl;

    int numData = arr.size();
    ///////////////////////////////////////////////////////////////
    /* Using for loop ********************************************/
    ///////////////////////////////////////////////////////////////
    hj::Zcurve<DIM, float, uint64_t> zcurve;
    vector<pair<int, item>> arr1;

    int idx = 0;
    arr1.reserve(numData);
    for (size_t i = 0; i < numData; ++i) arr1.push_back({idx++, arr[i]});

    if (print) {
        cout << "\nBefore ordering: \n";
        for (auto el : arr1) {
            cout << el.second;
        }
    }

    // WARNING! Exact boundary points (e.g. (0,0) or (5,5)) cause wrong ordering
    // so make sure that the range is larger than the actual values
    auto start = hclock::now();
    zcurve.order<pair<int, item>>(item(0.f) - eps, item(1.f) + eps, arr1,
                                  [](auto& el) -> item& { return el.second; });
    duration time = hclock::now() - start;
    if (print) {
        cout << "\nAfter ordering (by For Loop): \n";
        for (auto el : arr1) {
            cout << zcurve.getMortonKey(el.second) << "\t" << el.second;
        }
    }
    std::cout << "took " << time.count() << "s to make mortonkey and sort" << std::endl;

    ///////////////////////////////////////////////////////////////
    /* Using Magic Bits ******************************************/
    ///////////////////////////////////////////////////////////////
    vector<item> arr2;

    arr2.reserve(numData);
    for (size_t i = 0; i < numData; ++i) arr2.push_back(arr[i]);

    start = hclock::now();
    zcurve.order(item(0.f) - eps, item(1.f) + eps, arr2);
    time = hclock::now() - start;
    if (print) {
        cout << "\nAfter ordering (by Magic Bits): \n";
        for (auto el : arr2) {
            cout << zcurve.getMortonKey(el) << "\t" << el;
        }
    }
    std::cout << "took " << time.count() << "s to make mortonkey and sort" << std::endl;

    bool result = true;
    for (int i = 0; i < arr.size(); ++i) {
        if (arr1[i].second != arr2[i]) {
            result = false;
            std::cout << arr1[i].second << "\t" << zcurve.getMortonKey(arr1[i].second)
                      << "\n is different from " << arr2[i] << "\t" << zcurve.getMortonKey(arr2[i])
                      << std::endl;
            break;
        }
    }
    return result;
}

bool magicBitsTest(bool print) {
    /////////////////////////////////////////////////////////////
    /* TEST1. Magic bits *****************************************/
    /////////////////////////////////////////////////////////////
    const int DIM = 3;
    hj::Zcurve<DIM, float, uint64_t> zcurve0;
    const int totalBits = 32;
    const int fieldBits = 10;
    const int dimension = 3;
    const int numMagicBits = (int)log2(fieldBits) + 1;
    std::vector<uint64_t> magicBits(numMagicBits, 0);

    zcurve0.getMagicBits(&(magicBits[0]), totalBits, fieldBits, dimension);
    std::cout << "magic bits:\n";
    for (auto v : magicBits) {
        std::cout << v << std::endl;
    }
    return true;
}

int main() {
    // if (higherPrecisionTest(true)) {
    //     std::cout << "higherPrecisionTest passed." << std::endl;
    // } else {
    //     std::cout << "higherPrecisionTest failed." << std::endl;
    // }

    // magicBitsTest(true);

    dimensionTimeComparison(false);

    return 0;
}