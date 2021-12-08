#pragma once

#include <boost/multiprecision/cpp_dec_float.hpp>

#include "Vector.h"
#include "sfc.h"

using namespace boost::multiprecision;

namespace sfc {
template <typename DataType = float, typename UInt = uint64_t, int Dims = 3, int Bits = 8>
class Zcurve : public SFC<DataType, UInt, Dims, Bits> {
    using Base = SFC<DataType, UInt, Dims, Bits>;

    // Using base's types
    using typename Base::_Point;

    // Using base's member variables
    using Base::numBitsTotal;
    using Base::numStrataPerAxis;

    // Own member variables
    using bigfloat =
        typename std::conditional<std::is_same<UInt, uint128_t>::value,
                                  boost::multiprecision::cpp_dec_float_50, float>::type;
    const float eps = 1e-6f;
    const int numMagicBits = (int)log2(Bits) + 1;
    const UInt One = 1;
    const UInt oneShiftedByFieldBits = One << Bits;
    const size_t oneShiftedByLogFieldBits = size_t((Dims - 1) * (One << int(log2(Bits))));

    std::vector<UInt> magicBits;

   public:
    Zcurve() {
        std::cout << "\t[Initialization of Z-Curve]\n"
                  << "\t\tUse " << numBitsTotal << " bits for " << Dims << " dimension of data\n"
                  << "\t\tEach field uses: " << Bits << " bits" << std::endl;
        magicBits.resize(numMagicBits);
        getMagicBits(&(magicBits[0]), numBitsTotal, Bits, Dims);
    }

   public:
    UInt encode(const std::array<uint32_t, Dims>& x) const {
        UInt result = 0;
        for (int i = 0; i < Dims; i++) {
            assert(x[i] >= 0 && x[i] < numStrataPerAxis);
            result |= (getMortonKey(x[i]) << i);
        }
        return result;
    }
    UInt encode(const std::array<DataType, Dims>& x, const std::array<DataType, Dims>& pMin,
                const std::array<DataType, Dims>& pMax) const {
        std::array<uint32_t, Dims> uarr;
        for (int i = 0; i < Dims; ++i) uarr[i] = Base::normalize(x[i], pMin[i], pMax[i]);

        return encode(uarr);
    }
    void decode(const UInt v, std::array<uint32_t, Dims>& x) const { sfc::NotImplemented(); }
    void decode(const UInt v, std::array<DataType, Dims>& x, const std::array<DataType, Dims>& pMin,
                const std::array<DataType, Dims>& pMax) const {
        sfc::NotImplemented();
    }

   public:
    // For arbtrary type of element with accessor
    template <typename ArrayElementType>
    void order(const _Point& pmin, const _Point& pmax, std::vector<ArrayElementType>& arr,
               const std::function<_Point&(ArrayElementType&)>& accessor = 0,
               const std::function<void(void)>& progressUpdate = 0) const {
        std::vector<std::pair<_Point, UInt>> ordered(arr.size());

        for (size_t i = 0; i < arr.size(); ++i) {
            _Point& v = accessor(arr[i]);
            auto p = normalize(v, pmin - eps, pmax + eps);
            UInt key = getMortonKey(p);
            ordered[i] = std::make_pair(v, key);
        }

        std::sort(ordered.begin(), ordered.end(),
                  [&progressUpdate](std::pair<_Point, UInt>& lhs, std::pair<_Point, UInt>& rhs) {
                      if (progressUpdate) progressUpdate();
                      return lhs.second < rhs.second;
                  });

        for (auto i = 0u; i < arr.size(); ++i) {
            accessor(arr[i]) = ordered[i].first;
        }
    }

    // Overloading for sfc::Vector
    void order(const _Point& pmin, const _Point& pmax, std::vector<_Point>& arr,
               const std::function<_Point&(_Point&)>& accessor = 0,
               const std::function<void(void)>& progressUpdate = 0) const {
        std::vector<std::pair<_Point, UInt>> ordered(arr.size());

        for (size_t i = 0; i < arr.size(); ++i) {
            _Point& v = arr[i];
            auto p = normalize(v, pmin - eps, pmax + eps);
            UInt key = getMortonKey(p);
            ordered[i] = std::make_pair(v, key);
        }

        std::sort(ordered.begin(), ordered.end(),
                  [&progressUpdate](std::pair<_Point, UInt>& lhs, std::pair<_Point, UInt>& rhs) {
                      if (progressUpdate) progressUpdate();
                      return lhs.second < rhs.second;
                  });

        for (auto i = 0u; i < arr.size(); ++i) {
            arr[i] = ordered[i].first;
        }
    }

    // Overloading for std::array
    void order(const _Point& pmin, const _Point& pmax,
               std::vector<std::array<DataType, Dims>>& arr) const {
        std::vector<std::pair<_Point, UInt>> ordered(arr.size());

        for (size_t i = 0; i < arr.size(); ++i) {
            _Point v(arr[i]);
            auto p = Base::normalize(v, pmin - eps, pmax + eps);
            UInt key = getMortonKey(p);
            ordered[i] = std::make_pair(v, key);
        }

        std::sort(ordered.begin(), ordered.end(),
                  [](std::pair<_Point, UInt>& lhs, std::pair<_Point, UInt>& rhs) {
                      return lhs.second < rhs.second;
                  });

        for (auto i = 0u; i < arr.size(); ++i) {
            arr[i] = ordered[i].first;
        }
    }

    inline UInt getMortonKey(const _Point& p, const _Point& pMin, const _Point& pMax) const {
        return getMortonKey(normalize(p, pMin - eps, pMax + eps));
    }

    // input: _p normalized point in [0, 1]
    inline UInt getMortonKey(const _Point& _p) const {
        UInt result = 0;

        for (int i = 0; i < Dims; i++) {
            result |= (getMortonKey(UInt(bigfloat(_p[i]) * bigfloat(numStrataPerAxis))) << i);
        }

        return result;
    }

    inline UInt getMortonKey(const UInt& _v) const {
        static auto bit_and = std::bit_and<UInt>();

        auto x = bit_and(_v, oneShiftedByFieldBits - 1);  // take field bits
        size_t remainingBits = oneShiftedByLogFieldBits;
        int idx = 0;
        for (int i = 0; i < numMagicBits; ++i) {
            x = (x | (x << remainingBits)) & magicBits[idx++];
            remainingBits /= 2;
        }

        return x;
    }

    void getMagicBits(UInt* mBits, int totalBits, int bitsPerAxis, int dimension) {
        assert(bitsPerAxis <
               sizeof(size_t) * 8);  // sizeof(size_t)*8 = 64 is the maximum size of integer
                                     // that can be used as an operand for shift operator

        const UInt one = 1;

        // 8 bits
        size_t leftMostBit = size_t(one << size_t(log2(bitsPerAxis)));
        // 8 bits with all 1
        UInt a = (one << leftMostBit) - 1;

        UInt b = a | (a << (leftMostBit * dimension));
        int idx = 0;
        mBits[idx++] = b;
        // std::cout << b << std::endl;
        while (leftMostBit > 1) {
            UInt c = b & (b << (leftMostBit / 2));
            b = b ^ c ^ (c << (leftMostBit / 2) * (dimension - 1));
            mBits[idx++] = b;
            // std::cout << b << std::endl;
            leftMostBit /= 2;
        }
    }

    // For-loop basaed calculation with bugs (results are different from magic bits)
    UInt getMortonKey_for(const _Point& p, const _Point& pMin, const _Point& pMax) const {
        return getMortonKey_for(normalize(p, pMin - eps, pMax + eps));
    }

    UInt getMortonKey_for(const _Point& _p) const {
        UInt result = 0;

        Vector<UInt, Dims> floored;
        for (int i = 0; i < Dims; ++i) {
            floored[i] = UInt(bigfloat(_p[i]) * bigfloat(numStrataPerAxis));
        }

        for (int i = 0; i < numBitsTotal; i++) {
            for (int j = 0; j < Dims; ++j) {
                result = result | ((floored[j] & (1ULL << i)) << (i + j));
            }
        }

        return result;
    }
};

}  // namespace sfc