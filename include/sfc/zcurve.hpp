#pragma once

#include <boost/multiprecision/cpp_dec_float.hpp>

#include "sfc.h"

using namespace boost::multiprecision;

namespace sfc {
template <typename DataType = float, typename UInt = uint64_t, int Dims = 3, int Bits = 8>
class Zcurve : public SFC<DataType, UInt, Dims, Bits> {
    using Base = SFC<DataType, UInt, Dims, Bits>;

    // Using base's member variables
    using Base::NumBitsTotal;
    using Base::NumStrataPerAxis;

    // Own member variables
    using bigfloat =
        typename std::conditional<std::is_same<UInt, uint128_t>::value,
                                  boost::multiprecision::cpp_dec_float_50, float>::type;
    const float Eps = 1e-6f;
    const int NumMagicBits = (int)log2(Bits) + 1;
    const UInt One = 1;
    const UInt OneShiftedByFieldBits = One << Bits;
    const size_t OneShiftedByLogFieldBits = size_t((Dims - 1) * (One << int(log2(Bits))));

    std::vector<UInt> magicBits;

   public:
    Zcurve() {
        std::cout << "\t[Initialization of Z-Curve]\n"
                  << "\t\tUse " << NumBitsTotal << " bits for " << Dims << " dimension of data\n"
                  << "\t\tEach field uses: " << Bits << " bits" << std::endl;
        magicBits.resize(NumMagicBits);
        getMagicBits(&(magicBits[0]), NumBitsTotal, Bits, Dims);
    }

    UInt encode(const std::array<uint32_t, Dims>& x) const { return getMortonKey(x); }
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
    template <typename ArrayLike>
    void order(
        const ArrayLike& pmin, const ArrayLike& pmax, std::vector<ArrayLike>& arr,
        const std::function<ArrayLike&(ArrayLike&)>& accessor = [](ArrayLike& e) { return e; },
        const std::function<void(void)>& progressUpdate = 0) const {
        std::vector<std::pair<ArrayLike, UInt>> ordered(arr.size());

        for (size_t i = 0; i < arr.size(); ++i) {
            ArrayLike& v = accessor(arr[i]);
            auto p = normalize(v, pmin - Eps, pmax + Eps);
            UInt key = getMortonKey(p);
            ordered[i] = std::make_pair(v, key);
        }

        std::sort(
            ordered.begin(), ordered.end(),
            [&progressUpdate](std::pair<ArrayLike, UInt>& lhs, std::pair<ArrayLike, UInt>& rhs) {
                if (progressUpdate) progressUpdate();
                return lhs.second < rhs.second;
            });

        for (auto i = 0u; i < arr.size(); ++i) {
            accessor(arr[i]) = ordered[i].first;
        }
    }

    /**
     * @brief Get the Morton Key object
     *
     * @tparam ArrayLike
     * @param p
     * @param pMin
     * @param pMax
     * @return UInt
     */
    template <typename ArrayLike>
    inline UInt getMortonKey(const ArrayLike& p, const ArrayLike& pMin,
                             const ArrayLike& pMax) const {
        return getMortonKey(Base::normalize(p, pMin - Eps, pMax + Eps));
    }

    /**
     * @brief Get the Morton Key object
     *
     * @tparam ArrayLike
     * @param _p point in [0, NuMStrataPerAxis-1]
     * @return UInt Morton key
     */
    template <typename ArrayLike>
    inline UInt getMortonKey(const ArrayLike& _p) const {
        UInt result = 0;

        for (int i = 0; i < Dims; i++) {
            result |= (getMortonKey(_p[i]) << i);
        }

        return result;
    }

    inline UInt getMortonKey(const uint32_t& _v) const {
        static auto bit_and = std::bit_and<uint32_t>();

        auto x = bit_and(_v, OneShiftedByFieldBits - 1);  // take field bits
        size_t remainingBits = OneShiftedByLogFieldBits;
        int idx = 0;
        for (int i = 0; i < NumMagicBits; ++i) {
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
        while (leftMostBit > 1) {
            UInt c = b & (b << (leftMostBit / 2));
            b = b ^ c ^ (c << (leftMostBit / 2) * (dimension - 1));
            mBits[idx++] = b;
            leftMostBit /= 2;
        }
    }

    // For-loop basaed calculation with bugs (results are different from magic bits)
    template <typename ArrayLike>
    UInt getMortonKey_for(const ArrayLike& p, const ArrayLike& pMin, const ArrayLike& pMax) const {
        return getMortonKey_for(normalize(p, pMin - Eps, pMax + Eps));
    }

    template <typename ArrayLike>
    UInt getMortonKey_for(const ArrayLike& _p) const {
        UInt result = 0;

        ArrayLike floored;
        for (int i = 0; i < Dims; ++i) {
            floored[i] = UInt(bigfloat(_p[i]) * bigfloat(NumStrataPerAxis));
        }

        for (int i = 0; i < NumBitsTotal; i++) {
            for (int j = 0; j < Dims; ++j) {
                result = result | ((floored[j] & (1ULL << i)) << (i + j));
            }
        }

        return result;
    }
};

}  // namespace sfc