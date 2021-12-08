#pragma once

#include "Header.h"
#include "Vector.h"

namespace hj {
template <int Dimension, typename DataType, typename UnsignedIntegerType = uint64_t>
class Zcurve {
    static_assert(Dimension > 0, "Parameter 'Dimension' must be > 0.");
    static_assert(std::is_constructible<UnsignedIntegerType, DataType>::value,
                  "UnsignedIntegerType instance cannot be constructed using DataType.");

    using _Vector = Vector<DataType, Dimension>;
    using bigfloat =
        typename std::conditional<std::is_same<UnsignedIntegerType, uint128_t>::value,
                                  boost::multiprecision::cpp_dec_float_50, float>::type;

    // clang-format off
	const int numBits =
		(std::is_same<UnsignedIntegerType, uint8_t>::value) ? 8 :
		(std::is_same<UnsignedIntegerType, uint16_t>::value) ? 16 :
		(std::is_same<UnsignedIntegerType, uint32_t>::value) ? 32 :
		(std::is_same<UnsignedIntegerType, uint64_t>::value) ? 64 :
		(std::is_same<UnsignedIntegerType, uint128_t>::value) ? 128 : 0;
    // clang-format on

    const int halfNumBits = numBits / 2;
    const float eps = 1e-6f;
    const UnsignedIntegerType maxCells = UnsignedIntegerType(1) << (numBits / Dimension);
    const int fieldBits =
        numBits / Dimension;  // use possible maximum number of bits for each field
    const int numMagicBits = (int)log2(fieldBits) + 1;
    const UnsignedIntegerType One = 1;
    const UnsignedIntegerType oneShiftedByFieldBits = One << fieldBits;
    const size_t oneShiftedByLogFieldBits = size_t((Dimension - 1) * (One << int(log2(fieldBits))));

    std::vector<UnsignedIntegerType> magicBits;

   public:
    Zcurve() {
        std::cout << "\t[Initialization of Z-Curve]\n"
                  << "Use " << numBits << " bits for " << Dimension << " dimension of data\n"
                  << "Each field uses: " << fieldBits << " bits" << std::endl;
        magicBits.resize(numMagicBits);
        getMagicBits(&(magicBits[0]), numBits, fieldBits, Dimension);
    }

    // For arbtrary type of element with accessor
    template <typename ArrayElementType>
    void order(const _Vector& pmin, const _Vector& pmax, std::vector<ArrayElementType>& arr,
               const std::function<_Vector&(ArrayElementType&)>& accessor = 0,
               const std::function<void(void)>& progressUpdate = 0) const {
        std::vector<std::pair<_Vector, UnsignedIntegerType>> ordered(arr.size());

        for (size_t i = 0; i < arr.size(); ++i) {
            _Vector& v = accessor(arr[i]);
            auto p = normalize(v, pmin - eps, pmax + eps);
            UnsignedIntegerType key = getMortonKey(p);
            ordered[i] = std::make_pair(v, key);
        }

        std::sort(ordered.begin(), ordered.end(),
                  [&progressUpdate](std::pair<_Vector, UnsignedIntegerType>& lhs,
                                    std::pair<_Vector, UnsignedIntegerType>& rhs) {
                      if (progressUpdate) progressUpdate();
                      return lhs.second < rhs.second;
                  });

        for (auto i = 0u; i < arr.size(); ++i) {
            accessor(arr[i]) = ordered[i].first;
        }
    }

    // Overloading for hj::Vector
    void order(const _Vector& pmin, const _Vector& pmax, std::vector<_Vector>& arr,
               const std::function<_Vector&(_Vector&)>& accessor = 0,
               const std::function<void(void)>& progressUpdate = 0) const {
        std::vector<std::pair<_Vector, UnsignedIntegerType>> ordered(arr.size());

        for (size_t i = 0; i < arr.size(); ++i) {
            _Vector& v = arr[i];
            auto p = normalize(v, pmin - eps, pmax + eps);
            UnsignedIntegerType key = getMortonKey(p);
            ordered[i] = std::make_pair(v, key);
        }

        std::sort(ordered.begin(), ordered.end(),
                  [&progressUpdate](std::pair<_Vector, UnsignedIntegerType>& lhs,
                                    std::pair<_Vector, UnsignedIntegerType>& rhs) {
                      if (progressUpdate) progressUpdate();
                      return lhs.second < rhs.second;
                  });

        for (auto i = 0u; i < arr.size(); ++i) {
            arr[i] = ordered[i].first;
        }
    }

    // Overloading for std::array
    void order(const _Vector& pmin, const _Vector& pmax,
               std::vector<std::array<DataType, Dimension>>& arr) const {
        std::vector<std::pair<_Vector, UnsignedIntegerType>> ordered(arr.size());

        for (size_t i = 0; i < arr.size(); ++i) {
            _Vector v(arr[i]);
            auto p = normalize(v, pmin - eps, pmax + eps);
            UnsignedIntegerType key = getMortonKey(p);
            ordered[i] = std::make_pair(v, key);
        }

        std::sort(
            ordered.begin(), ordered.end(),
            [](std::pair<_Vector, UnsignedIntegerType>& lhs,
               std::pair<_Vector, UnsignedIntegerType>& rhs) { return lhs.second < rhs.second; });

        for (auto i = 0u; i < arr.size(); ++i) {
            arr[i] = ordered[i].first;
        }
    }

    inline UnsignedIntegerType getMortonKey(const _Vector& p, const _Vector& pMin,
                                            const _Vector& pMax) const {
        return getMortonKey(normalize(p, pMin - eps, pMax + eps));
    }

    // input: _p normalized point in [0, 1]
    inline UnsignedIntegerType getMortonKey(const _Vector& _p) const {
        UnsignedIntegerType result = 0;

        for (int i = 0; i < Dimension; i++) {
            result |=
                (getMortonKey(UnsignedIntegerType(bigfloat(_p[i]) * bigfloat(maxCells))) << i);
        }

        return result;
    }

    inline UnsignedIntegerType getMortonKey(const UnsignedIntegerType& _v) const {
        static auto bit_and = std::bit_and<UnsignedIntegerType>();

        auto x = bit_and(_v, oneShiftedByFieldBits - 1);  // take field bits
        size_t remainingBits = oneShiftedByLogFieldBits;
        int idx = 0;
        for (int i = 0; i < numMagicBits; ++i) {
            x = (x | (x << remainingBits)) & magicBits[idx++];
            remainingBits /= 2;
        }

        return x;
    }

    void getMagicBits(UnsignedIntegerType* mBits, int totalBits, int fieldBits, int dimension) {
        assert(totalBits / dimension == fieldBits);
        assert(totalBits <= numBits);
        assert(fieldBits <
               sizeof(size_t) * 8);  // sizeof(size_t)*8 = 64 is the maximum size of integer
                                     // that can be used as an operand for shift operator

        const UnsignedIntegerType one = 1;

        // 8 bits
        size_t leftMostBit = size_t(one << size_t(log2(fieldBits)));
        // 8 bits with all 1
        UnsignedIntegerType a = (one << leftMostBit) - 1;

        UnsignedIntegerType b = a | (a << (leftMostBit * dimension));
        int idx = 0;
        mBits[idx++] = b;
        // std::cout << b << std::endl;
        while (leftMostBit > 1) {
            UnsignedIntegerType c = b & (b << (leftMostBit / 2));
            b = b ^ c ^ (c << (leftMostBit / 2) * (dimension - 1));
            mBits[idx++] = b;
            // std::cout << b << std::endl;
            leftMostBit /= 2;
        }
    }

    // For-loop basaed calculation with bugs (results are different from magic bits)
    UnsignedIntegerType getMortonKey_for(const _Vector& p, const _Vector& pMin,
                                         const _Vector& pMax) const {
        return getMortonKey_for(normalize(p, pMin - eps, pMax + eps));
    }

    UnsignedIntegerType getMortonKey_for(const _Vector& _p) const {
        UnsignedIntegerType result = 0;

        Vector<UnsignedIntegerType, Dimension> floored;
        for (int i = 0; i < Dimension; ++i) {
            floored[i] = UnsignedIntegerType(bigfloat(_p[i]) * bigfloat(maxCells));
        }

        for (int i = 0; i < numBits; i++) {
            for (int j = 0; j < Dimension; ++j) {
                result = result | ((floored[j] & (1ULL << i)) << (i + j));
            }
        }

        return result;
    }

    inline static _Vector normalize(const _Vector& point, const _Vector& pmin,
                                    const _Vector& pmax) {
        auto origin = (point - pmin);
        auto size = (pmax - pmin);

        // [0,0] ~ [1,1]
        return (origin / size);
    }
};

}  // namespace hj