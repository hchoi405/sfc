#ifndef HJ_SFC
#define HJ_SFC

#include <stdint.h>

#include <vector>

#include "Vector.h"

namespace sfc {

/**
 * @brief An abstract Space-Filling Curve class that has common variables for
 * all SFC methods like number of bits per dimension.
 *
 * @tparam Dims Number of dimensions (e.g., 3 for 3-D space)
 * @tparam DataType Type of original data (e.g., float array)
 * @tparam UInt Unsigned integer type for key (e.g., uint64_t)
 */
template <int Dims = 3, typename DataType = float, typename UInt = uint64_t>
class SFC {
    ///////////////////////////////////////////////////////////////
    // CHECKS
    ///////////////////////////////////////////////////////////////
    static_assert(Dims > 0, "Parameter 'Dims' must be > 0.");
    static_assert(std::is_constructible<UInt, DataType>::value,
                  "UInt instance cannot be constructed using DataType.");

   protected:
    ///////////////////////////////////////////////////////////////
    // Type Definitions
    ///////////////////////////////////////////////////////////////
    // Point type used as a wrapper for various point classes (e.g., PBRT's Point3f)
    using _Point = Vector<DataType, Dims>;

   protected:
    ///////////////////////////////////////////////////////////////
    // Member variables
    ///////////////////////////////////////////////////////////////
    // clang-format off
	static const int numBitsTotal =
		(std::is_same<UInt, uint8_t>::value) ? 8 :
		(std::is_same<UInt, uint16_t>::value) ? 16 :
		(std::is_same<UInt, uint32_t>::value) ? 32 :
		(std::is_same<UInt, uint64_t>::value) ? 64 :
		(std::is_same<UInt, uint128_t>::value) ? 128 : 0;
    // clang-format on
    static const int numBitsPerAxis = numBitsTotal / Dims;  // Maximum number of bits per axis
    static const UInt numStrataPerAxis = UInt(1) << numBitsPerAxis;  // Number of strata per axis

   public:
    ///////////////////////////////////////////////////////////////
    // Interfaces
    ///////////////////////////////////////////////////////////////
    /**
     * @brief Encode given uint array as a 1-D Space-Filling key
     *
     * @param x Array of uint32_t type in range [0, numStrataPerAxis)
     * @return UInt
     */
    virtual UInt encode(const std::array<uint32_t, Dims>& x) const = 0;

    /**
     * @brief Encode given DataType array as a 1-D Space-Filling key
     *
     * @param x
     * @param pMin
     * @param pMax
     * @return UInt
     */
    virtual UInt encode(const std::array<DataType, Dims>& x, const std::array<DataType, Dims> pMin,
                        const std::array<DataType, Dims> pMax) const = 0;

    /**
     * @brief Decode given 1-D Space-Filling key to uint array
     *
     * @param v
     * @param x
     */
    virtual void decode(const UInt v, std::array<uint32_t, Dims>& x) const = 0;

    /**
     * @brief Decode given 1-D Space-Filling key to DataType array
     *
     * @param v
     * @param x
     * @param pMin
     * @param pMax
     */
    virtual void decode(const UInt v, std::array<DataType, Dims>& x,
                        const std::array<DataType, Dims> pMin,
                        const std::array<DataType, Dims> pMax) const = 0;

   protected:
    /**
     * @brief Normalization value in range [pMin, pMax] to [0, numStrataPerAxis-1]
     *
     * @param val
     * @param pMin
     * @param pMax
     * @return T
     */
    template <typename T>
    inline static T normalize(const T& val, const T& pMin, const T& pMax) {
        assert(val >= pMin && val <= pMax);
        auto origin = (val - pMin);
        auto size = (pMax - pMin);

        // [0,0] ~ [1,1]
        return std::min((origin / size) * numStrataPerAxis, T(numStrataPerAxis - 1));
    }

    /**
     * @brief De-normalize value in range [0, numStrataPerAxis - 1] to [pMin, pMax]
     *
     * @tparam T
     * @param val
     * @param pMin
     * @param pMax
     * @return T
     */
    // template <typename T>
    // inline static T denormalize(const T& val, const T& pMin, const T& pMax) {
    //     // Not impelmented
    //     return T(-1);
    // }
};

class NotImplemented {
   public:
    NotImplemented(const std::string& funcname = __builtin_FUNCTION()) {
        throw std::runtime_error("Not implemented: " + funcname);
    }
};

}  // namespace sfc

#endif