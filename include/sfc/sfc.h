#ifndef HJ_SFC
#define HJ_SFC

#include <boost/multiprecision/cpp_int.hpp>
#include <cstdint>
#include <vector>

using namespace boost::multiprecision;

namespace sfc {

class NotImplemented {
   public:
    NotImplemented(const std::string& funcname = __builtin_FUNCTION()) {
        throw std::runtime_error("Not implemented: " + funcname);
    }
};

/**
 * @brief An abstract Space-Filling Curve class that has common variables for
 * all SFC methods like number of bits per dimension.
 *
 * @tparam DataType Type of original data (e.g., float array)
 * @tparam UInt Unsigned integer type for key (e.g., uint64_t)
 * @tparam Dims Number of dimensions (e.g., 3 for 3-D space)
 * @tparam Bits Number of bits per dimension (e.g., 10 bits for 3-D space uses total 30 bits)
 */
template <typename DataType = float, typename UInt = uint64_t, int Dims = 3, int Bits = 8>
class SFC {
   protected:
    ///////////////////////////////////////////////////////////////
    // Type Definitions
    ///////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////
    // Member variables
    ///////////////////////////////////////////////////////////////
    // clang-format off
	static const int NumBitsTotal =
		(std::is_same<UInt, uint8_t>::value) ? 8 :
		(std::is_same<UInt, uint16_t>::value) ? 16 :
		(std::is_same<UInt, uint32_t>::value) ? 32 :
		(std::is_same<UInt, uint64_t>::value) ? 64 :
		(std::is_same<UInt, uint128_t>::value) ? 128 : 0;
    // clang-format on
    static const UInt NumStrataPerAxis = UInt(1) << Bits;  // Number of strata per axis

    ///////////////////////////////////////////////////////////////
    // CHECKS
    ///////////////////////////////////////////////////////////////
    static_assert(Dims > 0, "Parameter 'Dims' must be > 0.");
    static_assert(std::is_constructible<UInt, DataType>::value,
                  "UInt instance cannot be constructed using DataType.");
    static_assert(Bits <= 32, "Bits should be less than 32.");
    static_assert(Dims * Bits <= NumBitsTotal);

   public:
    ///////////////////////////////////////////////////////////////
    // Functions (Implementations)
    ///////////////////////////////////////////////////////////////
    /**
     * @brief Encode given point as a 1-D Space-Filling key. This is a template
     * funciton to handle arbitrary type of array by converting it ot std::array<uint32_t, Dims>.
     *
     * @param x which has an operator[] and in range [0, numStrataPerAxis)
     * @return UInt
     */
    template <typename UIntPoint>
    UInt encode(const UIntPoint& x) const {
        std::array<uint32_t, Dims> uarr;
        for (int i = 0; i < Dims; ++i) uarr[i] = x[i];

        return encode(uarr);
    }

    /**
     * @brief Encode given ArrayLike type of point as a 1-D Space-Filling key
     *
     * @param x ArrayLike type which has an operator[] and in range [pMin, pMax]
     * @return UInt
     */
    template <typename FloatPoint>
    UInt encode(const FloatPoint& x, const FloatPoint& pMin, const FloatPoint& pMax) const {
        std::array<uint32_t, Dims> uarr;
        for (int i = 0; i < Dims; ++i) uarr[i] = normalize(x[i], pMin[i], pMax[i]);

        return encode(uarr);
    }

    /**
     * @brief Decode for UIntPoint type
     *
     * @tparam UIntPoint array-like point type
     * @param v Code encoded by space-filling curve
     * @param x Output int [0, NumStrataPerAxis)
     */
    template <typename UIntPoint>
    void decode(const UInt v, UIntPoint& x) const {
        std::array<uint32_t, Dims> uarr{};
        decode(v, uarr);

        for (int i = 0; i < Dims; ++i) x[i] = uarr[i];
    }

    /**
     * @brief Decode for UIntPoint type
     *
     * @tparam UIntPoint array-like point type
     * @param v Code encoded by space-filling curve
     * @param x Output int [0, NumStrataPerAxis)
     */
    template <typename UIntPoint>
    void decode(const UInt v, UIntPoint* x) const {
        std::array<uint32_t, Dims> uarr{};
        decode(v, uarr);

        for (int i = 0; i < Dims; ++i) x[i] = uarr[i];
    }

    /**
     * @brief Decode for FloatPoint type
     *
     * @tparam FloatPoint point type with float elements
     * @param v Code encoded by space-filling curve
     * @param x Output int [pMin, pMax]
     * @param pMin Minimum value for output
     * @param pMax Maximum value for output
     */
    template <typename FloatPoint>
    void decode(const UInt v, FloatPoint& x, const FloatPoint& pMin, const FloatPoint& pMax) const {
        std::array<uint32_t, Dims> uarr{};
        decode(v, uarr);

        for (int i = 0; i < Dims; ++i) x[i] = denormalize(uarr[i], pMin[i], pMax[i]);
    }

    /**
     * @brief Normalization value in range [pMin, pMax] to [0, numStrataPerAxis-1]
     *
     * @param val
     * @param pMin
     * @param pMax
     * @return T
     */
    template <typename T>
    inline static uint32_t normalize(const T& val, const T& pMin, const T& pMax) {
        assert(val >= pMin && val <= pMax);
        auto origin = (val - pMin);
        auto size = (pMax - pMin);

        return std::min(uint32_t((origin / size) * NumStrataPerAxis),
                        uint32_t(NumStrataPerAxis - 1));
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
    template <typename T>
    inline static T denormalize(const uint32_t& val, const T& pMin, const T& pMax) {
        T v = T(val) / NumStrataPerAxis;
        v *= (pMax - pMin);
        v += pMin;
        return v;
    }

    ///////////////////////////////////////////////////////////////
    // Functions (Interfaces)
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
    virtual UInt encode(const std::array<DataType, Dims>& x, const std::array<DataType, Dims>& pMin,
                        const std::array<DataType, Dims>& pMax) const = 0;

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
                        const std::array<DataType, Dims>& pMin,
                        const std::array<DataType, Dims>& pMax) const = 0;
};

}  // namespace sfc

#endif