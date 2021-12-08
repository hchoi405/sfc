#ifndef HJ_HILBERT
#define HJ_HILBERT

#include <bitset>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include "sfc.h"

namespace sfc {
using namespace std;

template <int Dims, typename DataType, typename UInt = uint64_t>
class Hilbert : public SFC<Dims, DataType, UInt> {
    using Base = SFC<Dims, DataType, UInt>;

    // Using base's types
    using typename Base::_Point;

    // Using base's memeber variables
    using Base::numBitsPerAxis;
    using Base::numBitsTotal;
    using Base::numStrataPerAxis;

    // Own member variables

   public:
    Hilbert() {}

    UInt encode(const std::array<uint32_t, Dims>& x) const {
        FromAxes(const_cast<uint32_t*>(x.data()));
        // Interleaving the values of axes
        uint64_t sum = interleave(x.data());

        return sum;
    }
    UInt encode(const std::array<DataType, Dims>& x, const std::array<DataType, Dims> pMin,
                const std::array<DataType, Dims> pMax) const {
        std::array<uint32_t, Dims> uarr;
        for (int i = 0; i < Dims; ++i) uarr[i] = Base::normalize(x[i], pMin[i], pMax[i]);

        return encode(uarr);
    }
    void decode(const UInt v, std::array<uint32_t, Dims>& x) const {
        deinterleave(const_cast<uint32_t*>(x.data()), v);

        ToAxes(const_cast<uint32_t*>(x.data()));
    }
    void decode(const UInt v, std::array<DataType, Dims>& x, const std::array<DataType, Dims> pMin,
                const std::array<DataType, Dims> pMax) const {
        sfc::NotImplemented();
    }

   private:
    uint64_t interleave(const uint32_t* X) const {
        uint64_t start = 1ull << (Dims * numBitsPerAxis - 1), sum = 0;
        for (int i = numBitsPerAxis - 1; i >= 0; --i) {
            for (int j = 0; j < Dims; ++j) {
                // cout << (X[j] >> i & 1);
                sum += (X[j] >> i & 1) * start;
                start >>= 1;
            }
        }
        return sum;
    }

    void deinterleave(uint32_t* X, const uint64_t v) const {
        // De-interleaving the sum into axes values
        uint64_t src_pos = 1, dst_pos = 1;
        for (int i = 0; i < numBitsPerAxis; ++i) {
            for (int j = Dims - 1; j >= 0; --j) {
                X[j] ^= (-(uint32_t)((v & dst_pos) != 0) ^ X[j]) & src_pos;
                dst_pos <<= 1;
            }
            src_pos <<= 1;
        }
    }

    void FromAxes(uint32_t* X) const {
        uint32_t M = 1u << (numBitsPerAxis - 1), P, Q, t;
        int i;
        // Inverse undo
        for (Q = M; Q > 1; Q >>= 1) {
            P = Q - 1;
            for (i = 0; i < Dims; i++) {
                if (X[i] & Q) {
                    X[0] ^= P;
                }
                // invert
                else {
                    t = (X[0] ^ X[i]) & P;
                    X[0] ^= t;
                    X[i] ^= t;
                }
            }
        }  // exchange
        // Gray encode

        for (i = 1; i < Dims; i++) {
            X[i] ^= X[i - 1];
        }

        t = 0;
        for (Q = M; Q > 1; Q >>= 1) {
            if (X[Dims - 1] & Q) {
                t ^= Q - 1;
            }
        }

        for (i = 0; i < Dims; i++) {
            X[i] ^= t;
        }
    }

    void ToAxes(uint32_t* X) const {
        uint32_t N = 2u << (numBitsPerAxis - 1), P, Q, t;
        int i;
        // Gray decode by H ^ (H/2)

        t = X[Dims - 1] >> 1;
        for (i = Dims - 1; i > 0; i--) {
            X[i] ^= X[i - 1];
        }
        X[0] ^= t;

        // Undo excess work

        for (Q = 2; Q != N; Q <<= 1) {
            P = Q - 1;
            for (i = Dims - 1; i >= 0; i--) {
                if (X[i] & Q) {
                    X[0] ^= P;
                }
                // invert
                else {
                    t = (X[0] ^ X[i]) & P;
                    X[0] ^= t;
                    X[i] ^= t;
                }
            }
        }  // exchange
    }
};

};  // namespace sfc

#endif