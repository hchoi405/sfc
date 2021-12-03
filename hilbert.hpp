#include <bitset>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

namespace hj {
using namespace std;

class Hilbert {
   public:
    Hilbert() : dims(0), bits(0) {}
    Hilbert(int dims, int bits) : dims(dims), bits(bits) {}

    void init(int dims, int bits) {
        this->dims = dims;
        this->bits = bits;
    }

    // Decode de-interleaved number
    void decode(uint32_t* X) { ToAxes(X); }

    // De-interleave and decode it into axes
    void decode(const uint64_t v, uint32_t* X) {
        deinterleave(X, v);

        decode(X);
    }

    // Encode axes into hilbert code
    uint64_t encode(uint32_t* X) {
        FromAxes(X);

        // Interleaving the values of axes
        uint64_t sum = interleave(X);

        return sum;
    }

   private:
    uint64_t interleave(const uint32_t* X) {
        uint64_t start = 1ull << (dims * bits - 1), sum = 0;
        for (int i = bits - 1; i >= 0; --i) {
            for (int j = 0; j < dims; ++j) {
                // cout << (X[j] >> i & 1);
                sum += (X[j] >> i & 1) * start;
                start >>= 1;
            }
        }
        return sum;
    }

    void deinterleave(uint32_t* X, const uint64_t v) {
        // De-interleaving the sum into axes values
        uint64_t src_pos = 1, dst_pos = 1;
        for (int i = 0; i < bits; ++i) {
            for (int j = dims - 1; j >= 0; --j) {
                X[j] ^= (-(uint32_t)((v & dst_pos) != 0) ^ X[j]) & src_pos;
                dst_pos <<= 1;
            }
            src_pos <<= 1;
        }
    }

    void FromAxes(uint32_t* X) {
        uint32_t M = 1u << (bits - 1), P, Q, t;
        int i;
        // Inverse undo
        for (Q = M; Q > 1; Q >>= 1) {
            P = Q - 1;
            for (i = 0; i < dims; i++) {
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

        for (i = 1; i < dims; i++) {
            X[i] ^= X[i - 1];
        }

        t = 0;
        for (Q = M; Q > 1; Q >>= 1) {
            if (X[dims - 1] & Q) {
                t ^= Q - 1;
            }
        }

        for (i = 0; i < dims; i++) {
            X[i] ^= t;
        }
    }

    void ToAxes(uint32_t* X) {
        uint32_t N = 2u << (bits - 1), P, Q, t;
        int i;
        // Gray decode by H ^ (H/2)

        t = X[dims - 1] >> 1;
        for (i = dims - 1; i > 0; i--) {
            X[i] ^= X[i - 1];
        }
        X[0] ^= t;

        // Undo excess work

        for (Q = 2; Q != N; Q <<= 1) {
            P = Q - 1;
            for (i = dims - 1; i >= 0; i--) {
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

   private:
    int dims, bits;
};

};  // namespace hj