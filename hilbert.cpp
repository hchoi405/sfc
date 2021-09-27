#include <bitset>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

using namespace std;

#define DIMS 3
#define BITS 5
#define PRINT_BITS DIMS + BITS + 10
#define PRINT true

const int newline = 0;

template <typename Type>
void printWithBits(const string& name, Type v, int lv = 0, string prefix = "") {
    if (!PRINT) return;
    if (name == "newline") {
        cout << endl;
        return;
    }
    if (lv > 4) cerr << "ERROR: Max lv is 4" << endl;

    string str = prefix;
    str.resize(10, ' ');

    for (int i = 0; i < lv; ++i) str += "  ";
    str += name + ": " + bitset<PRINT_BITS>(v).to_string();

    str += ", ";
    for (int i = 0; i < 4 - lv; ++i) str += "  ";

    int r = name.size() - 1;
    for (int i = 0; i < 4 - r; ++i) str += " ";

    str += to_string((int)(bitset<PRINT_BITS>(v).to_ulong()));

    cout << str << endl;
}

#define print(v, ...) printWithBits(#v, (v), ##__VA_ARGS__)

void FromAxes(uint32_t* X, int b, int n)  // position, #bits, dimension
{
    uint32_t M = 1 << (b - 1), P, Q, t;
    print(M, 0);
    int i;
    // Inverse undo
    for (Q = M; Q > 1; Q >>= 1) {
        P = Q - 1;
        print(Q, 1);
        print(P, 1);
        for (i = 0; i < n; i++) {
            print(i, 2);
            if (X[i] & Q) {
                print(X[0], 3, "before");
                X[0] ^= P;
                print(X[0], 3, "after");
            }
            // invert
            else {
                t = (X[0] ^ X[i]) & P;
                print(t, 3);
                print(X[0], 3, "before");
                X[0] ^= t;
                print(X[0], 3, "after");
                print(X[i], 3, "before");
                X[i] ^= t;
                print(X[i], 3, "after");
            }
        }
    }  // exchange
    // Gray encode

    print(newline);
    for (i = 1; i < n; i++) {
        print(i, 1);
        print(X[i], 1, "before");
        X[i] ^= X[i - 1];
        print(X[i], 1, "after");
    }

    print(newline);

    t = 0;
    for (Q = M; Q > 1; Q >>= 1) {
        print(Q, 1);
        if (X[n - 1] & Q) {
            print(t, 2, "before");
            t ^= Q - 1;
            print(t, 2, "after");
        }
    }

    print(newline);

    for (i = 0; i < n; i++) {
        print(i, 1);
        print(X[i], 1, "before");
        X[i] ^= t;
        print(X[i], 1, "after");
    }
}

uint64_t encode(uint32_t* X) {
    print(newline);

    cout << "--------------[ Encode ]--------------" << endl;
    print(X[0], 0, "initial");
    print(X[1], 0, "initial");
#if DIMS == 3
    print(X[2], 0, "initial");
#endif
    print(newline);

    FromAxes(X, BITS, DIMS);

    print(newline);

    print(X[0], 0, "final");
    print(X[1], 0, "final");
#if DIMS == 3
    print(X[2], 0, "final");
#endif

    // Interleaving the values of axes
    uint64_t start = 1 << (DIMS * BITS - 1), sum = 0;
    for (int i = BITS - 1; i >= 0; --i) {
        for (int j = 0; j < DIMS; ++j) {
            // cout << (X[j] >> i & 1);
            sum += (X[j] >> i & 1) * start;
            start >>= 1;
        }
    }

    print(newline);

    print(sum, 0);

    return sum;
}

void ToAxes(uint32_t* X, int b, int n)  // position, #bits, dimension
{
    uint32_t N = 2 << (b - 1), P, Q, t;
    int i;
    // Gray decode by H ^ (H/2)
    print(newline);

    t = X[n - 1] >> 1;
    print(N, 0);
    print(t, 0);
    for (i = n - 1; i >= 0; i--) {
        print(X[i], 1, "before");
        X[i] ^= X[i - 1];
        print(X[i], 1, "after");
    }
    print(X[0], 0, "before");
    X[0] ^= t;
    print(X[0], 0, "after");

    // Undo excess work
    print(newline);

    for (Q = 2; Q != N; Q <<= 1) {
        print(Q, 1);
        P = Q - 1;
        print(P, 1);
        for (i = n - 1; i >= 0; i--) {
            print(i, 2);
            if (X[i] & Q) {
                print(X[0], 3, "before");
                X[0] ^= P;
                print(X[0], 3, "after");
            }
            // invert
            else {
                t = (X[0] ^ X[i]) & P;
                X[0] ^= t;
                X[i] ^= t;

                print(t, 3);
                print(X[0], 3, "before");
                print(X[0], 3, "after");
                print(X[i], 3, "before");
                print(X[i], 3, "after");
            }
        }
    }  // exchange
}

void decode(uint32_t* X) {
    print(newline);

    cout << "--------------[ Decode ]--------------" << endl;
    print(X[0], 0, "initial");
    print(X[1], 0, "initial");
#if DIMS == 3
    print(X[2], 0, "initial");
#endif

    ToAxes(X, BITS, DIMS);

    // Bits
    // auto b =
    //     bitset<32>(bitset<32 - BITS * DIMS>(-1).to_ullong()) & bitset<32>(-1);
    // print(b.to_ullong());
    // X[0] &= b.to_ullong();
    // X[1] &= b.to_ullong();
    // X[2] &= b.to_ullong();

    print(newline);

    print(X[0], 0, "final");
    print(X[1], 0, "final");
#if DIMS == 3
    print(X[2], 0, "final");
#endif
}

void decode(uint64_t sum, uint32_t* X) {
    //
}

// rotate/flip a quadrant appropriately
void rot(int n, int* x, int* y, int rx, int ry) {
    if (ry == 0) {
        if (rx == 1) {
            *x = n - 1 - *x;
            *y = n - 1 - *y;
        }

        // Swap x and y
        int t = *x;
        *x = *y;
        *y = t;
    }
}

// convert (x,y) to d
int xy2d(int n, int x, int y) {
    int rx, ry, s, d = 0;
    for (s = n / 2; s > 0; s /= 2) {
        rx = (x & s) > 0;
        ry = (y & s) > 0;
        d += s * s * ((3 * rx) ^ ry);
        rot(n, &x, &y, rx, ry);
    }
    return d;
}

// convert d to (x,y)
void d2xy(int n, int d, int* x, int* y) {
    int rx, ry, s, t = d;
    *x = *y = 0;
    for (s = 1; s < n; s *= 2) {
        rx = 1 & (t / 2);
        ry = 1 & (t ^ rx);
        rot(s, x, y, rx, ry);
        *x += s * rx;
        *y += s * ry;
        t /= 4;
    }
}

int main() {
    uint32_t X[DIMS] = {5, 10, 20};

    uint64_t sum = encode(X);
    cout << "encoded: " << sum << endl;

    // sum to X

    // decod
    decode(X);
    cout << "decoded: ";
    for (int i = 0; i < DIMS; ++i) {
        cout << X[i];
        if (i < DIMS) cout << ", ";
    }

    // int ret = xy2d(8, 3, 5);
    // print(ret, 0);

    // vector<vector<int>> arr(8, vector<int>(8));
    // for (int i = 0; i < 8; ++i) {
    //     for (int j = 0; j < 8; ++j) {
    //         arr[i][j] = xy2d(8, i, j);
    //     }
    // }

    // for (int i = 0; i < 8; ++i) {
    //     for (int j = 0; j < 8; ++j) {
    //         cout << setw(2) << arr[i][j] << " ";
    //     }
    //     cout << endl;
    // }

    // X[0] = sum;
    // uint32_t* ret = decode(DIMS, BITS, X);

    // printf("%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d = 7865 check\n",
    //        X[0] >> 4 & 1, X[1] >> 4 & 1, X[2] >> 4 & 1,  // a
    //        X[0] >> 3 & 1, X[1] >> 3 & 1, X[2] >> 3 & 1,  // b
    //        X[0] >> 2 & 1, X[1] >> 2 & 1, X[2] >> 2 & 1,  // c
    //        X[0] >> 1 & 1, X[1] >> 1 & 1, X[2] >> 1 & 1,  // d
    //        X[0] >> 0 & 1, X[1] >> 0 & 1, X[2] >> 0 & 1);

    return 0;
}
