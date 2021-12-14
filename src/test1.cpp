#include <bitset>
#include <chrono>
#include <cmath>
#include <numeric>
#include <random>

#include "sfc/hilbert.hpp"
#include "sfc/utils.hpp"

using namespace sfc;
using namespace std;

using hclock = chrono::high_resolution_clock;
using duration = chrono::duration<double>;

const int DIMS = 2;
const int BITS = 10;
using CodeType = uint32_t;

using Curve = SFC<float, CodeType, DIMS, BITS>;

void benchmark() {
    // uint32_t X[DIMS];
    std::array<uint32_t, 2> X;
    int maxValue = pow(2, BITS);
    for (int i = 0; i < DIMS; ++i) {
        X[i] = rand() % maxValue;
    }

    const size_t num_tests = 100000000;
    vector<CodeType> results(num_tests);

    std::unique_ptr<Curve> hilbert(new Hilbert<float, CodeType, DIMS, BITS>());

    auto start = hclock::now();
    for (size_t i = 0; i < num_tests; ++i) {
        CodeType sum = hilbert->encode(X);
        results[i] = sum;
    }
    duration dur = hclock::now() - start;

    auto sum = accumulate(results.begin(), results.end(), 0);

    cout << "Took total " << dur.count() << "s for " << num_tests << " tests (" << sum << ")"
         << endl;
    cout << "Took avg.  " << dur.count() / num_tests << "s " << endl;
    cout << endl << endl;
}

void printEncoded(const Curve* hilbert, std::vector<uint32_t> X) {
    printArray(X, DIMS);
    std::cout << " -> " << hilbert->encode(X.data()) << ",\t"
              << std::bitset<64>(hilbert->encode(X.data())) << std::endl;
}

void printDecoded(const Curve* hilbert, CodeType sum) {
    std::vector<uint32_t> X(DIMS);
    cout << sum << ",\t" << std::bitset<64>(sum) << " -> ";
    hilbert->decode(sum, X.data());
    printArray(X, DIMS);
    std::cout << endl;
}

int main() {
    // uint32_t X[DIMS] = {5, 10};
    // uint32_t X[DIMS] = {2, 0, 0, 1, 5};  // (5 dims, 8 bits): 1504
    // uint32_t X[DIMS] = {0,3,5,1,0};  // (5 dims, 8 bits): 1504
    const int MAX_VAL = pow(2, BITS);
    std::unique_ptr<Curve> hilbert(new Hilbert<float, CodeType, DIMS, BITS>());
    // uint32_t X[DIMS];
    std::vector<uint32_t> X(DIMS);

    std::mt19937 rand;
    std::uniform_int_distribution<uint32_t> dist(0, MAX_VAL);

    for (int i = 0; i < DIMS; ++i) {
        X[i] = dist(rand);
    }

    cout << "--------------[ Encode ]--------------" << endl;
    cout << "encoded: from [";
    for (int i = 0; i < DIMS; ++i) {
        cout << X[i];
        if (i < DIMS - 1) cout << ", ";
    }
    CodeType sum = hilbert->encode(X.data());
    cout << "] to " << sum << endl;
    // printEncoded(hilbert.get(), {7, 5});

    cout << "--------------[ Decode ]--------------" << endl;
    uint32_t X2[DIMS] = {0};
    hilbert->decode(sum, X2);
    cout << "decoded: from " << sum << " to [";
    for (int i = 0; i < DIMS; ++i) {
        cout << X2[i];
        if (i < DIMS - 1) cout << ", ";
    }
    cout << "]" << endl;
    // for (int i = 0; i < 16; ++i) printDecoded(hilbert.get(), i);
    // printDecoded(hilbert, 68719476736ull);
    // printDecoded(hilbert, 137438953472ull);
    // printDecoded(hilbert, 206158430208ull);
    // printDecoded(hilbert, 1099511627776ull-1);
    // printDecoded(hilbert, 9674);
    // printDecoded(hilbert, 16107);
    // printDecoded(hilbert, 18131);
    // printDecoded(hilbert, 22734);
    // printDecoded(hilbert, 24919);
    // printDecoded(hilbert, 30550);
    // printDecoded(hilbert, 34554);
    // printDecoded(hilbert, 40510);
    // printDecoded(hilbert, 43305);
    // printDecoded(hilbert, 48162);
    // printDecoded(hilbert, 49923);
    // printDecoded(hilbert, 54480);
    // printDecoded(hilbert, 59378);
    // printDecoded(hilbert, 62076);

    // cout << "--------------[ Benchmark ]--------------" << endl;
    // benchmark();

    return 0;
}
