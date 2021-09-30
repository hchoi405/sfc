#include "hilbert.hpp"

#include <chrono>
#include <cmath>
#include <numeric>

using namespace hj;
using namespace std;

using hclock = chrono::high_resolution_clock;
using duration = chrono::duration<double>;

const int DIMS = 3;
const int BITS = 5;

void benchmark() {
    uint32_t X[DIMS];
    int maxValue = pow(2, BITS);
    for (int i = 0; i < DIMS; ++i) {
        X[i] = rand() % maxValue;
    }

    const size_t num_tests = 100000000;
    vector<uint64_t> results(num_tests);

    Hilbert hilbert(DIMS, BITS);

    auto start = hclock::now();
    for (size_t i = 0; i < num_tests; ++i) {
        uint64_t sum = hilbert.encode(X);
        results[i] = sum;
    }
    duration dur = hclock::now() - start;

    auto sum = accumulate(results.begin(), results.end(), 0);

    cout << "Took total " << dur.count() << "s for " << num_tests << " tests (" << sum << ")"
         << endl;
    cout << "Took avg.  " << dur.count() / num_tests << "s " << endl;
    cout << endl << endl;
}

int main() {
    uint32_t X[DIMS] = {5, 10, 20};

    Hilbert hilbert(DIMS, BITS);

    cout << "--------------[ Encode ]--------------" << endl;
    cout << "encoded: from [";
    for (int i = 0; i < DIMS; ++i) {
        cout << X[i];
        if (i < DIMS - 1) cout << ", ";
    }

    uint64_t sum = hilbert.encode(X);
    cout << "] to " << sum << endl;

    cout << "--------------[ Decode ]--------------" << endl;
    uint32_t X2[DIMS] = {0};
    hilbert.decode(sum, X2);
    cout << "decoded: from " << sum << " to [";
    for (int i = 0; i < DIMS; ++i) {
        cout << X2[i];
        if (i < DIMS - 1) cout << ", ";
    }
    cout << "]" << endl;

    cout << "--------------[ Benchmark ]--------------" << endl;
    benchmark();

    return 0;
}
