#include <algorithm>
#include <array>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <vector>

#include "sfc/hilbert.hpp"
#include "sfc/svg.h"
#include "sfc/utils.hpp"
#include "sfc/zcurve.hpp"
#include "zorder_knn/less.hpp"

using hclock = std::chrono::high_resolution_clock;
using duration = std::chrono::duration<float>;

// Configs
using PointType = float;
using Pt2 = std::array<PointType, 2>;
using CodeType = uint32_t;
const PointType MAX_VAL = 1;
const PointType MIN_VAL = -1;
const int nDims = 2;
const int nBits = 5;

// Type alias
using Pt2ui = std::array<uint32_t, 2>;
using Elem = std::pair<Pt2, CodeType>;

using SFC23 = sfc::SFC<PointType, CodeType, nDims, nBits>;
using Hilbert23 = sfc::Hilbert<PointType, CodeType, nDims, nBits>;
using Zcurve23 = sfc::Zcurve<PointType, CodeType, nDims, nBits>;

void normalizePoints(std::unique_ptr<SFC23>& sfc, const std::vector<Pt2>& pts,
                     std::vector<Pt2ui>& out) {
    const uint32_t maxVal = uint32_t(1) << nBits;

    // Discretize points from [MIN_VAL, MAX_VAL] to [0, maxVal-1]
    const float extent = MAX_VAL - MIN_VAL;
    out.clear();
    out.reserve(out.size());
    for (const auto& p : pts) {
        auto x = sfc->normalize(p[0], MIN_VAL, MAX_VAL);
        auto y = sfc->normalize(p[1], MIN_VAL, MAX_VAL);
        out.push_back({(uint32_t)x, (uint32_t)y});
    }
}

void pointToHilbert(const std::vector<Pt2>& pts, std::vector<CodeType>& codes) {
    std::unique_ptr<SFC23> hilbert = std::make_unique<Hilbert23>();

    std::vector<Pt2ui> pts_ui;
    normalizePoints(hilbert, pts, pts_ui);

    // Convert points to Hilbert codes
    codes.clear();
    codes.reserve(codes.size());
    for (auto& p : pts_ui) {
        sfc::printArray(p, 2);
        std::cout << std::endl;
        CodeType code = hilbert->encode(p);
        codes.push_back(code);
    }
}

void hilbertToPoint(const std::vector<CodeType>& codes, std::vector<Pt2>& pts) {
    std::unique_ptr<SFC23> hilbert = std::make_unique<Hilbert23>();

    pts.resize(codes.size());

    for (size_t i = 0; i < codes.size(); ++i) {
        Pt2ui up{};
        hilbert->decode(codes[i], up);
        Pt2 p{(float)up[0], (float)up[1]};
        pts[i][0] = hilbert->denormalize(p[0], MIN_VAL, MAX_VAL);
        pts[i][1] = hilbert->denormalize(p[1], MIN_VAL, MAX_VAL);

        sfc::printArray(up, 2);
        std::cout << std::endl;
    }
}

void pointToMorton(const std::vector<Pt2>& pts, std::vector<CodeType>& codes) {
    std::unique_ptr<SFC23> zcurve = std::make_unique<Zcurve23>();

    std::vector<Pt2ui> pts_ui;
    normalizePoints(zcurve, pts, pts_ui);

    // Convert points to Morton codes
    codes.clear();
    codes.reserve(codes.size());
    for (auto& p : pts_ui) {
        CodeType code = zcurve->encode(p);
        codes.push_back(code);
    }
}

void sortByCodes(const std::vector<Pt2>& pts, const std::vector<CodeType>& codes,
                 std::vector<Pt2>& out) {
    std::vector<Elem> pairs;
    pairs.clear();
    pairs.reserve(pts.size());
    for (size_t i = 0; i < codes.size(); ++i) {
        pairs.push_back({pts[i], codes[i]});
    }
    std::sort(pairs.begin(), pairs.end(),
              [](const auto& lhs, const auto& rhs) { return lhs.second < rhs.second; });

    out.resize(pts.size());
    for (size_t i = 0; i < pairs.size(); ++i) {
        out[i] = pairs[i].first;
    }
}

void testEncode() {
    std::random_device rd;
    std::mt19937 e2(0);
    std::uniform_real_distribution<float> dist(MIN_VAL, MAX_VAL);

    // Generate random points in float
    std::vector<Pt2> pts(100);
    std::generate(pts.begin(), pts.end(), [&] {
        Pt2 p;
        std::generate(p.begin(), p.end(), [&] { return dist(e2); });
        return p;
    });
    svg_save("example_random.svg", pts);

    std::vector<CodeType> codes;

    std::vector<Pt2> sorted_hilbert, sorted_zcurve;
    sorted_hilbert.reserve(pts.size());
    sorted_zcurve.reserve(pts.size());

    /////////////////////////////////////////////////////////////
    // Hilbert Curve
    /////////////////////////////////////////////////////////////
    auto start = hclock::now();
    pointToHilbert(pts, codes);
    sortByCodes(pts, codes, sorted_hilbert);
    duration dur = hclock::now() - start;

    std::cout << "Took " << dur.count() << "s for Hilbert Curve" << std::endl;
    svg_save("example_hilbert.svg", sorted_hilbert);

    /////////////////////////////////////////////////////////////
    // Z-Curve
    /////////////////////////////////////////////////////////////
    start = hclock::now();
    pointToMorton(pts, codes);
    sortByCodes(pts, codes, sorted_zcurve);
    dur = hclock::now() - start;

    std::cout << "Took " << dur.count() << "s for zCurve" << std::endl;
    svg_save("example_morton.svg", sorted_zcurve);
}

void testEncodeDecode() {
    std::random_device rd;
    std::mt19937 e2(0);
    std::uniform_real_distribution<float> dist(MIN_VAL, MAX_VAL);

    // Generate random points in float
    // std::vector<Pt2> pts(4);
    // std::generate(pts.begin(), pts.end(), [&] {
    //     Pt2 p;
    //     std::generate(p.begin(), p.end(), [&] { return dist(e2); });
    //     return p;
    // });

    std::vector<Pt2> pts;
    pts.push_back({-1.0f, -1.0f});
    pts.push_back({1.0f, -1.0f});
    pts.push_back({1.0f, 1.0f});
    pts.push_back({-1.0f, 1.0f});
    svg_save("example_random.svg", pts);

    std::vector<CodeType> codes;

    std::vector<Pt2> sorted, decoded;
    sorted.reserve(pts.size());

    /////////////////////////////////////////////////////////////
    // Hilbert Curve
    /////////////////////////////////////////////////////////////
    auto start = hclock::now();
    duration dur = hclock::now() - start;

    pointToHilbert(pts, codes);
    sortByCodes(pts, codes, sorted);
    hilbertToPoint(codes, decoded);

    for (size_t i = 0; i < sorted.size(); ++i) {
        sfc::printArray(pts[i], nDims);
        std::cout << " -> ";
        sfc::printArray(decoded[i], nDims);
        std::cout << std::endl;
    }

    std::cout << "Took " << dur.count() << "s for Hilbert Curve" << std::endl;
    svg_save("example_hilbert.svg", sorted);
}

void testDecode() {
    std::random_device rd;
    std::mt19937 e2(0);

    CodeType maxVal = -1;
    std::uniform_int_distribution<CodeType> dist(0, maxVal);

    // Generate random points in float
    std::vector<CodeType> codes(100);
    std::generate(codes.begin(), codes.end(), [&] { return dist(e2); });
    std::sort(codes.begin(), codes.end());
    // svg_save("example_random.svg", pts);

    // std::vector<uint64_t> codes;

    std::vector<Pt2> hilbert, zcurve;
    hilbert.reserve(codes.size());
    zcurve.reserve(codes.size());

    /////////////////////////////////////////////////////////////
    // Hilbert Curve
    /////////////////////////////////////////////////////////////
    auto start = hclock::now();
    hilbertToPoint(codes, hilbert);
    duration dur = hclock::now() - start;

    std::cout << "Took " << dur.count() << "s for Hilbert Curve" << std::endl;
    svg_save("example_hilbert.svg", hilbert);

    // /////////////////////////////////////////////////////////////
    // // Z-Curve
    // /////////////////////////////////////////////////////////////
    // start = hclock::now();
    // pointToMorton(pts, codes);
    // sortByCodes(pts, codes, sorted_zcurve);
    // dur = hclock::now() - start;

    // std::cout << "Took " << dur.count() << "s for zCurve" << std::endl;
    // svg_save("example_morton.svg", sorted_zcurve);
}

int main(int argc, char* argv[]) {
    // testDecode();
    testEncodeDecode();

    return EXIT_SUCCESS;
}