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

#include "hilbert.hpp"
#include "utils.hpp"
#include "zorder_knn/less.hpp"

using hclock = std::chrono::high_resolution_clock;
using duration = std::chrono::duration<float>;

template <typename T, std::size_t d>
using Pt = std::array<T, d>;
using Pt2 = Pt<float, 2>;
using Pt2ui = Pt<uint32_t, 2>;
const float MAX_VAL = 1;
const float MIN_VAL = -1;

template <typename... Args>
std::string string_format(const std::string& format, Args... args) {
    size_t size = snprintf(nullptr, 0, format.c_str(), args...) + 1;
    if (size <= 0) {
        throw std::runtime_error("Error during formatting.");
    }
    std::unique_ptr<char[]> buf(new char[size]);
    snprintf(buf.get(), size, format.c_str(), args...);
    return std::string(buf.get(), buf.get() + size - 1);
}

std::string svg_circle(std::size_t id, Pt2 const& c_xy, float r = 0.2f) {
    std::stringstream circle;
    circle << "<circle id=\"circle" << id << "\" cx=\"" << c_xy[0] << "\" cy=\"" << c_xy[1]
           << "\" r=\"" << r << "\"/>";
    return circle.str();
}

std::string svg_path(std::size_t id, float s, Pt2 const& p0, Pt2 const& p1, float w = 0.2) {
    std::stringstream path;
    path << "<path id=\"path" << id << "\" style=\"stroke:#" << std::setfill('0') << std::hex
         << std::setw(2) << static_cast<int>(255.0f * s) << std::setw(2) << 0u << std::setw(2)
         << static_cast<int>(255.0f * (1.0 - s)) << ";stroke-width:" << w << "\" d=\"M " << p0[0]
         << "," << p0[1] << " " << p1[0] << "," << p1[1] << "\"/>";
    return path.str();
}

void svg_save(std::string filename, std::vector<Pt2> const& pts) {
    auto it = std::max_element(pts.begin(), pts.end(), [](const Pt2& lhs, const Pt2& rhs) {
        return *std::max_element(lhs.begin(), lhs.end()) <
               *std::max_element(rhs.begin(), rhs.end());
    });
    auto maxVal = *std::max_element(it->begin(), it->end());
    auto radius = maxVal * 0.01f;
    auto width = maxVal * 0.02f;
    std::cout << "maxVal: " << maxVal << std::endl;
    std::cout << "circle radius: " << radius << std::endl;
    std::cout << "path width: " << width << std::endl;

    const std::string header = string_format(
        "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n"
        "<svg id=\"svg8\""
        " xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\""
        " xmlns=\"http://www.w3.org/2000/svg\""
        " width=\"1024px\" height=\"1024px\""
        " version=\"1.1\""
        " xmlns:cc=\"http://creativecommons.org/ns#\""
        " xmlns:dc=\"http://purl.org/dc/elements/1.1/\""
        " viewBox=\"%f %f %f %f\">\n"
        "  <metadata id=\"metadata5\">\n"
        "    <rdf:RDF>\n"
        "      <cc:Work rdf:about=\"\">\n"
        "        <dc:format>image/svg+xml</dc:format>\n"
        "        <dc:type"
        " rdf:resource=\"http://purl.org/dc/dcmitype/StillImage\" />\n"
        "        <dc:title></dc:title>\n"
        "      </cc:Work>\n"
        "    </rdf:RDF>\n"
        "  </metadata>\n"
        "<rect x=\"%f\" y=\"%f\" width=\"200%\" height=\"200%\" "
        "fill=\"white\"/>\n"
        "<g id=\"layer1\">\n",
        -1.1f * maxVal, -1.1f * maxVal, 2 * 1.1f * maxVal, 2 * 1.1f * maxVal);

    const std::string footer =
        "</g>\n"
        "</svg>\n";

    std::ofstream svg(filename);
    svg << header;
    for (std::size_t i(0); i < pts.size() - 1; ++i) {
        svg << "  "
            << svg_path(i, static_cast<float>(i) / (pts.size() - 2), pts[i], pts[i + 1], width)
            << "\n";
    }
    for (std::size_t i(0); i < pts.size(); ++i) {
        svg << "  " << svg_circle(i, pts[i], radius) << "\n";
    }
    svg << footer;
    svg.close();
}

void convertToHilbert(const std::vector<Pt2>& pts, std::vector<uint64_t>& codes, const int nDims,
                      const int nBits) {
    hj::Hilbert hilbert(nDims, nBits);
    const uint32_t maxVal = pow(2, nBits);
    hj::print(maxVal);

    // Discretize points from [MIN_VAL, MAX_VAL] to [0, maxVal-1]
    const float extent = MAX_VAL - MIN_VAL;
    std::vector<Pt2ui> pts_ui;
    for (const auto& p : pts) {
        auto x = (p[0] - MIN_VAL) / extent;              // [0, 1]
        x = std::min(uint32_t(x * maxVal), maxVal - 1);  // [0, maxVal-1]
        auto y = (p[1] - MIN_VAL) / extent;              // [0, 1]
        y = std::min(uint32_t(y * maxVal), maxVal - 1);  // [0, maxVal-1]
        pts_ui.push_back({(uint32_t)x, (uint32_t)y});
    }

    // Convert points to Hilbert codes
    codes.clear();
    codes.reserve(codes.size());
    for (auto& p : pts_ui) {
        uint64_t code = hilbert.encode(p.data());
        codes.push_back(code);
    }
}

int main(int argc, char* argv[]) {
    std::random_device rd;
    std::mt19937 e2(0);
    std::uniform_real_distribution<float> dist(MIN_VAL, MAX_VAL);

    // Generate random points in float
    std::vector<Pt2> pts(10000);
    std::generate(pts.begin(), pts.end(), [&] {
        Pt2 p;
        std::generate(p.begin(), p.end(), [&] { return dist(e2); });
        return p;
    });

    std::vector<uint64_t> codes;
    const int nDims = 2;
    const int nBits = 16;

    using Elem = std::pair<Pt2, uint64_t>;
    std::vector<Elem> pairs;
    pairs.reserve(pts.size());

    auto start = hclock::now();
    convertToHilbert(pts, codes, nDims, nBits);
    for (size_t i = 0; i < codes.size(); ++i) {
        pairs.push_back({pts[i], codes[i]});
    }
    std::sort(pairs.begin(), pairs.end(),
              [](const auto& lhs, const auto& rhs) { return lhs.second < rhs.second; });

    std::vector<Pt2> sorted_pts(pts.size());
    for (size_t i = 0; i < pairs.size(); ++i) {
        sorted_pts[i] = pairs[i].first;
    }
    duration dur = hclock::now() - start;
    std::cout << "Took " << dur.count() << "s for Hilbert Curve" << std::endl;

    svg_save("example_random.svg", pts);
    svg_save("example_hilbert.svg", sorted_pts);

    start = hclock::now();
    std::sort(pts.begin(), pts.end(), zorder_knn::Less<Pt2, 2>());
    dur = hclock::now() - start;
    std::cout << "Took " << dur.count() << "s for zCurve" << std::endl;
    svg_save("example_morton.svg", pts);

    return EXIT_SUCCESS;
}