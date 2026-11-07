#pragma once

#include <string>
#include <vector>

namespace lnmix {

struct MarketPoint {
    double k = 0.0;              // K / Forward
    double strike = 0.0;         // physical strike if available, otherwise k
    double marketVol = 0.0;      // annualized Black implied vol
    double sourcePrice = 0.0;    // source option mid, NaN if unavailable
    double sourceIv = 0.0;       // source-displayed IV, NaN if unavailable
    std::string pillar;          // strike or FX delta label
};

struct MarketSlice {
    std::string instrument;
    std::string asOf;
    std::string maturity;
    double t = 0.0;
    double forward = 1.0;
    std::vector<MarketPoint> points;
};

struct MarketSurface {
    std::string instrument;
    std::string asOf;
    std::string provenance;
    std::vector<MarketSlice> slices;
};

MarketSurface loadStrikeSurfaceCsv(const std::string& path);
MarketSurface loadFxDeltaSurfaceCsv(const std::string& path);

} // namespace lnmix
