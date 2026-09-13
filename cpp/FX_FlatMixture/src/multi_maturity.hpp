#pragma once

#include "lognormal_mixture.hpp"
#include "market_data.hpp"

#include <cstddef>
#include <vector>

namespace lnmix {

struct CalibratedSlice {
    MarketSlice market;
    CalibrationResult calibration;
};

class MultiMaturitySurface {
public:
    explicit MultiMaturitySurface(MarketSurface market);

    void calibrate(double eta = 0.10, int maxIter = 600, double tol = 1e-14);

    const MarketSurface& market() const { return market_; }
    const std::vector<CalibratedSlice>& slices() const { return slices_; }

    double call(double k, double t) const;
    double density(double k, double t) const;
    double impliedVol(double k, double t) const;
    double localVol(double k, double t) const;
    bool calendarArbitrageFree(double kMin = 0.75, double kMax = 1.25,
                               int nK = 151, double tol = 1e-10) const;

private:
    MarketSurface market_;
    std::vector<CalibratedSlice> slices_;

    std::size_t leftIndex(double t) const;
};

} // namespace lnmix
