#pragma once

#include "market_data.hpp"

#include <cstddef>
#include <vector>

namespace lnmix {

// Natural cubic spline y(x), with linear extrapolation outside the knot range.
class NaturalCubicSpline {
public:
    void fit(const std::vector<double>& x, const std::vector<double>& y);
    bool empty() const { return x_.empty(); }
    double value(double x) const;
    double firstDerivative(double x) const;
    double secondDerivative(double x) const;

private:
    std::vector<double> x_;
    std::vector<double> a_; // interval constant
    std::vector<double> b_; // interval linear
    std::vector<double> c_; // interval quadratic
    std::vector<double> d_; // interval cubic

    std::size_t interval(double x) const;
};

struct CubicVarianceSlice {
    MarketSlice market;
    NaturalCubicSpline spline; // w(x) = sigma_imp^2 * T, x = log(K/F)

    double rawTotalVariance(double k) const;
    double totalVariance(double k) const;
    double annualVol(double k) const;
    double call(double k) const;
    double butterflyFactor(double k) const;
    double density(double k) const; // d^2 normalized call / dk^2
};

class CubicVarianceSurface {
public:
    explicit CubicVarianceSurface(MarketSurface market);

    void calibrate(); // exact interpolation of each maturity slice

    const MarketSurface& market() const { return market_; }
    const std::vector<CubicVarianceSlice>& slices() const { return slices_; }

    double call(double k, double t) const;
    double density(double k, double t) const;
    double impliedVol(double k, double t) const;
    double localVol(double k, double t) const;
    bool calendarArbitrageFree(double kMin = 0.75, double kMax = 1.25,
                               int nK = 151, double tol = 1e-10) const;

private:
    MarketSurface market_;
    std::vector<CubicVarianceSlice> slices_;
    std::size_t leftIndex(double t) const;
};

} // namespace lnmix
