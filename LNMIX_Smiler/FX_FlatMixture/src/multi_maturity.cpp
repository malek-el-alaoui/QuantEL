#include "multi_maturity.hpp"

#include <cstddef>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace lnmix {

MultiMaturitySurface::MultiMaturitySurface(MarketSurface market)
    : market_(std::move(market)) {
    std::sort(market_.slices.begin(), market_.slices.end(),
              [](const MarketSlice& a, const MarketSlice& b) { return a.t < b.t; });
}

void MultiMaturitySurface::calibrate(double eta, int maxIter, double tol) {
    slices_.clear();
    LognormalMixtureCalibrator calibrator;
    for (const auto& s : market_.slices) {
        std::vector<std::pair<double, double>> q;
        q.reserve(s.points.size());
        for (const auto& p : s.points) q.emplace_back(p.k, p.marketVol * std::sqrt(s.t));
        auto result = calibrator.calibrate(q, eta, maxIter, tol);
        slices_.push_back({s, std::move(result)});
    }
}

std::size_t MultiMaturitySurface::leftIndex(double t) const {
    if (slices_.empty()) throw std::runtime_error("Surface is not calibrated");
    if (slices_.size() == 1 || t <= slices_.front().market.t) return 0;
    if (t >= slices_.back().market.t) return slices_.size() - 1;
    for (std::size_t i = 0; i + 1 < slices_.size(); ++i)
        if (t >= slices_[i].market.t && t <= slices_[i + 1].market.t) return i;
    return slices_.size() - 1;
}

double MultiMaturitySurface::call(double k, double t) const {
    const std::size_t i = leftIndex(t);
    if (slices_.size() == 1 || i == slices_.size() - 1 || t <= slices_.front().market.t)
        return LognormalMixtureCalibrator::mixtureCall(slices_[i].calibration.components, k);
    const auto& a = slices_[i];
    const auto& b = slices_[i + 1];
    const double w = (t - a.market.t) / (b.market.t - a.market.t);
    const double ca = LognormalMixtureCalibrator::mixtureCall(a.calibration.components, k);
    const double cb = LognormalMixtureCalibrator::mixtureCall(b.calibration.components, k);
    return (1.0 - w) * ca + w * cb;
}

double MultiMaturitySurface::density(double k, double t) const {
    const std::size_t i = leftIndex(t);
    if (slices_.size() == 1 || i == slices_.size() - 1 || t <= slices_.front().market.t)
        return LognormalMixtureCalibrator::density(slices_[i].calibration.components, k);
    const auto& a = slices_[i];
    const auto& b = slices_[i + 1];
    const double w = (t - a.market.t) / (b.market.t - a.market.t);
    const double da = LognormalMixtureCalibrator::density(a.calibration.components, k);
    const double db = LognormalMixtureCalibrator::density(b.calibration.components, k);
    return (1.0 - w) * da + w * db;
}

double MultiMaturitySurface::impliedVol(double k, double t) const {
    if (t <= 0.0) return std::numeric_limits<double>::quiet_NaN();
    return LognormalMixtureCalibrator::impliedTotalVol(k, call(k, t)) / std::sqrt(t);
}

double MultiMaturitySurface::localVol(double k, double t) const {
    if (slices_.size() < 2 || t < slices_.front().market.t || t >= slices_.back().market.t)
        return std::numeric_limits<double>::quiet_NaN();
    const std::size_t i = leftIndex(t);
    if (i >= slices_.size() - 1) return std::numeric_limits<double>::quiet_NaN();
    const auto& a = slices_[i];
    const auto& b = slices_[i + 1];
    const double dt = b.market.t - a.market.t;
    const double ca = LognormalMixtureCalibrator::mixtureCall(a.calibration.components, k);
    const double cb = LognormalMixtureCalibrator::mixtureCall(b.calibration.components, k);
    const double dcdt = (cb - ca) / dt;
    const double d = density(k, t);
    const double denom = k * k * d;
    if (denom <= 1e-14 || dcdt < -1e-12) return std::numeric_limits<double>::quiet_NaN();
    return std::sqrt(std::max(0.0, 2.0 * dcdt / denom));
}

bool MultiMaturitySurface::calendarArbitrageFree(double kMin, double kMax, int nK, double tol) const {
    if (slices_.size() < 2) return true;
    for (std::size_t j = 0; j + 1 < slices_.size(); ++j) {
        for (int i = 0; i < nK; ++i) {
            const double k = kMin + (kMax - kMin) * i / static_cast<double>(nK - 1);
            const double a = LognormalMixtureCalibrator::mixtureCall(slices_[j].calibration.components, k);
            const double b = LognormalMixtureCalibrator::mixtureCall(slices_[j + 1].calibration.components, k);
            if (b + tol < a) return false;
        }
    }
    return true;
}

} // namespace lnmix
