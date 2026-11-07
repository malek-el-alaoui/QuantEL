#include "cubic_variance.hpp"
#include "lognormal_mixture.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace lnmix {

void NaturalCubicSpline::fit(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() != y.size() || x.size() < 2)
        throw std::runtime_error("NaturalCubicSpline needs at least two x/y points");

    const std::size_t n = x.size();
    for (std::size_t i = 1; i < n; ++i) {
        if (!(x[i] > x[i - 1]))
            throw std::runtime_error("NaturalCubicSpline x values must be strictly increasing");
    }

    x_ = x;
    a_.assign(n - 1, 0.0);
    b_.assign(n - 1, 0.0);
    c_.assign(n - 1, 0.0);
    d_.assign(n - 1, 0.0);

    std::vector<double> h(n - 1, 0.0);
    for (std::size_t i = 0; i + 1 < n; ++i) h[i] = x[i + 1] - x[i];

    std::vector<double> alpha(n, 0.0);
    for (std::size_t i = 1; i + 1 < n; ++i) {
        alpha[i] = 3.0 * (y[i + 1] - y[i]) / h[i]
                 - 3.0 * (y[i] - y[i - 1]) / h[i - 1];
    }

    std::vector<double> l(n, 1.0), mu(n, 0.0), z(n, 0.0), cc(n, 0.0);
    for (std::size_t i = 1; i + 1 < n; ++i) {
        l[i] = 2.0 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
        if (std::abs(l[i]) < 1e-18) throw std::runtime_error("Degenerate cubic spline system");
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }

    // Natural boundary conditions: second derivative is zero at both ends.
    cc[n - 1] = 0.0;
    for (std::size_t jj = n - 1; jj-- > 0;) {
        const std::size_t j = jj;
        cc[j] = z[j] - mu[j] * cc[j + 1];
        a_[j] = y[j];
        b_[j] = (y[j + 1] - y[j]) / h[j] - h[j] * (cc[j + 1] + 2.0 * cc[j]) / 3.0;
        c_[j] = cc[j];
        d_[j] = (cc[j + 1] - cc[j]) / (3.0 * h[j]);
    }
}

std::size_t NaturalCubicSpline::interval(double x) const {
    if (x_.empty()) throw std::runtime_error("Spline has not been fitted");
    if (x <= x_.front()) return 0;
    if (x >= x_.back()) return x_.size() - 2;
    const std::vector<double>::const_iterator it = std::upper_bound(x_.begin(), x_.end(), x);
    return static_cast<std::size_t>(it - x_.begin() - 1);
}

double NaturalCubicSpline::value(double x) const {
    const std::size_t i = interval(x);
    if (x < x_.front()) {
        return a_.front() + b_.front() * (x - x_.front());
    }
    if (x > x_.back()) {
        const std::size_t j = x_.size() - 2;
        const double h = x_.back() - x_[j];
        const double yEnd = a_[j] + b_[j] * h + c_[j] * h * h + d_[j] * h * h * h;
        const double slopeEnd = b_[j] + 2.0 * c_[j] * h + 3.0 * d_[j] * h * h;
        return yEnd + slopeEnd * (x - x_.back());
    }
    const double dx = x - x_[i];
    return a_[i] + b_[i] * dx + c_[i] * dx * dx + d_[i] * dx * dx * dx;
}

double NaturalCubicSpline::firstDerivative(double x) const {
    const std::size_t i = interval(x);
    if (x < x_.front()) return b_.front();
    if (x > x_.back()) {
        const std::size_t j = x_.size() - 2;
        const double h = x_.back() - x_[j];
        return b_[j] + 2.0 * c_[j] * h + 3.0 * d_[j] * h * h;
    }
    const double dx = x - x_[i];
    return b_[i] + 2.0 * c_[i] * dx + 3.0 * d_[i] * dx * dx;
}

double NaturalCubicSpline::secondDerivative(double x) const {
    const std::size_t i = interval(x);
    if (x < x_.front() || x > x_.back()) return 0.0; // linear extrapolation
    const double dx = x - x_[i];
    return 2.0 * c_[i] + 6.0 * d_[i] * dx;
}

double CubicVarianceSlice::rawTotalVariance(double k) const {
    if (!(k > 0.0)) return std::numeric_limits<double>::quiet_NaN();
    return spline.value(std::log(k));
}

double CubicVarianceSlice::totalVariance(double k) const {
    const double w = rawTotalVariance(k);
    if (!std::isfinite(w)) return w;
    return std::max(1e-14, w);
}

double CubicVarianceSlice::annualVol(double k) const {
    if (!(market.t > 0.0)) return std::numeric_limits<double>::quiet_NaN();
    return std::sqrt(totalVariance(k) / market.t);
}

double CubicVarianceSlice::call(double k) const {
    return LognormalMixtureCalibrator::bsCall(k, std::sqrt(totalVariance(k)));
}

double CubicVarianceSlice::butterflyFactor(double k) const {
    if (!(k > 0.0)) return std::numeric_limits<double>::quiet_NaN();
    const double x = std::log(k);
    const double w = spline.value(x);
    if (!(w > 1e-14) || !std::isfinite(w)) return std::numeric_limits<double>::quiet_NaN();
    const double wp = spline.firstDerivative(x);
    const double wpp = spline.secondDerivative(x);
    const double a = 1.0 - x * wp / (2.0 * w);
    return a * a - 0.25 * wp * wp * (1.0 / w + 0.25) + 0.5 * wpp;
}

double CubicVarianceSlice::density(double k) const {
    if (!(k > 0.0)) return std::numeric_limits<double>::quiet_NaN();
    const double x = std::log(k);
    const double w = spline.value(x);
    if (!(w > 1e-14) || !std::isfinite(w)) return std::numeric_limits<double>::quiet_NaN();
    const double sqrtW = std::sqrt(w);
    const double d2 = -x / sqrtW - 0.5 * sqrtW;
    const double g = butterflyFactor(k);
    return LognormalMixtureCalibrator::normPdf(d2) * g / (k * sqrtW);
}

CubicVarianceSurface::CubicVarianceSurface(MarketSurface market)
    : market_(std::move(market)) {
    std::sort(market_.slices.begin(), market_.slices.end(),
              [](const MarketSlice& a, const MarketSlice& b) { return a.t < b.t; });
}

void CubicVarianceSurface::calibrate() {
    slices_.clear();
    for (std::size_t s = 0; s < market_.slices.size(); ++s) {
        MarketSlice slice = market_.slices[s];
        std::sort(slice.points.begin(), slice.points.end(),
                  [](const MarketPoint& a, const MarketPoint& b) { return a.k < b.k; });
        if (slice.points.size() < 2) throw std::runtime_error("Cubic variance model needs at least two quotes per maturity");

        std::vector<double> x, w;
        x.reserve(slice.points.size());
        w.reserve(slice.points.size());
        for (std::size_t i = 0; i < slice.points.size(); ++i) {
            if (!(slice.points[i].k > 0.0) || !(slice.points[i].marketVol > 0.0) || !(slice.t > 0.0))
                throw std::runtime_error("Invalid quote for cubic variance model");
            x.push_back(std::log(slice.points[i].k));
            w.push_back(slice.points[i].marketVol * slice.points[i].marketVol * slice.t);
        }
        CubicVarianceSlice out;
        out.market = slice;
        out.spline.fit(x, w);
        slices_.push_back(out);
    }
}

std::size_t CubicVarianceSurface::leftIndex(double t) const {
    if (slices_.empty()) throw std::runtime_error("Cubic variance surface is not calibrated");
    if (slices_.size() == 1 || t <= slices_.front().market.t) return 0;
    if (t >= slices_.back().market.t) return slices_.size() - 1;
    for (std::size_t i = 0; i + 1 < slices_.size(); ++i)
        if (t >= slices_[i].market.t && t <= slices_[i + 1].market.t) return i;
    return slices_.size() - 1;
}

double CubicVarianceSurface::call(double k, double t) const {
    const std::size_t i = leftIndex(t);
    if (slices_.size() == 1 || i == slices_.size() - 1 || t <= slices_.front().market.t)
        return slices_[i].call(k);
    const CubicVarianceSlice& a = slices_[i];
    const CubicVarianceSlice& b = slices_[i + 1];
    const double q = (t - a.market.t) / (b.market.t - a.market.t);
    return (1.0 - q) * a.call(k) + q * b.call(k);
}

double CubicVarianceSurface::density(double k, double t) const {
    const std::size_t i = leftIndex(t);
    if (slices_.size() == 1 || i == slices_.size() - 1 || t <= slices_.front().market.t)
        return slices_[i].density(k);
    const CubicVarianceSlice& a = slices_[i];
    const CubicVarianceSlice& b = slices_[i + 1];
    const double q = (t - a.market.t) / (b.market.t - a.market.t);
    const double da = a.density(k);
    const double db = b.density(k);
    if (!std::isfinite(da) || !std::isfinite(db)) return std::numeric_limits<double>::quiet_NaN();
    return (1.0 - q) * da + q * db;
}

double CubicVarianceSurface::impliedVol(double k, double t) const {
    if (!(t > 0.0)) return std::numeric_limits<double>::quiet_NaN();
    return LognormalMixtureCalibrator::impliedTotalVol(k, call(k, t)) / std::sqrt(t);
}

double CubicVarianceSurface::localVol(double k, double t) const {
    if (slices_.size() < 2 || t < slices_.front().market.t || t >= slices_.back().market.t)
        return std::numeric_limits<double>::quiet_NaN();
    const std::size_t i = leftIndex(t);
    if (i >= slices_.size() - 1) return std::numeric_limits<double>::quiet_NaN();
    const CubicVarianceSlice& a = slices_[i];
    const CubicVarianceSlice& b = slices_[i + 1];
    const double dcdt = (b.call(k) - a.call(k)) / (b.market.t - a.market.t);
    const double den = density(k, t);
    const double denom = k * k * den;
    if (!std::isfinite(denom) || denom <= 1e-14 || dcdt < -1e-12)
        return std::numeric_limits<double>::quiet_NaN();
    return std::sqrt(std::max(0.0, 2.0 * dcdt / denom));
}

bool CubicVarianceSurface::calendarArbitrageFree(double kMin, double kMax, int nK, double tol) const {
    if (slices_.size() < 2) return true;
    for (std::size_t j = 0; j + 1 < slices_.size(); ++j) {
        for (int i = 0; i < nK; ++i) {
            const double k = kMin + (kMax - kMin) * i / static_cast<double>(nK - 1);
            const double a = slices_[j].call(k);
            const double b = slices_[j + 1].call(k);
            if (b + tol < a) return false;
        }
    }
    return true;
}

} // namespace lnmix
