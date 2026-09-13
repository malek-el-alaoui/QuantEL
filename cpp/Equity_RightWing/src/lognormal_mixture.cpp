#include "lognormal_mixture.hpp"

#include <cstddef>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace lnmix {
namespace {
constexpr double PI = 3.141592653589793238462643383279502884;

double intrinsic(double center, double k) {
    return std::max(center - k, 0.0);
}

} // namespace

double LognormalMixtureCalibrator::normCdf(double x) {
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

double LognormalMixtureCalibrator::normPdf(double x) {
    return std::exp(-0.5 * x * x) / std::sqrt(2.0 * PI);
}

// Peter J. Acklam inverse normal CDF approximation.
double LognormalMixtureCalibrator::invNormCdf(double p) {
    if (!(p > 0.0 && p < 1.0)) {
        if (p == 0.0) return -std::numeric_limits<double>::infinity();
        if (p == 1.0) return std::numeric_limits<double>::infinity();
        throw std::runtime_error("invNormCdf requires 0 < p < 1");
    }

    const double a1 = -3.969683028665376e+01;
    const double a2 =  2.209460984245205e+02;
    const double a3 = -2.759285104469687e+02;
    const double a4 =  1.383577518672690e+02;
    const double a5 = -3.066479806614716e+01;
    const double a6 =  2.506628277459239e+00;
    const double b1 = -5.447609879822406e+01;
    const double b2 =  1.615858368580409e+02;
    const double b3 = -1.556989798598866e+02;
    const double b4 =  6.680131188771972e+01;
    const double b5 = -1.328068155288572e+01;
    const double c1 = -7.784894002430293e-03;
    const double c2 = -3.223964580411365e-01;
    const double c3 = -2.400758277161838e+00;
    const double c4 = -2.549732539343734e+00;
    const double c5 =  4.374664141464968e+00;
    const double c6 =  2.938163982698783e+00;
    const double d1 =  7.784695709041462e-03;
    const double d2 =  3.224671290700398e-01;
    const double d3 =  2.445134137142996e+00;
    const double d4 =  3.754408661907416e+00;
    const double pLow = 0.02425;
    const double pHigh = 1.0 - pLow;

    double x;
    if (p < pLow) {
        const double q = std::sqrt(-2.0 * std::log(p));
        x = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) /
            ((((d1*q+d2)*q+d3)*q+d4)*q+1.0);
    } else if (p <= pHigh) {
        const double q = p - 0.5;
        const double r = q * q;
        x = (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q /
            (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1.0);
    } else {
        const double q = std::sqrt(-2.0 * std::log(1.0 - p));
        x = -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) /
             ((((d1*q+d2)*q+d3)*q+d4)*q+1.0);
    }

    // One Halley refinement step.
    const double e = normCdf(x) - p;
    const double u = e / normPdf(x);
    x -= u / (1.0 + 0.5 * x * u);
    return x;
}

double LognormalMixtureCalibrator::bsCall(double k, double v) {
    if (v <= 1e-14) return std::max(1.0 - k, 0.0);
    const double d1 = -std::log(k) / v + 0.5 * v;
    const double d2 = d1 - v;
    return normCdf(d1) - k * normCdf(d2);
}

double LognormalMixtureCalibrator::bsVega(double k, double v) {
    if (v <= 1e-12) {
        return std::abs(k - 1.0) < 1e-12 ? 1.0 / std::sqrt(2.0 * PI) : 0.0;
    }
    const double d1 = -std::log(k) / v + 0.5 * v;
    return normPdf(d1);
}

double LognormalMixtureCalibrator::shiftedCall(double center, double k, double sigma) {
    if (sigma <= 1e-14) return intrinsic(center, k);
    const double d1 = -std::log(k / center) / sigma + 0.5 * sigma;
    const double d2 = d1 - sigma;
    return center * normCdf(d1) - k * normCdf(d2);
}

double LognormalMixtureCalibrator::shiftedVega(double center, double k, double sigma) {
    if (sigma <= 1e-12) {
        return std::abs(center - k) < 1e-12 ? center / std::sqrt(2.0 * PI) : 0.0;
    }
    const double d1 = -std::log(k / center) / sigma + 0.5 * sigma;
    return center * normPdf(d1);
}

double LognormalMixtureCalibrator::mixtureCall(const std::vector<Component>& cs, double k) {
    double x = 0.0;
    for (const auto& q : cs) x += q.weight * shiftedCall(q.center, k, q.sigma);
    return x;
}

double LognormalMixtureCalibrator::density(const std::vector<Component>& cs, double k) {
    if (k <= 0.0) return 0.0;
    double x = 0.0;
    for (const auto& q : cs) {
        if (q.sigma <= 0.0 || q.weight <= 0.0) continue;
        const double z = std::log(k / q.center) / q.sigma;
        x += q.weight * std::exp(-0.5 * z * z) /
             (k * q.sigma * std::sqrt(2.0 * PI));
    }
    return x;
}


double LognormalMixtureCalibrator::wingSigmaMultiplier(const CalibrationResult& r,
                                                        double k) {
    if (!(k > 0.0)) return 0.0;
    const double power = r.wingPower > 0.0 ? r.wingPower : 3.0;

    // All kernels are unchanged on the observed market interval. Outside it,
    // every kernel receives the same smooth side-specific multiplicative
    // decay.  Using a common multiplier preserves the calibrated relative
    // kernel widths while allowing the whole mixture to contract in the wing.
    if (r.marketKMin > 0.0 && k < r.marketKMin && r.leftWingAlpha > 0.0) {
        const double d = std::log(r.marketKMin / k);
        return std::exp(-r.leftWingAlpha * std::pow(d, power));
    }
    if (r.marketKMax > 0.0 && k > r.marketKMax && r.rightWingAlpha > 0.0) {
        const double d = std::log(k / r.marketKMax);
        return std::exp(-r.rightWingAlpha * std::pow(d, power));
    }
    return 1.0;
}

double LognormalMixtureCalibrator::wingSigma(const CalibrationResult& r,
                                              std::size_t componentIndex,
                                              double k) {
    if (r.components.empty() || componentIndex >= r.components.size() || !(k > 0.0))
        return 0.0;

    // Component-weighted wing decay.  All calibrated sigmas remain EXACTLY
    // unchanged on the observed market interval.  Outside it, components
    // whose centers lie farther from the active wing are damped more strongly.
    // This removes cross-wing leakage while keeping only ONE alpha per side.
    const double base = r.components[componentIndex].sigma;
    if (!(r.marketKMin > 0.0) || !(r.marketKMax > r.marketKMin))
        return base;

    const double x = std::log(k);
    const double xL = std::log(r.marketKMin);
    const double xR = std::log(r.marketKMax);
    const double power = r.wingPower > 0.0 ? r.wingPower : 3.0;

    double multiplier = 1.0;
    if (x < xL && r.leftWingAlpha > 0.0) {
        const double d = xL - x;
        const double sigmaRef = std::max(1e-8, r.atmTotalVol > 0.0 ? r.atmTotalVol : 0.5 * (r.market.front().totalVol + r.market.back().totalVol));
        const double volRatio = base / sigmaRef;
        const double componentFactor = std::max(0.01, std::min(9.0, volRatio * volRatio));
        multiplier = std::exp(-r.leftWingAlpha * componentFactor * std::pow(d, power));
    } else if (x > xR && r.rightWingAlpha > 0.0) {
        const double d = x - xR;
        const double sigmaRef = std::max(1e-8, r.atmTotalVol > 0.0 ? r.atmTotalVol : 0.5 * (r.market.front().totalVol + r.market.back().totalVol));
        const double volRatio = base / sigmaRef;
        const double componentFactor = std::max(0.01, std::min(9.0, volRatio * volRatio));
        multiplier = std::exp(-r.rightWingAlpha * componentFactor * std::pow(d, power));
    }

    return std::max(1e-10, base * multiplier);
}

double LognormalMixtureCalibrator::wingedMixtureCall(const CalibrationResult& r, double k) {
    double x = 0.0;
    for (std::size_t i = 0; i < r.components.size(); ++i) {
        const Component& q = r.components[i];
        if (q.weight <= 0.0) continue;
        x += q.weight * shiftedCall(q.center, k, wingSigma(r, i, k));
    }
    return x;
}

double LognormalMixtureCalibrator::wingedImpliedTotalVariance(const CalibrationResult& r, double k) {
    const double v = impliedTotalVol(k, wingedMixtureCall(r, k));
    return v * v;
}

double LognormalMixtureCalibrator::wingedButterflyFactor(const CalibrationResult& r, double k) {
    if (!(k > 0.0)) return std::numeric_limits<double>::quiet_NaN();
    const double x = std::log(k);
    const double h = 5e-4;
    const double wm = wingedImpliedTotalVariance(r, std::exp(x - h));
    const double w0 = wingedImpliedTotalVariance(r, k);
    const double wpv = wingedImpliedTotalVariance(r, std::exp(x + h));
    if (!(w0 > 1e-12) || !std::isfinite(wm) || !std::isfinite(w0) || !std::isfinite(wpv))
        return std::numeric_limits<double>::quiet_NaN();
    const double w1 = (wpv - wm) / (2.0 * h);
    const double w2 = (wpv - 2.0 * w0 + wm) / (h * h);
    const double a = 1.0 - x * w1 / (2.0 * w0);
    return a * a - 0.25 * w1 * w1 * (1.0 / w0 + 0.25) + 0.5 * w2;
}

double LognormalMixtureCalibrator::wingedDensity(const CalibrationResult& r, double k) {
    if (!(k > 0.0)) return std::numeric_limits<double>::quiet_NaN();
    const double w = wingedImpliedTotalVariance(r, k);
    if (!(w > 1e-12) || !std::isfinite(w)) {
        // In a numerically exhausted far tail, the density tends to zero.
        return 0.0;
    }
    const double g = wingedButterflyFactor(r, k);
    if (!std::isfinite(g)) return std::numeric_limits<double>::quiet_NaN();
    const double x = std::log(k);
    const double sqrtW = std::sqrt(w);
    const double d2 = -x / sqrtW - 0.5 * sqrtW;
    return normPdf(d2) * g / (k * sqrtW);
}

WingArbitrageMetrics LognormalMixtureCalibrator::wingArbitrageMetrics(
    const CalibrationResult& r,
    double xMin,
    double xMax,
    int points) {

    WingArbitrageMetrics out;
    if (!(xMax > xMin) || points < 11) return out;

    std::vector<double> x(points), k(points), call(points), w(points), timeValue(points);
    for (int i = 0; i < points; ++i) {
        x[i] = xMin + (xMax - xMin) * i / static_cast<double>(points - 1);
        k[i] = std::exp(x[i]);
        call[i] = wingedMixtureCall(r, k[i]);
        const double intrinsicValue = std::max(1.0 - k[i], 0.0);
        timeValue[i] = call[i] - intrinsicValue;
        const double v = impliedTotalVol(k[i], call[i]);
        w[i] = v * v;
    }

    // Direct normalized-call conditions are the definitive strike-static test.
    out.directCallArbitrageFree = discreteArbitrageFree(k, call, 2e-8);

    const double h = (xMax - xMin) / static_cast<double>(points - 1);
    double minG = std::numeric_limits<double>::infinity();
    double minD = std::numeric_limits<double>::infinity();
    int negative = 0;
    int invalid = 0;

    // Ignore points where option time value is numerically exhausted.  There
    // the Black inverse is ill-conditioned even though the call itself can be
    // perfectly monotone and convex.
    for (int i = 1; i + 1 < points; ++i) {
        if (!(w[i] > 1e-8) || !(timeValue[i] > 1e-9) ||
            !std::isfinite(w[i - 1]) || !std::isfinite(w[i]) || !std::isfinite(w[i + 1])) {
            ++invalid;
            continue;
        }
        const double w1 = (w[i + 1] - w[i - 1]) / (2.0 * h);
        const double w2 = (w[i + 1] - 2.0 * w[i] + w[i - 1]) / (h * h);
        const double a = 1.0 - x[i] * w1 / (2.0 * w[i]);
        const double g = a * a - 0.25 * w1 * w1 * (1.0 / w[i] + 0.25) + 0.5 * w2;
        if (!std::isfinite(g)) {
            ++invalid;
            continue;
        }
        minG = std::min(minG, g);
        const double sqrtW = std::sqrt(w[i]);
        const double d2 = -x[i] / sqrtW - 0.5 * sqrtW;
        const double den = normPdf(d2) * g / (k[i] * sqrtW);
        if (std::isfinite(den)) {
            minD = std::min(minD, den);
            if (den < -1e-8) ++negative;
        }
    }

    out.minButterflyG = std::isfinite(minG) ? minG : std::numeric_limits<double>::quiet_NaN();
    out.minDensity = std::isfinite(minD) ? minD : std::numeric_limits<double>::quiet_NaN();
    out.negativeDensityPoints = negative;
    out.invalidVariancePoints = invalid;
    return out;
}

void LognormalMixtureCalibrator::configureWingDecay(CalibrationResult& r) {
    if (r.components.empty() || r.market.empty()) return;

    r.marketKMin = r.market.front().k;
    r.marketKMax = r.market.front().k;
    std::size_t right = 0;
    for (std::size_t i = 1; i < r.market.size(); ++i) {
        if (r.market[i].k < r.marketKMin) { r.marketKMin = r.market[i].k; }
        if (r.market[i].k > r.marketKMax) { r.marketKMax = r.market[i].k; right = i; }
    }
    r.leftWingBaseSigma = r.components.front().sigma;
    r.rightWingBaseSigma = r.components.back().sigma;
    r.wingPower = 2.0;

    // Equity-only policy: preserve the left wing exactly as calibrated and
    // control only the problematic right wing.  A quadratic attachment is C1
    // at x_R and becomes effective much faster than the previous cubic onset.
    r.leftWingAlpha = 0.0;
    r.leftWingCriticalAlpha = 0.0;

    const double deltaSafety = 0.001; // ~0.1-delta safety tail
    const double d1RightSafety = invNormCdf(deltaSafety);
    const double vRight = std::max(1e-6, r.market[right].totalVol);
    const double xDeltaRightSafety = vRight * (0.5 * vRight - d1RightSafety);
    const double xR = std::log(r.marketKMax);
    const double xL = std::log(r.marketKMin);
    r.wingSafetyXMin = xL;
    r.wingSafetyXMax = std::max(xR + 0.04, xDeltaRightSafety);

    // Wing diagnostic/shape target is measured to about 0.5-delta, which is
    // far enough to reveal the extrapolation without relying on numerically
    // exhausted option time value.
    r.wingTargetDelta = 0.005;
    const double d1RightTarget = invNormCdf(r.wingTargetDelta);
    const double xRightTarget = std::max(xR + 0.015,
        vRight * (0.5 * vRight - d1RightTarget));
    const double xEvalMax = std::min(r.wingSafetyXMax, xRightTarget);

    const double alphaCap = 2.0e6;
    const double gTol = -5e-4; // FD tolerance; call convexity is definitive.

    auto feasible = [&](double aR) {
        CalibrationResult trial = r;
        trial.leftWingAlpha = 0.0;
        trial.rightWingAlpha = aR;
        const WingArbitrageMetrics m = wingArbitrageMetrics(
            trial, trial.wingSafetyXMin, trial.wingSafetyXMax, 1801);
        return m.directCallArbitrageFree &&
               (m.negativeDensityPoints == 0) &&
               (!std::isfinite(m.minButterflyG) || m.minButterflyG >= gTol);
    };

    // Largest statically admissible right-wing decay parameter.
    double lo = 0.0, hi = 1.0;
    while (hi < alphaCap && feasible(hi)) {
        lo = hi;
        hi *= 2.0;
    }
    if (hi >= alphaCap && feasible(alphaCap)) {
        r.rightWingCriticalAlpha = alphaCap;
    } else {
        hi = std::min(hi, alphaCap);
        for (int it = 0; it < 48; ++it) {
            const double mid = 0.5 * (lo + hi);
            if (feasible(mid)) lo = mid;
            else hi = mid;
        }
        r.rightWingCriticalAlpha = lo;
    }

    // Diagnose aggregate implied-vol shape on the right wing.  Because the
    // quadratic sigma attachment has zero first derivative exactly at x_R,
    // the very first infinitesimal segment may inherit a positive slope from
    // the calibrated mixture.  We therefore require non-increase after a tiny
    // attachment buffer (2% of the tested log-wing span, at least 2e-4).
    auto shapeMetrics = [&](double aR,
                            double& boundaryIv,
                            double& endIv,
                            double& maxOvershootBp,
                            double& maxStepRiseBp,
                            bool& nonIncreasing) {
        CalibrationResult trial = r;
        trial.leftWingAlpha = 0.0;
        trial.rightWingAlpha = aR;
        const int nGrid = 801;
        const double span = std::max(1e-6, xEvalMax - xR);
        const double buffer = std::max(2e-4, 0.02 * span);
        const double xStartMono = xR + buffer;
        boundaryIv = impliedTotalVol(r.marketKMax,
                    wingedMixtureCall(trial, r.marketKMax));
        double prevIv = boundaryIv;
        double maxIv = boundaryIv;
        maxStepRiseBp = 0.0;
        nonIncreasing = true;
        bool haveMonoPrev = false;
        double monoPrev = 0.0;
        for (int i = 1; i < nGrid; ++i) {
            const double x = xR + span * i / static_cast<double>(nGrid - 1);
            const double k = std::exp(x);
            const double iv = impliedTotalVol(k, wingedMixtureCall(trial, k));
            if (std::isfinite(iv)) {
                maxIv = std::max(maxIv, iv);
                const double rise = (iv - prevIv) * 10000.0;
                maxStepRiseBp = std::max(maxStepRiseBp, rise);
                prevIv = iv;
                if (x >= xStartMono) {
                    if (haveMonoPrev && iv > monoPrev + 1e-7) nonIncreasing = false;
                    monoPrev = iv;
                    haveMonoPrev = true;
                }
            }
        }
        endIv = prevIv;
        maxOvershootBp = std::max(0.0, (maxIv - boundaryIv) * 10000.0);
    };

    // Search for the SMALLEST alpha that achieves a non-increasing aggregate
    // right wing after the tiny C1 attachment buffer.  If static arbitrage
    // prevents that, use 95% of the critical alpha and report the residual
    // overshoot explicitly rather than hiding the limitation.
    const double usableCritical = 0.95 * r.rightWingCriticalAlpha;
    double bIv=0.0, eIv=0.0, over=0.0, rise=0.0;
    bool mono=false;
    shapeMetrics(usableCritical, bIv, eIv, over, rise, mono);
    r.rightWingMonotoneTargetFeasible = mono && feasible(usableCritical);

    if (r.rightWingMonotoneTargetFeasible) {
        double aLo = 0.0, aHi = usableCritical;
        for (int it = 0; it < 48; ++it) {
            const double mid = 0.5 * (aLo + aHi);
            double tb=0.0, te=0.0, to=0.0, tr=0.0; bool tm=false;
            shapeMetrics(mid, tb, te, to, tr, tm);
            if (tm && feasible(mid)) aHi = mid;
            else aLo = mid;
        }
        r.rightWingAlpha = aHi;
    } else {
        r.rightWingAlpha = usableCritical;
    }

    r.wingSafetyScale = 1.0;
    shapeMetrics(r.rightWingAlpha,
                 r.rightWingBoundaryIv,
                 r.rightWingEndIv,
                 r.rightWingMaxOvershootBp,
                 r.rightWingMaxStepRiseBp,
                 r.rightWingNonIncreasing);

    // A convenient diagnostic ratio: multiplier of an ATM-width kernel at the
    // end of the tested right-wing target.  Actual components scale according
    // to (sigma_j/sigma_ATM)^2 in wingSigma().
    const double dTarget = std::max(0.0, xEvalMax - xR);
    r.wingTargetSigmaRatio = std::exp(-r.rightWingAlpha * dTarget * dTarget);

    const WingArbitrageMetrics finalMetrics = wingArbitrageMetrics(
        r, r.wingSafetyXMin, r.wingSafetyXMax, 2201);
    r.wingDirectCallArbitrageFree = finalMetrics.directCallArbitrageFree;
    r.wingMinButterflyG = finalMetrics.minButterflyG;
    r.wingMinDensity = finalMetrics.minDensity;
    r.wingNegativeDensityPoints = finalMetrics.negativeDensityPoints;
    r.wingInvalidVariancePoints = finalMetrics.invalidVariancePoints;
}

double LognormalMixtureCalibrator::impliedTotalVol(double k, double call) {
    const double intrinsicValue = std::max(1.0 - k, 0.0);
    call = std::max(call, intrinsicValue);
    call = std::min(call, 1.0 - 1e-14);

    if (call <= intrinsicValue + 1e-14) return 0.0;

    // Adaptive upper bracket.  This matters for the deliberately remote
    // artificial quote used to make its boundary mass very small.
    double lo = 1e-10;
    double hi = 1.0;
    while (hi < 32.0 && bsCall(k, hi) < call) hi *= 2.0;
    if (bsCall(k, hi) < call) return hi;

    for (int i = 0; i < 180; ++i) {
        const double m = 0.5 * (lo + hi);
        if (bsCall(k, m) < call) lo = m;
        else hi = m;
    }
    return 0.5 * (lo + hi);
}

bool LognormalMixtureCalibrator::discreteArbitrageFree(const std::vector<double>& k,
                                                        const std::vector<double>& c,
                                                        double eps) {
    if (k.size() != c.size() || k.size() < 2) return false;
    double prevSlope = -1.0;
    for (std::size_t i = 0; i < k.size(); ++i) {
        if (!(k[i] > 0.0)) return false;
        if (i > 0 && !(k[i] > k[i-1])) return false;
        if (c[i] < std::max(1.0 - k[i], 0.0) - eps || c[i] > 1.0 + eps) return false;
        if (i > 0) {
            const double slope = (c[i] - c[i-1]) / (k[i] - k[i-1]);
            if (slope < prevSlope - eps || slope > eps) return false;
            prevSlope = slope;
        }
    }
    return true;
}

CalibrationResult LognormalMixtureCalibrator::calibrateBalancedBoundary(
    const std::vector<std::pair<double, double>>& kv,
    double eta,
    int maxIter,
    double tol,
    double artificialQuoteCalibrationWeight) const {

    if (kv.size() < 3) throw std::runtime_error("Need at least 3 smile quotes");

    std::vector<std::pair<double, double>> q = kv;
    std::sort(q.begin(), q.end(),
              [](const std::pair<double, double>& a,
                 const std::pair<double, double>& b) {
                  return a.first < b.first;
              });

    const int nMarket = static_cast<int>(q.size());
    std::vector<double> km(nMarket), cm(nMarket), vm(nMarket);
    for (int i = 0; i < nMarket; ++i) {
        km[i] = q[i].first;
        vm[i] = q[i].second;
        if (km[i] <= 0.0 || vm[i] <= 0.0)
            throw std::runtime_error("k and totalVol must be positive");
        cm[i] = bsCall(km[i], vm[i]);
    }

    const bool marketArb = discreteArbitrageFree(km, cm);
    if (!marketArb) {
        throw std::runtime_error(
            "Input call quotes are not discretely arbitrage-free. See diagnostics/output.");
    }

    const double p1 = cm.front() + km.front() - 1.0;
    const double cN = cm.back();
    const double compareTol = 1e-13;
    const double ARTIFICIAL_QUOTE_CALIBRATION_WEIGHT =
        std::max(0.0, std::min(1.0, artificialQuoteCalibrationWeight));

    if (p1 < -1e-12 || cN < -1e-12)
        throw std::runtime_error("Invalid boundary option prices p1 or cN");

    std::vector<double> k;
    std::vector<double> c;
    std::vector<double> fitVol;
    k.reserve(nMarket + 1);
    c.reserve(nMarket + 1);
    fitVol.reserve(nMarket + 1);

    std::string artificialSide = "none";
    double artificialK = 0.0;
    double artificialCall = 0.0;
    double theoreticalMinArtificialWeight = 0.0;
    double boundarySlope = 0.0;

    // ------------------------------------------------------------------
    // BALANCED ONE-POINT CONSTRUCTION.
    //
    // The boundary quote is still the one-point construction of the paper:
    //   left  if p1 > cN : p0 = cN,
    //   right if cN > p1 : cN+1 = p1.
    //
    // What changes is HOW its strike is chosen.  The total butterfly mass
    // carried by the boundary knot and the neighboring market knot is fixed.
    // We choose the artificial strike so that this fixed mass is split as
    // evenly as possible.  This simultaneously minimizes
    //      max(w_boundary, w_neighbor)
    // and
    //      w_boundary^2 + w_neighbor^2,
    // which is the natural local choice for suppressing an endpoint hump.
    // ------------------------------------------------------------------

    if (p1 > cN + compareTol) {
        // LEFT case.  Let d1 be the slope from the artificial quote to the
        // first market quote, and d2 the first market-market slope.
        // Then
        //      w0 = 1 + d1,
        //      w1 = d2 - d1,
        // and w0+w1 = 1+d2 is fixed.  The minimax/L2 optimum is w0=w1,
        // hence d1* = (d2-1)/2.
        const double d2 = (cm[1] - cm[0]) / (km[1] - km[0]);
        if (!(d2 > -1.0 + 1e-14 && d2 < 1e-14))
            throw std::runtime_error("Unexpected first market slope in left boundary construction");

        const double A = p1 - cN;
        const double k0Lower = std::max(cN, 1e-14); // c0<1 and k0>0
        const double k0Upper = km.front() - A / (1.0 + d2); // d1<=d2
        if (!(k0Upper > k0Lower))
            throw std::runtime_error("Cannot place balanced left artificial quote while preserving convexity");

        const double guard = std::max(1e-12, 1e-12 * std::max(1.0, km.front()));
        const double lowerGuard = k0Lower + guard;
        const double upperGuard = k0Upper - guard;
        if (!(upperGuard > lowerGuard))
            throw std::runtime_error("Left artificial quote feasible interval collapsed numerically");

        const double dTarget = 0.5 * (d2 - 1.0);
        const double dLower = -1.0 + A / (km.front() - lowerGuard);
        const double dUpper = -1.0 + A / (km.front() - upperGuard);
        const double dChosen = std::max(dLower, std::min(dUpper, dTarget));

        const double k0 = km.front() - A / (1.0 + dChosen);
        const double c0 = 1.0 - k0 + cN;

        artificialSide = "left";
        artificialK = k0;
        artificialCall = c0;
        boundarySlope = dChosen;
        theoreticalMinArtificialWeight = A / (km.front() - k0Lower);

        k.push_back(k0);
        c.push_back(c0);
        fitVol.push_back(impliedTotalVol(k0, c0));
        k.insert(k.end(), km.begin(), km.end());
        c.insert(c.end(), cm.begin(), cm.end());
        fitVol.insert(fitVol.end(), vm.begin(), vm.end());
    }
    else if (cN > p1 + compareTol) {
        // RIGHT case.  Let dN be the last market-market slope and dR the
        // slope from the last market quote to the artificial quote.  Then
        //      wN   = dR - dN,
        //      wN+1 = -dR,
        // and their sum -dN is fixed.  The optimum is dR=dN/2.
        const double dN =
            (cm[nMarket - 1] - cm[nMarket - 2]) /
            (km[nMarket - 1] - km[nMarket - 2]);
        if (!(dN < -1e-14))
            throw std::runtime_error("Last market slope too close to zero for right boundary construction");

        const double dR = 0.5 * dN;
        const double kR = km.back() + (p1 - cN) / dR;
        const double cR = p1;
        if (!(kR > km.back()) || !std::isfinite(kR))
            throw std::runtime_error("Invalid balanced right artificial quote location");

        artificialSide = "right";
        artificialK = kR;
        artificialCall = cR;
        boundarySlope = dR;
        theoreticalMinArtificialWeight = 0.0;

        k.insert(k.end(), km.begin(), km.end());
        c.insert(c.end(), cm.begin(), cm.end());
        fitVol.insert(fitVol.end(), vm.begin(), vm.end());
        k.push_back(kR);
        c.push_back(cR);
        fitVol.push_back(impliedTotalVol(kR, cR));
    }
    else {
        k = km;
        c = cm;
        fitVol = vm;
    }

    if (!discreteArbitrageFree(k, c))
        throw std::runtime_error("Balanced one-point smile is not arbitrage-free");

    // Mixture weights = jumps of piecewise-linear call slopes.
    const int nKnots = static_cast<int>(k.size());
    std::vector<double> slope(nKnots + 1, 0.0);
    slope[0] = -1.0;
    for (int i = 1; i < nKnots; ++i)
        slope[i] = (c[i] - c[i - 1]) / (k[i] - k[i - 1]);
    slope[nKnots] = 0.0;

    std::vector<Component> comp;
    comp.reserve(nKnots);
    for (int i = 0; i < nKnots; ++i) {
        const double weight = slope[i + 1] - slope[i];
        if (weight < -1e-9)
            throw std::runtime_error("Negative mixture weight after balanced one-point construction");
        if (weight > 1e-12) comp.push_back({k[i], weight, 0.0});
    }
    const int n = static_cast<int>(comp.size());
    if (n == 0) throw std::runtime_error("No positive mixture components");

    double sumW = 0.0, firstMoment = 0.0;
    for (const auto& z : comp) {
        sumW += z.weight;
        firstMoment += z.weight * z.center;
    }
    if (std::abs(sumW - 1.0) > 1e-7 || std::abs(firstMoment - 1.0) > 1e-7)
        throw std::runtime_error("Internal mixture moment constraint failure");

    int atmIdx = 0;
    for (int i = 1; i < nMarket; ++i)
        if (std::abs(km[i] - 1.0) < std::abs(km[atmIdx] - 1.0)) atmIdx = i;

    // User-selectable lower bound eta*sigma_ATM (e.g. eta=0.10 or 0.02).
    const double effectiveEta = std::max(1e-6, eta);
    const double atmTotalVol = vm[atmIdx];
    const double sigmaMin = effectiveEta * atmTotalVol;
    const double floorSigma = std::max(1e-8, sigmaMin * (1.0 + 1e-12));

    // Equity-tail stabilization: prevent the calibration from using extremely
    // broad kernels to fit a modest market smile.  Such kernels dominate the
    // remote right tail and can force an unrealistic ascending equity wing.
    // The bound is on TOTAL kernel volatility, consistently with vm[].
    const double sigmaCap = std::max(floorSigma * (1.0 + 1e-10),
                                     3.0 * atmTotalVol * (1.0 - 1e-12));

    const double smoothLevel = std::max(1e-12, std::min(p1, cN));
    for (auto& z : comp) {
        const double firstGuess = smoothLevel * std::sqrt(2.0 * PI) /
                                  std::max(1e-12, z.weight * z.center);
        z.sigma = std::min(sigmaCap, std::max(floorSigma, firstGuess));
    }

    // Fit all market quotes and give the single synthetic boundary quote an
    // almost negligible objective weight.
    const int nFit = static_cast<int>(k.size());
    std::vector<double> targetVega(nFit, 0.0);
    std::vector<double> fitObjectiveWeight(nFit, 1.0);
    for (int i = 0; i < nFit; ++i)
        targetVega[i] = std::max(1e-7, bsVega(k[i], fitVol[i]));

    int artificialFitIndex = -1;
    if (artificialSide == "left") artificialFitIndex = 0;
    else if (artificialSide == "right") artificialFitIndex = nFit - 1;
    if (artificialFitIndex >= 0)
        fitObjectiveWeight[artificialFitIndex] = ARTIFICIAL_QUOTE_CALIBRATION_WEIGHT;

    auto objective = [&](const std::vector<Component>& cc) {
        double s = 0.0;
        for (int i = 0; i < nFit; ++i) {
            const double r = (mixtureCall(cc, k[i]) - c[i]) / targetVega[i];
            s += fitObjectiveWeight[i] * r * r;
        }
        return s;
    };

    double obj = objective(comp);
    double lambda = 1e-5;
    int it = 0;
    for (; it < maxIter && obj > tol; ++it) {
        std::vector<double> residual(nFit, 0.0);
        std::vector<std::vector<double>> J(nFit, std::vector<double>(n, 0.0));

        for (int i = 0; i < nFit; ++i) {
            const double rowScale = std::sqrt(fitObjectiveWeight[i]);
            residual[i] = rowScale *
                (mixtureCall(comp, k[i]) - c[i]) / targetVega[i];
            for (int j = 0; j < n; ++j) {
                J[i][j] = rowScale * comp[j].weight *
                    shiftedVega(comp[j].center, k[i], comp[j].sigma) /
                    targetVega[i];
            }
        }

        std::vector<std::vector<double>> A(n, std::vector<double>(n, 0.0));
        std::vector<double> b(n, 0.0);
        for (int a = 0; a < n; ++a) {
            for (int i = 0; i < nFit; ++i) {
                b[a] -= J[i][a] * residual[i];
                for (int bb = 0; bb < n; ++bb) A[a][bb] += J[i][a] * J[i][bb];
            }
            A[a][a] += lambda;
        }

        std::vector<double> dx(n, 0.0), lower(n, 0.0), upper(n, 0.0);
        for (int j = 0; j < n; ++j) {
            lower[j] = floorSigma - comp[j].sigma;
            upper[j] = sigmaCap - comp[j].sigma;
        }

        for (int sweep = 0; sweep < 2000; ++sweep) {
            double maxChange = 0.0;
            for (int j = 0; j < n; ++j) {
                double rhs = b[j];
                for (int l = 0; l < n; ++l)
                    if (l != j) rhs -= A[j][l] * dx[l];
                double candidate = rhs / std::max(1e-18, A[j][j]);
                candidate = std::max(lower[j], std::min(upper[j], candidate));
                maxChange = std::max(maxChange, std::abs(candidate - dx[j]));
                dx[j] = candidate;
            }
            if (maxChange < 1e-13) break;
        }

        std::vector<Component> trial = comp;
        double maxStep = 0.0;
        for (int j = 0; j < n; ++j) {
            trial[j].sigma = comp[j].sigma + dx[j];
            if (trial[j].sigma < floorSigma) trial[j].sigma = floorSigma;
            if (trial[j].sigma > sigmaCap) trial[j].sigma = sigmaCap;
            maxStep = std::max(maxStep, std::abs(dx[j]));
        }

        const double newObj = objective(trial);
        if (newObj < obj) {
            comp = trial;
            obj = newObj;
            lambda = std::max(1e-12, lambda * 0.25);
            if (maxStep < 1e-12) break;
        } else {
            lambda = std::min(1e12, lambda * 10.0);
        }
    }

    CalibrationResult out;
    out.components = comp;
    out.objective = obj;
    out.converged = obj <= tol;
    out.iterations = it;
    out.inputArbitrageFree = marketArb;
    out.artificialSide = artificialSide;
    out.artificialK = artificialK;
    out.artificialCall = artificialCall;
    out.theoreticalMinArtificialWeight = theoreticalMinArtificialWeight;
    out.artificialCalibrationWeight =
        artificialSide == "none" ? 0.0 : ARTIFICIAL_QUOTE_CALIBRATION_WEIGHT;
    out.sigmaFloor = floorSigma;
    out.sigmaCap = sigmaCap;
    out.atmTotalVol = atmTotalVol;
    out.p1 = p1;
    out.cN = cN;
    out.boundarySlope = boundarySlope;

    if (artificialSide == "left") {
        out.artificialWeight = slope[1] - slope[0];
        out.adjacentMarketWeight = slope[2] - slope[1];
    } else if (artificialSide == "right") {
        out.adjacentMarketWeight = slope[nKnots - 1] - slope[nKnots - 2];
        out.artificialWeight = slope[nKnots] - slope[nKnots - 1];
    } else {
        out.artificialWeight = 0.0;
        out.adjacentMarketWeight = 0.0;
    }
    out.boundaryPairWeightSum = out.artificialWeight + out.adjacentMarketWeight;
    out.boundaryPairMaxWeight = std::max(out.artificialWeight, out.adjacentMarketWeight);
    out.boundaryPairSquaredObjective =
        out.artificialWeight * out.artificialWeight +
        out.adjacentMarketWeight * out.adjacentMarketWeight;

    for (int i = 0; i < nMarket; ++i) {
        out.market.push_back({km[i], vm[i], cm[i]});
        out.maxPriceError = std::max(
            out.maxPriceError,
            std::abs(mixtureCall(comp, km[i]) - cm[i]));
    }

    // The wing functions are attached only outside the first/last observed
    // market strikes, so this step does not alter the calibrated market fit.
    configureWingDecay(out);
    return out;
}

CalibrationResult LognormalMixtureCalibrator::calibrate(
    const std::vector<std::pair<double, double>>& kv,
    double eta,
    int maxIter,
    double tol) const {

    // The synthetic quote is needed for the mixture moment constraint, but
    // the market fit should be driven almost entirely by real observations.
    const double ARTIFICIAL_CALIBRATION_WEIGHT = 1e-10;
    return calibrateBalancedBoundary(
        kv, eta, maxIter, tol, ARTIFICIAL_CALIBRATION_WEIGHT);
}

} // namespace lnmix
