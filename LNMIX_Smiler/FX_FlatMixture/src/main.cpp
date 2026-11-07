#include "cubic_variance.hpp"
#include "lognormal_mixture.hpp"
#include "market_data.hpp"
#include "multi_maturity.hpp"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <vector>
#ifdef _WIN32
#include <direct.h>
#endif

using namespace lnmix;

namespace {

// Dense smile/density outputs use a market-proportional wing range rather
// than a fixed log-moneyness box.  The left edge is extended to an
// approximate |put forward delta| of 2%, and the right edge to an
// approximate call forward delta of 2%, using the nearest observed market
// wing volatility.  Change this to 0.05 for 5-delta output.
const double OUTPUT_DELTA_CUTOFF = 0.02;
const int DENSE_GRID_POINTS = 401;

bool ensureOutputDirectory() {
#ifdef _WIN32
    const int rc = _mkdir("output");
#else
    const int rc = mkdir("output", 0755);
#endif
    return rc == 0 || errno == EEXIST;
}

bool fileExists(const std::string& path) {
    std::ifstream f(path.c_str(), std::ios::binary);
    return static_cast<bool>(f);
}

void appendCsv(const std::string& inputPath, std::ofstream& out, bool& wroteHeader) {
    if (!fileExists(inputPath)) return;
    std::ifstream in(inputPath.c_str());
    std::string line;
    bool first = true;
    while (std::getline(in, line)) {
        if (first) {
            if (!wroteHeader) {
                out << line << '\n';
                wroteHeader = true;
            }
            first = false;
        } else if (!line.empty()) {
            out << line << '\n';
        }
    }
}

double model1AnnualIv(const std::vector<Component>& c, double k, double t) {
    const double call = LognormalMixtureCalibrator::mixtureCall(c, k);
    return LognormalMixtureCalibrator::impliedTotalVol(k, call) / std::sqrt(t);
}

double forwardAt(const MarketSurface& m, double t) {
    if (m.slices.empty()) return 1.0;
    if (t <= m.slices.front().t) return m.slices.front().forward;
    if (t >= m.slices.back().t) return m.slices.back().forward;
    for (std::size_t i = 0; i + 1 < m.slices.size(); ++i) {
        if (t >= m.slices[i].t && t <= m.slices[i + 1].t) {
            const double q = (t - m.slices[i].t) / (m.slices[i + 1].t - m.slices[i].t);
            return (1.0 - q) * m.slices[i].forward + q * m.slices[i + 1].forward;
        }
    }
    return m.slices.back().forward;
}

std::string safeFinite(double x) {
    if (!std::isfinite(x)) return "";
    std::ostringstream ss;
    ss << std::setprecision(15) << x;
    return ss.str();
}

void dumpMarketVsModels(const std::string& prefix,
                        const MultiMaturitySurface& model1,
                        const CubicVarianceSurface& model2) {
    std::ofstream f(("output/" + prefix + "_market_vs_models.csv").c_str());
    f << std::setprecision(15);
    f << "instrument,asof,maturity,ttm,forward,strike,k,log_moneyness,pillar,"
         "market_iv,model1_lnmix_iv,model2_cubic_iv,model1_iv_error_bp,model2_iv_error_bp,"
         "market_call_normalized,model1_call_normalized,model2_call_normalized,"
         "model1_price_error_normalized,model2_price_error_normalized,"
         "model1_density,model2_density,model2_butterfly_g,model2_butterfly_ok,"
         "source_call_mid,source_iv\n";

    const std::vector<CalibratedSlice>& s1 = model1.slices();
    const std::vector<CubicVarianceSlice>& s2 = model2.slices();
    if (s1.size() != s2.size()) throw std::runtime_error("Internal slice mismatch between models");

    for (std::size_t j = 0; j < s1.size(); ++j) {
        const CalibratedSlice& a = s1[j];
        const CubicVarianceSlice& b = s2[j];
        for (std::size_t i = 0; i < a.market.points.size(); ++i) {
            const MarketPoint& p = a.market.points[i];
            const double marketCall = LognormalMixtureCalibrator::bsCall(p.k, p.marketVol * std::sqrt(a.market.t));
            const double call1 = LognormalMixtureCalibrator::mixtureCall(a.calibration.components, p.k);
            const double iv1 = LognormalMixtureCalibrator::impliedTotalVol(p.k, call1) / std::sqrt(a.market.t);
            const double call2 = b.call(p.k);
            const double iv2 = b.annualVol(p.k);
            const double d1 = LognormalMixtureCalibrator::density(a.calibration.components, p.k);
            const double d2 = b.density(p.k);
            const double g2 = b.butterflyFactor(p.k);
            const bool ok2 = std::isfinite(d2) && std::isfinite(g2) && d2 >= -1e-10 && g2 >= -1e-10;

            f << a.market.instrument << ',' << a.market.asOf << ',' << a.market.maturity << ','
              << a.market.t << ',' << a.market.forward << ',' << p.strike << ',' << p.k << ',' << std::log(p.k) << ',' << p.pillar << ','
              << p.marketVol << ',' << iv1 << ',' << iv2 << ','
              << (iv1 - p.marketVol) * 10000.0 << ',' << (iv2 - p.marketVol) * 10000.0 << ','
              << marketCall << ',' << call1 << ',' << call2 << ','
              << (call1 - marketCall) << ',' << (call2 - marketCall) << ','
              << d1 << ',' << safeFinite(d2) << ',' << safeFinite(g2) << ',' << (ok2 ? 1 : 0) << ',';
            if (std::isfinite(p.sourcePrice)) f << p.sourcePrice;
            f << ',';
            if (std::isfinite(p.sourceIv)) f << p.sourceIv;
            f << '\n';
        }
    }
}

void dumpModel1Components(const std::string& prefix, const MultiMaturitySurface& model1) {
    std::ofstream f(("output/" + prefix + "_model1_components.csv").c_str());
    f << std::setprecision(15);
    f << "instrument,maturity,ttm,component,center_k,weight,kernel_total_vol,kernel_annualized_vol,is_artificial_component\n";
    for (std::size_t j = 0; j < model1.slices().size(); ++j) {
        const CalibratedSlice& s = model1.slices()[j];
        for (std::size_t i = 0; i < s.calibration.components.size(); ++i) {
            const Component& c = s.calibration.components[i];
            const bool artificial = s.calibration.artificialSide != "none" &&
                std::abs(c.center - s.calibration.artificialK) <=
                1e-10 * std::max(1.0, std::abs(s.calibration.artificialK));
            f << s.market.instrument << ',' << s.market.maturity << ',' << s.market.t << ',' << i << ','
              << c.center << ',' << c.weight << ',' << c.sigma << ',' << c.sigma / std::sqrt(s.market.t) << ','
              << (artificial ? 1 : 0) << '\n';
        }
    }
}

void dumpArtificialQuoteDiagnostics(const std::string& prefix,
                                    const MultiMaturitySurface& model1) {
    std::ofstream f(("output/" + prefix + "_model1_artificial_quote_diagnostics.csv").c_str());
    f << std::setprecision(15);
    f << "instrument,maturity,ttm,forward,p1,cN,artificial_side,artificial_k,artificial_strike_equivalent,"
         "artificial_call,artificial_weight,adjacent_market_weight,boundary_pair_weight_sum,"
         "boundary_pair_max_weight,boundary_pair_squared_objective,boundary_slope,"
         "theoretical_min_artificial_weight,weight_over_theoretical_min,"
         "artificial_calibration_objective_weight,sigma_floor\n";
    for (std::size_t j = 0; j < model1.slices().size(); ++j) {
        const CalibratedSlice& s = model1.slices()[j];
        const CalibrationResult& r = s.calibration;
        f << s.market.instrument << ',' << s.market.maturity << ',' << s.market.t << ',' << s.market.forward << ','
          << r.p1 << ',' << r.cN << ',' << r.artificialSide << ',' << r.artificialK << ','
          << r.artificialK * s.market.forward << ',' << r.artificialCall << ',' << r.artificialWeight << ','
          << r.adjacentMarketWeight << ',' << r.boundaryPairWeightSum << ',' << r.boundaryPairMaxWeight << ','
          << r.boundaryPairSquaredObjective << ',' << r.boundarySlope << ','
          << r.theoreticalMinArtificialWeight << ',';
        if (r.theoreticalMinArtificialWeight > 0.0)
            f << r.artificialWeight / r.theoreticalMinArtificialWeight;
        f << ',' << r.artificialCalibrationWeight << ',' << r.sigmaFloor << '\n';
    }
}

void sliceGridBounds(const MarketSlice& s, double& kMin, double& kMax) {
    if (s.points.empty()) throw std::runtime_error("Empty market slice");
    kMin = s.points.front().k;
    kMax = s.points.front().k;
    for (std::size_t i = 1; i < s.points.size(); ++i) {
        kMin = std::min(kMin, s.points[i].k);
        kMax = std::max(kMax, s.points[i].k);
    }
}

struct DeltaRangeInfo {
    double marketXMin = 0.0;
    double marketXMax = 0.0;
    double xMin = 0.0;
    double xMax = 0.0;
    double marketKMin = 0.0;
    double marketKMax = 0.0;
    double rangeKMin = 0.0;
    double rangeKMax = 0.0;
    double leftPutDeltaProxy = 0.0;
    double rightCallDeltaProxy = 0.0;
};

DeltaRangeInfo deltaRangeBounds(const MarketSlice& s,
                                double deltaCutoff = OUTPUT_DELTA_CUTOFF) {
    if (s.points.empty()) throw std::runtime_error("Empty market slice");
    if (!(deltaCutoff > 0.0 && deltaCutoff < 0.5))
        throw std::runtime_error("Delta cutoff must lie in (0,0.5)");

    std::size_t left = 0, right = 0;
    for (std::size_t i = 1; i < s.points.size(); ++i) {
        if (s.points[i].k < s.points[left].k) left = i;
        if (s.points[i].k > s.points[right].k) right = i;
    }

    DeltaRangeInfo out;
    out.marketKMin = s.points[left].k;
    out.marketKMax = s.points[right].k;
    out.marketXMin = std::log(out.marketKMin);
    out.marketXMax = std::log(out.marketKMax);

    const double rootT = std::sqrt(s.t);
    const double vLeft = std::max(1e-10, s.points[left].marketVol * rootT);
    const double vRight = std::max(1e-10, s.points[right].marketVol * rootT);

    // Forward-delta convention:
    // call delta = N(d1), |put delta| = 1-N(d1),
    // d1 = -log(k)/v + v/2.
    const double d1Left = LognormalMixtureCalibrator::invNormCdf(1.0 - deltaCutoff);
    const double d1Right = LognormalMixtureCalibrator::invNormCdf(deltaCutoff);
    const double kAtLeftDelta = std::exp(0.5 * vLeft * vLeft - vLeft * d1Left);
    const double kAtRightDelta = std::exp(0.5 * vRight * vRight - vRight * d1Right);

    // Never cut away an observed market quote.  If the market already extends
    // beyond 2-delta under this proxy, retain the full observed range.
    out.rangeKMin = std::min(out.marketKMin, kAtLeftDelta);
    out.rangeKMax = std::max(out.marketKMax, kAtRightDelta);
    out.xMin = std::log(out.rangeKMin);
    out.xMax = std::log(out.rangeKMax);

    const double edgeD1Left = -std::log(out.rangeKMin) / vLeft + 0.5 * vLeft;
    const double edgeD1Right = -std::log(out.rangeKMax) / vRight + 0.5 * vRight;
    out.leftPutDeltaProxy = 1.0 - LognormalMixtureCalibrator::normCdf(edgeD1Left);
    out.rightCallDeltaProxy = LognormalMixtureCalibrator::normCdf(edgeD1Right);
    return out;
}

DeltaRangeInfo deltaRangeBoundsAtTime(const MultiMaturitySurface& model1, double t) {
    const std::vector<CalibratedSlice>& slices = model1.slices();
    if (slices.empty()) throw std::runtime_error("Empty model surface");
    if (t <= slices.front().market.t) return deltaRangeBounds(slices.front().market);
    if (t >= slices.back().market.t) return deltaRangeBounds(slices.back().market);

    for (std::size_t j = 0; j + 1 < slices.size(); ++j) {
        if (t >= slices[j].market.t && t <= slices[j + 1].market.t) {
            const DeltaRangeInfo a = deltaRangeBounds(slices[j].market);
            const DeltaRangeInfo b = deltaRangeBounds(slices[j + 1].market);
            const double q = (t - slices[j].market.t) /
                             (slices[j + 1].market.t - slices[j].market.t);
            DeltaRangeInfo r;
            r.xMin = (1.0 - q) * a.xMin + q * b.xMin;
            r.xMax = (1.0 - q) * a.xMax + q * b.xMax;
            r.rangeKMin = std::exp(r.xMin);
            r.rangeKMax = std::exp(r.xMax);
            r.marketXMin = (1.0 - q) * a.marketXMin + q * b.marketXMin;
            r.marketXMax = (1.0 - q) * a.marketXMax + q * b.marketXMax;
            r.marketKMin = std::exp(r.marketXMin);
            r.marketKMax = std::exp(r.marketXMax);
            r.leftPutDeltaProxy = OUTPUT_DELTA_CUTOFF;
            r.rightCallDeltaProxy = OUTPUT_DELTA_CUTOFF;
            return r;
        }
    }
    return deltaRangeBounds(slices.back().market);
}

void dumpContinuousSmiles(const std::string& prefix,
                          const MultiMaturitySurface& model1,
                          const CubicVarianceSurface& model2,
                          int gridPoints = DENSE_GRID_POINTS) {
    std::ofstream f(("output/" + prefix + "_continuous_smiles.csv").c_str());
    f << std::setprecision(15);
    f << "instrument,asof,maturity,ttm,forward,k,log_moneyness,strike_equivalent,"
         "market_k_min,market_k_max,range_k_min,range_k_max,range_delta_cutoff,"
         "left_put_delta_proxy,right_call_delta_proxy,is_extrapolation,"
         "model1_lnmix_iv,model2_cubic_iv,iv_difference_bp,"
         "model1_call_normalized,model2_call_normalized,call_difference_normalized,"
         "model2_raw_total_variance,model2_variance_valid\n";

    const std::vector<CalibratedSlice>& s1 = model1.slices();
    const std::vector<CubicVarianceSlice>& s2 = model2.slices();
    for (std::size_t j = 0; j < s1.size(); ++j) {
        const DeltaRangeInfo range = deltaRangeBounds(s1[j].market);

        for (int i = 0; i < gridPoints; ++i) {
            const double x = range.xMin + (range.xMax - range.xMin) * i / static_cast<double>(gridPoints - 1);
            const double k = std::exp(x);
            const bool extrap = x < range.marketXMin - 1e-14 || x > range.marketXMax + 1e-14;

            const double call1 = LognormalMixtureCalibrator::mixtureCall(s1[j].calibration.components, k);
            const double iv1 = LognormalMixtureCalibrator::impliedTotalVol(k, call1) / std::sqrt(s1[j].market.t);

            const double rawW2 = s2[j].rawTotalVariance(k);
            const bool validW2 = std::isfinite(rawW2) && rawW2 > 0.0;
            const double iv2 = validW2 ? std::sqrt(rawW2 / s1[j].market.t) : std::numeric_limits<double>::quiet_NaN();
            const double call2 = validW2 ? LognormalMixtureCalibrator::bsCall(k, std::sqrt(rawW2))
                                         : std::numeric_limits<double>::quiet_NaN();

            f << s1[j].market.instrument << ',' << s1[j].market.asOf << ',' << s1[j].market.maturity << ','
              << s1[j].market.t << ',' << s1[j].market.forward << ',' << k << ',' << x << ','
              << k * s1[j].market.forward << ',' << range.marketKMin << ',' << range.marketKMax << ','
              << range.rangeKMin << ',' << range.rangeKMax << ',' << OUTPUT_DELTA_CUTOFF << ','
              << range.leftPutDeltaProxy << ',' << range.rightCallDeltaProxy << ',' << (extrap ? 1 : 0) << ','
              << safeFinite(iv1) << ',' << safeFinite(iv2) << ',';
            if (std::isfinite(iv1) && std::isfinite(iv2)) f << (iv1 - iv2) * 10000.0;
            f << ',' << safeFinite(call1) << ',' << safeFinite(call2) << ',';
            if (std::isfinite(call1) && std::isfinite(call2)) f << (call1 - call2);
            f << ',' << safeFinite(rawW2) << ',' << (validW2 ? 1 : 0) << '\n';
        }
    }
}

void dumpDensities(const std::string& prefix,
                   const MultiMaturitySurface& model1,
                   const CubicVarianceSurface& model2,
                   int gridPoints = DENSE_GRID_POINTS) {
    std::ofstream f(("output/" + prefix + "_densities.csv").c_str());
    f << std::setprecision(15);
    f << "instrument,asof,maturity,ttm,forward,k,log_moneyness,strike_equivalent,"
         "market_k_min,market_k_max,range_k_min,range_k_max,range_delta_cutoff,"
         "left_put_delta_proxy,right_call_delta_proxy,is_extrapolation,"
         "model1_density,model2_density,density_difference,model2_butterfly_g,"
         "model2_butterfly_ok,model2_raw_total_variance,model2_variance_valid\n";

    const std::vector<CalibratedSlice>& s1 = model1.slices();
    const std::vector<CubicVarianceSlice>& s2 = model2.slices();
    for (std::size_t j = 0; j < s1.size(); ++j) {
        const DeltaRangeInfo range = deltaRangeBounds(s1[j].market);
        for (int i = 0; i < gridPoints; ++i) {
            const double x = range.xMin + (range.xMax - range.xMin) * i / static_cast<double>(gridPoints - 1);
            const double k = std::exp(x);
            const bool extrap = x < range.marketXMin - 1e-14 || x > range.marketXMax + 1e-14;
            const double d1 = LognormalMixtureCalibrator::density(s1[j].calibration.components, k);
            const double d2 = s2[j].density(k);
            const double g2 = s2[j].butterflyFactor(k);
            const double rawW2 = s2[j].rawTotalVariance(k);
            const bool validW2 = std::isfinite(rawW2) && rawW2 > 0.0;
            const bool ok2 = validW2 && std::isfinite(d2) && std::isfinite(g2) && d2 >= -1e-10 && g2 >= -1e-10;
            f << s1[j].market.instrument << ',' << s1[j].market.asOf << ',' << s1[j].market.maturity << ','
              << s1[j].market.t << ',' << s1[j].market.forward << ',' << k << ',' << x << ',' << k * s1[j].market.forward << ','
              << range.marketKMin << ',' << range.marketKMax << ',' << range.rangeKMin << ',' << range.rangeKMax << ','
              << OUTPUT_DELTA_CUTOFF << ',' << range.leftPutDeltaProxy << ',' << range.rightCallDeltaProxy << ','
              << (extrap ? 1 : 0) << ',' << d1 << ',' << safeFinite(d2) << ',';
            if (std::isfinite(d2)) f << (d1 - d2);
            f << ',' << safeFinite(g2) << ',' << (ok2 ? 1 : 0) << ',' << safeFinite(rawW2) << ',' << (validW2 ? 1 : 0) << '\n';
        }
    }
}

void dumpMultiMaturityComparison(const std::string& prefix,
                                 const MultiMaturitySurface& model1,
                                 const CubicVarianceSurface& model2,
                                 int gridPoints = 161) {
    if (model1.slices().size() < 2) return;
    std::ofstream f(("output/" + prefix + "_multi_maturity_comparison.csv").c_str());
    f << std::setprecision(15);
    f << "instrument,t,forward,k,strike_equivalent,range_k_min,range_k_max,range_delta_cutoff,"
         "model1_iv,model2_iv,iv_difference_bp,"
         "model1_density,model2_density,model1_local_vol,model2_local_vol,"
         "model1_call,model2_call,call_difference\n";

    std::vector<double> times;
    for (std::size_t j = 0; j < model1.slices().size(); ++j) {
        times.push_back(model1.slices()[j].market.t);
        if (j + 1 < model1.slices().size())
            times.push_back(0.5 * (model1.slices()[j].market.t + model1.slices()[j + 1].market.t));
    }
    std::sort(times.begin(), times.end());

    for (std::size_t ti = 0; ti < times.size(); ++ti) {
        const double t = times[ti];
        const double F = forwardAt(model1.market(), t);
        const DeltaRangeInfo range = deltaRangeBoundsAtTime(model1, t);
        for (int i = 0; i < gridPoints; ++i) {
            const double x = range.xMin + (range.xMax - range.xMin) * i / static_cast<double>(gridPoints - 1);
            const double k = std::exp(x);
            const double iv1 = model1.impliedVol(k, t);
            const double iv2 = model2.impliedVol(k, t);
            const double d1 = model1.density(k, t);
            const double d2 = model2.density(k, t);
            const double lv1 = model1.localVol(k, t);
            const double lv2 = model2.localVol(k, t);
            const double c1 = model1.call(k, t);
            const double c2 = model2.call(k, t);
            f << model1.market().instrument << ',' << t << ',' << F << ',' << k << ',' << k * F << ','
              << range.rangeKMin << ',' << range.rangeKMax << ',' << OUTPUT_DELTA_CUTOFF << ','
              << safeFinite(iv1) << ',' << safeFinite(iv2) << ',';
            if (std::isfinite(iv1) && std::isfinite(iv2)) f << (iv1 - iv2) * 10000.0;
            f << ',' << safeFinite(d1) << ',' << safeFinite(d2) << ',' << safeFinite(lv1) << ',' << safeFinite(lv2) << ','
              << safeFinite(c1) << ',' << safeFinite(c2) << ',';
            if (std::isfinite(c1) && std::isfinite(c2)) f << (c1 - c2);
            f << '\n';
        }
    }
}

void dumpCalendarChecks(const std::string& prefix,
                        const MultiMaturitySurface& model1,
                        const CubicVarianceSurface& model2) {
    if (model1.slices().size() < 2) return;
    std::ofstream f(("output/" + prefix + "_calendar_checks.csv").c_str());
    f << std::setprecision(15);
    f << "instrument,left_maturity,right_maturity,k,model1_left_call,model1_right_call,model1_difference,model1_ok,"
         "model2_left_call,model2_right_call,model2_difference,model2_ok\n";

    const std::vector<CalibratedSlice>& a = model1.slices();
    const std::vector<CubicVarianceSlice>& b = model2.slices();
    double kMin = 1e100, kMax = 0.0;
    for (std::size_t j = 0; j < a.size(); ++j) {
        double lo, hi;
        sliceGridBounds(a[j].market, lo, hi);
        kMin = std::min(kMin, lo);
        kMax = std::max(kMax, hi);
    }
    kMin = std::max(kMin, 0.6);
    kMax = std::min(kMax, 1.5);

    for (std::size_t j = 0; j + 1 < a.size(); ++j) {
        for (int i = 0; i <= 180; ++i) {
            const double k = kMin + (kMax - kMin) * i / 180.0;
            const double a1 = LognormalMixtureCalibrator::mixtureCall(a[j].calibration.components, k);
            const double a2 = LognormalMixtureCalibrator::mixtureCall(a[j + 1].calibration.components, k);
            const double b1 = b[j].call(k);
            const double b2 = b[j + 1].call(k);
            f << model1.market().instrument << ',' << a[j].market.maturity << ',' << a[j + 1].market.maturity << ',' << k << ','
              << a1 << ',' << a2 << ',' << (a2 - a1) << ',' << ((a2 + 1e-10 >= a1) ? 1 : 0) << ','
              << b1 << ',' << b2 << ',' << (b2 - b1) << ',' << ((b2 + 1e-10 >= b1) ? 1 : 0) << '\n';
        }
    }
}

struct SummaryRow {
    std::string dataset;
    std::string maturity;
    double t;
    int points;
    bool inputArb;
    bool model1Converged;
    int model1Iterations;
    double model1MaxIvBp;
    double model2MaxIvBp;
    double model2MinDensity;
    double model2MinG;
    int model2NegativeDensityPoints;
    int model2InvalidVariancePoints;
    std::string artificialSide;
    double artificialWeight;
    double adjacentMarketWeight;
    double boundaryPairMaxWeight;
    double boundaryPairSquaredObjective;
    double theoreticalMinArtificialWeight;
    double artificialCalibrationWeight;
    double artificialK;
    double rangeKMin;
    double rangeKMax;
    bool model1Calendar;
    bool model2Calendar;
};

std::vector<SummaryRow> runDataset(const std::string& prefix, const MarketSurface& market, double eta) {
    MultiMaturitySurface model1(market);
    CubicVarianceSurface model2(market);
    model1.calibrate(eta);
    model2.calibrate();

    dumpMarketVsModels(prefix, model1, model2);
    dumpModel1Components(prefix, model1);
    dumpArtificialQuoteDiagnostics(prefix, model1);
    dumpContinuousSmiles(prefix, model1, model2);
    dumpDensities(prefix, model1, model2);
    dumpMultiMaturityComparison(prefix, model1, model2);
    dumpCalendarChecks(prefix, model1, model2);

    const bool cal1 = model1.calendarArbitrageFree(0.75, 1.25, 201);
    const bool cal2 = model2.calendarArbitrageFree(0.75, 1.25, 201);
    std::vector<SummaryRow> rows;

    for (std::size_t j = 0; j < model1.slices().size(); ++j) {
        const CalibratedSlice& a = model1.slices()[j];
        const CubicVarianceSlice& b = model2.slices()[j];
        double max1 = 0.0, max2 = 0.0;
        for (std::size_t i = 0; i < a.market.points.size(); ++i) {
            const MarketPoint& p = a.market.points[i];
            const double iv1 = model1AnnualIv(a.calibration.components, p.k, a.market.t);
            const double iv2 = b.annualVol(p.k);
            max1 = std::max(max1, std::abs(iv1 - p.marketVol) * 10000.0);
            max2 = std::max(max2, std::abs(iv2 - p.marketVol) * 10000.0);
        }

        const DeltaRangeInfo range = deltaRangeBounds(a.market);
        double minD2 = std::numeric_limits<double>::infinity();
        double minG2 = std::numeric_limits<double>::infinity();
        int neg = 0;
        int invalidW = 0;
        for (int i = 0; i < DENSE_GRID_POINTS; ++i) {
            const double k = std::exp(range.xMin + (range.xMax - range.xMin) * i /
                                      static_cast<double>(DENSE_GRID_POINTS - 1));
            const double rawW = b.rawTotalVariance(k);
            if (!std::isfinite(rawW) || rawW <= 0.0) ++invalidW;
            const double d2 = b.density(k);
            const double g2 = b.butterflyFactor(k);
            if (std::isfinite(d2)) {
                minD2 = std::min(minD2, d2);
                if (d2 < -1e-10) ++neg;
            }
            if (std::isfinite(g2)) minG2 = std::min(minG2, g2);
        }
        if (!std::isfinite(minD2)) minD2 = std::numeric_limits<double>::quiet_NaN();
        if (!std::isfinite(minG2)) minG2 = std::numeric_limits<double>::quiet_NaN();

        rows.push_back(SummaryRow{prefix, a.market.maturity, a.market.t,
                                  static_cast<int>(a.market.points.size()), a.calibration.inputArbitrageFree,
                                  a.calibration.converged, a.calibration.iterations, max1, max2,
                                  minD2, minG2, neg, invalidW,
                                  a.calibration.artificialSide, a.calibration.artificialWeight,
                                  a.calibration.adjacentMarketWeight, a.calibration.boundaryPairMaxWeight,
                                  a.calibration.boundaryPairSquaredObjective,
                                  a.calibration.theoreticalMinArtificialWeight, a.calibration.artificialCalibrationWeight,
                                  a.calibration.artificialK, range.rangeKMin, range.rangeKMax, cal1, cal2});

        std::cout << std::left << std::setw(13) << prefix << "  " << std::setw(12) << a.market.maturity
                  << " points=" << std::setw(2) << a.market.points.size()
                  << "  M1 maxIV=" << std::fixed << std::setprecision(4) << max1 << " bp"
                  << "  M2 maxIV=" << max2 << " bp"
                  << "  M2 minDensity=" << std::scientific << minD2
                  << "  M2 butterfly=" << (neg == 0 ? "OK" : "NEGATIVE")
                  << "  boundary=" << a.calibration.artificialSide
                  << " mixW=" << a.calibration.artificialWeight
                  << " calibW=" << a.calibration.artificialCalibrationWeight
                  << "  rangeK=[" << range.rangeKMin << ',' << range.rangeKMax << "]\n";
    }

    if (model1.slices().size() > 1) {
        std::cout << "  calendar Model 1: " << (cal1 ? "PASS" : "FAIL")
                  << " | Model 2: " << (cal2 ? "PASS" : "FAIL") << "\n";
    }
    std::cout << '\n';
    return rows;
}

void writeCombinedFiles() {
    const std::vector<std::string> prefixes = {"eurostoxx50", "spx", "eurusd", "usdhkd", "aapl", "nvda", "tsla"};
    struct Pair { const char* suffix; const char* output; };
    const Pair jobs[] = {
        {"_market_vs_models.csv", "output/all_market_vs_models.csv"},
        {"_model1_components.csv", "output/all_model1_components.csv"},
        {"_model1_artificial_quote_diagnostics.csv", "output/all_model1_artificial_quote_diagnostics.csv"},
        {"_continuous_smiles.csv", "output/all_continuous_smiles.csv"},
        {"_densities.csv", "output/all_densities.csv"},
        {"_multi_maturity_comparison.csv", "output/all_multi_maturity_comparison.csv"}
    };

    for (std::size_t j = 0; j < sizeof(jobs) / sizeof(jobs[0]); ++j) {
        std::ofstream out(jobs[j].output);
        bool wroteHeader = false;
        for (std::size_t i = 0; i < prefixes.size(); ++i)
            appendCsv("output/" + prefixes[i] + jobs[j].suffix, out, wroteHeader);
    }
}

} // namespace

int main(int argc, char** argv) {
    try {
        if (!ensureOutputDirectory()) throw std::runtime_error("Cannot create output directory");
        const std::string mode = argc > 1 ? argv[1] : "all";
        const double eta = argc > 2 ? std::stod(argv[2]) : 0.10;
        if (!(eta > 0.0 && eta < 1.0)) throw std::runtime_error("eta must be in (0,1), e.g. 0.10 or 0.02");
        std::vector<SummaryRow> summary;

        std::cout << "Two-model volatility-smile comparison\n"
                     "Model 1: balanced-boundary arbitrage-free log-normal mixture\n"
                     "Model 2: natural cubic spline of total implied variance vs log-moneyness\n"
                     "==========================================================================\n\n";
        std::cout << "Kernel lower-bound eta = " << eta << "\n\n";

        if (mode == "all" || mode == "eurostoxx") {
            const MarketSurface m = loadStrikeSurfaceCsv("data/eurostoxx50_20140930.csv");
            const std::vector<SummaryRow> r = runDataset("eurostoxx50", m, eta);
            summary.insert(summary.end(), r.begin(), r.end());
        }
        if (mode == "all" || mode == "spx") {
            const MarketSurface m = loadStrikeSurfaceCsv("data/spx_20260730.csv");
            const std::vector<SummaryRow> r = runDataset("spx", m, eta);
            summary.insert(summary.end(), r.begin(), r.end());
        }
        if (mode == "all" || mode == "eurusd") {
            const MarketSurface m = loadFxDeltaSurfaceCsv("data/eurusd_20240515_delta.csv");
            const std::vector<SummaryRow> r = runDataset("eurusd", m, eta);
            summary.insert(summary.end(), r.begin(), r.end());
        }
        if (mode == "all" || mode == "usdhkd") {
            const MarketSurface m = loadFxDeltaSurfaceCsv("data/usdhkd_20250521_delta.csv");
            const std::vector<SummaryRow> r = runDataset("usdhkd", m, eta);
            summary.insert(summary.end(), r.begin(), r.end());
        }
        if (mode == "all" || mode == "aapl") {
            const MarketSurface m = loadStrikeSurfaceCsv("data/aapl_cboe_20260910.csv");
            const std::vector<SummaryRow> r = runDataset("aapl", m, eta);
            summary.insert(summary.end(), r.begin(), r.end());
        }
        if (mode == "all" || mode == "nvda") {
            const MarketSurface m = loadStrikeSurfaceCsv("data/nvda_cboe_20260910.csv");
            const std::vector<SummaryRow> r = runDataset("nvda", m, eta);
            summary.insert(summary.end(), r.begin(), r.end());
        }
        if (mode == "all" || mode == "tsla") {
            const MarketSurface m = loadStrikeSurfaceCsv("data/tsla_cboe_20260911.csv");
            const std::vector<SummaryRow> r = runDataset("tsla", m, eta);
            summary.insert(summary.end(), r.begin(), r.end());
        }
        if (summary.empty()) throw std::runtime_error("Unknown mode. Use all, eurostoxx, spx, eurusd, usdhkd, aapl, nvda, or tsla. Optional second argument: eta (0.10 or 0.02).");

        std::ofstream f("output/summary.csv");
        f << std::setprecision(15);
        f << "dataset,maturity,ttm,points,input_static_arb_free,model1_converged,model1_iterations,"
             "model1_max_iv_error_bp,model2_max_iv_error_bp,model2_min_density_delta_range,model2_min_butterfly_g_delta_range,"
             "model2_negative_density_delta_grid_points,model2_invalid_variance_delta_grid_points,"
             "model1_artificial_side,model1_artificial_weight,model1_adjacent_market_weight,"
             "model1_boundary_pair_max_weight,model1_boundary_pair_squared_objective,model1_theoretical_min_artificial_weight,"
             "model1_artificial_calibration_weight,model1_artificial_k,range_k_min,range_k_max,range_delta_cutoff,"
             "model1_calendar_ok,model2_calendar_ok\n";
        for (std::size_t i = 0; i < summary.size(); ++i) {
            const SummaryRow& r = summary[i];
            f << r.dataset << ',' << r.maturity << ',' << r.t << ',' << r.points << ',' << (r.inputArb ? 1 : 0) << ','
              << (r.model1Converged ? 1 : 0) << ',' << r.model1Iterations << ',' << r.model1MaxIvBp << ',' << r.model2MaxIvBp << ','
              << safeFinite(r.model2MinDensity) << ',' << safeFinite(r.model2MinG) << ',' << r.model2NegativeDensityPoints << ','
              << r.model2InvalidVariancePoints << ',' << r.artificialSide << ',' << r.artificialWeight << ','
              << r.adjacentMarketWeight << ',' << r.boundaryPairMaxWeight << ',' << r.boundaryPairSquaredObjective << ','
              << r.theoreticalMinArtificialWeight << ',' << r.artificialCalibrationWeight << ',' << r.artificialK << ','
              << r.rangeKMin << ',' << r.rangeKMax << ',' << OUTPUT_DELTA_CUTOFF << ','
              << (r.model1Calendar ? 1 : 0) << ',' << (r.model2Calendar ? 1 : 0) << '\n';
        }

        writeCombinedFiles();

        std::ofstream readme("output/README_RESULTS.txt");
        readme <<
            "MAIN FILES\n"
            "==========\n"
            "all_market_vs_models.csv       Market quotes vs Model 1 and Model 2.\n"
            "all_model1_artificial_quote_diagnostics.csv  Boundary quote / minimum-weight diagnostics.\n"
            "all_continuous_smiles.csv       Delta-range dense smile grids incl. extrapolation: Model 1 vs Model 2.\n"
            "all_densities.csv               Delta-range dense risk-neutral density grids incl. extrapolation.\n"
            "all_multi_maturity_comparison.csv Both models at quoted and intermediate maturities.\n"
            "summary.csv                     Calibration and arbitrage diagnostics.\n\n"
            "MODEL DEFINITIONS\n"
            "=================\n"
            "Model 1 = one-point paper log-normal mixture with a balanced boundary mass split.\n"
            "The artificial strike is chosen to minimize the larger of the boundary and adjacent-market butterfly masses; equivalently it minimizes their squared sum.\n"
            "Because their sum is fixed, the unconstrained optimum equalizes the two local butterfly masses.\n"
            "The artificial quote has only 1e-10 objective weight in calibration; every real market quote has weight 1.\n"
            "All Model 1 kernel total vols satisfy sigma_i > 0.10*sigma_ATM inside the constrained linear solve.\n"
            "Model 2 = natural cubic spline of total implied variance w(x)=IV^2*T vs x=log(K/F).\n"
            "The Model 2 density is computed analytically from w, w' and w'' using the standard butterfly factor g(x).\n"
            "Unlike Model 1, an unconstrained cubic spline is NOT guaranteed to be butterfly-arbitrage-free.\n"
            "Negative Model 2 densities are intentionally left visible in the CSV instead of being clipped.\n"
            "Dense outputs use a market-proportional 2-delta range based on the nearest observed wing volatilities.\n"
            "Change OUTPUT_DELTA_CUTOFF in src/main.cpp from 0.02 to 0.05 for a 5-delta range.\n"
            "Rows outside the observed quote interval are marked is_extrapolation=1.\n"
            "If cubic extrapolation makes raw total variance non-positive, Model 2 IV is left blank and model2_variance_valid=0.\n\n"
            "PLOTTING\n"
            "========\n"
            "Run python/plot_two_models.py after installing pandas and matplotlib.\n";

        std::cout << "CSV results written to ./output\n"
                     "  - all_market_vs_models.csv\n"
                     "  - all_model1_artificial_quote_diagnostics.csv\n"
                     "  - all_continuous_smiles.csv\n"
                     "  - all_densities.csv\n"
                     "  - summary.csv\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << '\n';
        return 1;
    }
}
