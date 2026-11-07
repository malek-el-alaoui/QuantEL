#pragma once
#include <string>
#include <utility>
#include <vector>

namespace lnmix {

struct Quote {
    double k = 0.0;          // strike / forward
    double totalVol = 0.0;   // annualized vol * sqrt(T)
    double call = 0.0;       // normalized undiscounted call
};

struct Component {
    double center = 0.0;
    double weight = 0.0;
    double sigma = 0.0;      // lognormal kernel total vol
};

struct WingArbitrageMetrics {
    bool directCallArbitrageFree = false;
    double minButterflyG = 0.0;
    double minDensity = 0.0;
    int negativeDensityPoints = 0;
    int invalidVariancePoints = 0;
};

struct CalibrationResult {
    std::vector<Quote> market;
    std::vector<Component> components;
    double objective = 0.0;
    double maxPriceError = 0.0;
    bool converged = false;
    bool inputArbitrageFree = false;
    int iterations = 0;

    // Diagnostics for the one-artificial-point construction.
    std::string artificialSide = "none";
    double artificialK = 0.0;
    double artificialCall = 0.0;
    double artificialWeight = 0.0;
    double theoreticalMinArtificialWeight = 0.0;
    double adjacentMarketWeight = 0.0;
    double boundaryPairWeightSum = 0.0;
    double boundaryPairMaxWeight = 0.0;
    double boundaryPairSquaredObjective = 0.0;
    double boundarySlope = 0.0;
    double artificialCalibrationWeight = 0.0;
    double sigmaFloor = 0.0;
    double sigmaCap = 0.0;
    double atmTotalVol = 0.0;
    double p1 = 0.0;
    double cN = 0.0;

    // Equity right-wing parametrization.
    // Every calibrated component sigma remains exactly flat on the observed
    // market interval [marketKMin, marketKMax].  The left wing is left flat.
    // On the right, each component receives a quadratic exponential decay
    // whose strength scales with (sigma_j / sigma_ATM)^2.  This targets broad
    // kernels that would otherwise dominate the remote equity right tail.
    double marketKMin = 0.0;
    double marketKMax = 0.0;
    double leftWingBaseSigma = 0.0;
    double rightWingBaseSigma = 0.0;
    double leftWingAlpha = 0.0;
    double rightWingAlpha = 0.0;
    double leftWingCriticalAlpha = 0.0;
    double rightWingCriticalAlpha = 0.0;
    double wingPower = 3.0;
    double wingSafetyXMin = 0.0;
    double wingSafetyXMax = 0.0;
    double wingSafetyScale = 1.0;
    double wingTargetSigmaRatio = 0.70;
    double rightWingMaxOvershootBp = 0.0;
    double rightWingMaxStepRiseBp = 0.0;
    double rightWingEndIv = 0.0;
    double rightWingBoundaryIv = 0.0;
    bool rightWingNonIncreasing = false;
    bool rightWingMonotoneTargetFeasible = false;
    double wingTargetDelta = 0.005;
    double wingMinButterflyG = 0.0;
    double wingMinDensity = 0.0;
    int wingNegativeDensityPoints = 0;
    int wingInvalidVariancePoints = 0;
    bool wingDirectCallArbitrageFree = false;
};

class LognormalMixtureCalibrator {
public:
    CalibrationResult calibrate(const std::vector<std::pair<double, double>>& kTotalVol,
                                double eta = 0.10,
                                int maxIter = 500,
                                double tol = 1e-14) const;

    static double normCdf(double x);
    static double normPdf(double x);
    static double invNormCdf(double p);

    static double bsCall(double k, double totalVol);
    static double bsVega(double k, double totalVol);
    static double shiftedCall(double center, double k, double sigma);
    static double shiftedVega(double center, double k, double sigma);
    static double impliedTotalVol(double k, double call);
    static double density(const std::vector<Component>& components, double k);
    static double mixtureCall(const std::vector<Component>& components, double k);

    // Wing-controlled Model 1 evaluators.
    static double wingSigma(const CalibrationResult& result,
                            std::size_t componentIndex,
                            double k);
    static double wingSigmaMultiplier(const CalibrationResult& result, double k);
    static double wingedMixtureCall(const CalibrationResult& result, double k);
    static double wingedImpliedTotalVariance(const CalibrationResult& result, double k);
    static double wingedButterflyFactor(const CalibrationResult& result, double k);
    static double wingedDensity(const CalibrationResult& result, double k);
    static WingArbitrageMetrics wingArbitrageMetrics(const CalibrationResult& result,
                                                     double xMin,
                                                     double xMax,
                                                     int points = 1601);

    static bool discreteArbitrageFree(const std::vector<double>& k,
                                      const std::vector<double>& c,
                                      double eps = 1e-10);

private:
    CalibrationResult calibrateBalancedBoundary(
        const std::vector<std::pair<double, double>>& kTotalVol,
        double eta,
        int maxIter,
        double tol,
        double artificialQuoteCalibrationWeight) const;

    static void configureWingDecay(CalibrationResult& result);
};

} // namespace lnmix
