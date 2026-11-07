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

struct CalibrationResult {
    std::vector<Quote> market;
    std::vector<Component> components;
    double objective = 0.0;
    double maxPriceError = 0.0;
    bool converged = false;
    bool inputArbitrageFree = false;
    int iterations = 0;

    // Diagnostics for the one-artificial-point construction.
    std::string artificialSide = "none"; // "left", "right", or "none"
    double artificialK = 0.0;
    double artificialCall = 0.0;
    double artificialWeight = 0.0;
    double theoreticalMinArtificialWeight = 0.0; // retained for compatibility/diagnostics
    double adjacentMarketWeight = 0.0; // butterfly mass at the neighboring market knot
    double boundaryPairWeightSum = 0.0; // fixed local mass split across boundary + neighbor
    double boundaryPairMaxWeight = 0.0; // minimized by the balanced construction
    double boundaryPairSquaredObjective = 0.0; // w_boundary^2 + w_neighbor^2
    double boundarySlope = 0.0; // slope between artificial quote and nearest market quote
    double artificialCalibrationWeight = 0.0; // objective weight of the artificial quote; market quotes use 1
    double sigmaFloor = 0.0;
    double p1 = 0.0;
    double cN = 0.0;
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
};

} // namespace lnmix
