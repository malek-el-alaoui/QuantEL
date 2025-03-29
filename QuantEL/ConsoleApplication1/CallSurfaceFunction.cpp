#include "CallSurfaceFunction.h"
#include <cmath>
#include <stdexcept>

// CDF of standard normal
double norm_cdf(double x) {
	return 0.5 * std::erfc(-x / std::sqrt(2));
}

CallSurfaceFunction::CallSurfaceFunction(const VolSurfaceInterpolator& surface, double spot_, double rate_)
	: volSurface(surface), spot(spot_), rate(rate_) {}

double CallSurfaceFunction::Eval(double K, double t) const {
	if (K <= 0 || t <= 0)
		return 0.0;

	double sigma = volSurface.GetImpliedVol(K, t);
	if (sigma <= 0.0)
		return 0.0;

	double sqrtT = std::sqrt(t);
	double d1 = (std::log(spot / K) + (rate + 0.5 * sigma * sigma) * t) / (sigma * sqrtT);
	double d2 = d1 - sigma * sqrtT;

	double callPrice = spot * norm_cdf(d1) - K * std::exp(-rate * t) * norm_cdf(d2);
	return callPrice;
}
