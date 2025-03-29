#include "VolSurfaceInterpolator.h"
#include <stdexcept>
#include <cmath>
#include <algorithm> // for std::min and std::max
#include <cstddef>   // for size_t

VolSurfaceInterpolator::VolSurfaceInterpolator(
	const std::vector<double>& strikes_,
	const std::vector<double>& expiries_,
	const std::vector<std::vector<double>>& vols_)
	: strikes(strikes_), expiries(expiries_), impliedVols(vols_) {

	if (strikes.empty() || expiries.empty() || impliedVols.empty())
		throw std::runtime_error("Invalid input grid");

	if (impliedVols.size() != strikes.size() || impliedVols[0].size() != expiries.size())
		throw std::runtime_error("Vol grid size mismatch");
}

double VolSurfaceInterpolator::CubicInterp1D(const std::vector<double>& x, const std::vector<double>& y, double x0) const {
	size_t n = x.size();

	if (n < 2)
		throw std::runtime_error("At least two points required for interpolation/extrapolation.");

	// Left extrapolation
	if (x0 <= x.front())
		return y[0] + (x0 - x[0]) * (y[1] - y[0]) / (x[1] - x[0]);

	// Right extrapolation
	if (x0 >= x.back())
		return y[n - 2] + (x0 - x[n - 2]) * (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]);

	// Now find interval and do cubic interpolation
	size_t i = 1;
	while (i + 2 < n && x0 > x[i + 1])
		++i;
	i = std::max<size_t>(1, std::min(i, n - 3));

	double x1 = x[i - 1], x2 = x[i], x3 = x[i + 1], x4 = x[i + 2];
	double y1 = y[i - 1], y2 = y[i], y3 = y[i + 1], y4 = y[i + 2];

	double t = (x0 - x2) / (x3 - x2);

	double a0 = y2;
	double a1 = -0.5 * y1 + 0.5 * y3;
	double a2 = y1 - 2.5 * y2 + 2 * y3 - 0.5 * y4;
	double a3 = -0.5 * y1 + 1.5 * y2 - 1.5 * y3 + 0.5 * y4;

	return ((a3 * t + a2) * t + a1) * t + a0;
}


double VolSurfaceInterpolator::BicubicInterp(double K, double T) const {
	std::vector<double> tempVols;

	for (const auto& volRow : impliedVols)
		tempVols.push_back(CubicInterp1D(expiries, volRow, T));

	return CubicInterp1D(strikes, tempVols, K);
}

double VolSurfaceInterpolator::GetImpliedVol(double K, double T) const {
	return BicubicInterp(K, T);
}
