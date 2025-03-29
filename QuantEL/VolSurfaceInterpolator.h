#ifndef VOL_SURFACE_INTERPOLATOR_HPP
#define VOL_SURFACE_INTERPOLATOR_HPP

#include <vector>

class VolSurfaceInterpolator {
private:
	std::vector<double> strikes;
	std::vector<double> expiries;
	std::vector<std::vector<double>> impliedVols; // vol[i][j] = vol at (strike[i], expiry[j])

public:
	VolSurfaceInterpolator(
		const std::vector<double>& strikes,
		const std::vector<double>& expiries,
		const std::vector<std::vector<double>>& vols);

	double GetImpliedVol(double K, double T) const;

private:
	double CubicInterp1D(const std::vector<double>& x, const std::vector<double>& y, double x0) const;
	double BicubicInterp(double K, double T) const;
};

#endif
#pragma once
