#include "LocalVolatilityFunction.h"
#include <cmath>
#include <stdexcept>

LocalVolatilityFunction::LocalVolatilityFunction(std::shared_ptr<RealFunction2D> callSurface_, double dK_, double dt_)
	: callSurface(callSurface_), dK(dK_), dt(dt_) {}

Complex LocalVolatilityFunction::Evaluate(double K, double T) const {
	if (K <= 0 || T <= 0)
		return 0.0;

	// First time derivative
	double Ct_plus = callSurface->Eval(K, T + dt);
	double Ct_minus = callSurface->Eval(K, T - dt);
	double dC_dt = (Ct_plus - Ct_minus) / (2.0 * dt);

	// Second strike derivative
	double C_up = callSurface->Eval(K + dK, T);
	double C_mid = callSurface->Eval(K, T);
	double C_dn = callSurface->Eval(K - dK, T);

	double d2C_dK2 = (C_up - 2 * C_mid + C_dn) / (dK * dK);

	// Dupire formula
	double denom = 0.5 * K * K * d2C_dK2;
	if (denom <= 1e-10)
		return 0.0;

	double localVar = dC_dt / denom;
	if (localVar <= 0.0)
		return 0.0;

	return Complex(dC_dt / denom).sqrt();
}
