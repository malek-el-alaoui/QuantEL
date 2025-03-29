#ifndef LOCAL_VOLATILITY_FUNCTION_HPP
#define LOCAL_VOLATILITY_FUNCTION_HPP

#include "CallSurfaceFunction.h"

class LocalVolatilityFunction : public ComplexFunction2D {
private:
	std::shared_ptr<RealFunction2D> callSurface;
	double dK; // small shift in strike
	double dt; // small shift in time

public:
	LocalVolatilityFunction(std::shared_ptr<RealFunction2D> callSurface, double dK = 1e-2, double dt = 1e-3);

	Complex Evaluate(double strike, double t) const override;
};

#endif
#pragma once
