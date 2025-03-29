#ifndef ANALYTIC_PAYOFF_HPP
#define ANALYTIC_PAYOFF_HPP

#include "ComplexFramework.h"
#include <cmath>

class ComplexCallPayoff : public ComplexFunction {
private:
	double K;
	double sigma;

public:
	ComplexCallPayoff(double K_, double sigma_);
	Complex Evaluate(const Complex& z) const override;
};

#endif
