#include "ComplexCallPayoff.h"
#include "Complex.h"
#include <cmath>
#include <complex>


ComplexCallPayoff::ComplexCallPayoff(double K_, double sigma_)
	: K(K_), sigma(sigma_) {}

Complex ComplexCallPayoff::Evaluate(const Complex& z) const {
	Complex ratio = (z / K).log();
	Complex d = ratio / sigma;

	// Here, we define N(d) = norm_cdf(Re(d))
	Complex kernel = Complex::norm_cdf(d);

	return (z - K) * kernel;
}
