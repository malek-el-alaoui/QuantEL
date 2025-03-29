#include "ComplexMonteCarloPricer.h"
#include <random>
#include <cmath>

ComplexMonteCarloPricer::ComplexMonteCarloPricer(std::shared_ptr<ComplexFunction2D> sigmaLoc_,
	std::shared_ptr<ComplexFunction> payoff_,
	Complex z0_,
	double r_,
	double T_,
	int steps_,
	int paths_)
	: sigmaLoc(sigmaLoc_), payoff(payoff_), z0(z0_),
	r(r_), T(T_), steps(steps_), paths(paths_) {}

Complex ComplexMonteCarloPricer::Price() {
	std::mt19937 gen(std::random_device{}());
	std::normal_distribution<> norm(0.0, 1.0);

	double dt = T / steps;
	Complex sum = Complex(0.0);

	for (int p = 0; p < paths; ++p) {
		Complex z = z0;
		double t = 0.0;

		for (int i = 0; i < steps; ++i) {
			Complex sigma = sigmaLoc->Evaluate(z, t);
			double dW = norm(gen) * std::sqrt(dt);

			// Euler-Maruyama step for complex SDE
			Complex drift = Complex(r * dt);
			Complex diffusion = sigma * dW;

			z = z * (Complex(1.0) + drift + diffusion);
			t += dt;
		}

		sum = sum + payoff->Evaluate(z);
	}

	return sum / double(paths);  // No discounting since r is inside drift
}
