#ifndef COMPLEX_MONTE_CARLO_PRICER_HPP
#define COMPLEX_MONTE_CARLO_PRICER_HPP

#include "ComplexFramework.h"
#include <memory>

class ComplexMonteCarloPricer {
private:
	std::shared_ptr<ComplexFunction2D> sigmaLoc;
	std::shared_ptr<ComplexFunction> payoff;
	Complex z0;
	double r;
	double T;
	int steps;
	int paths;

public:
	ComplexMonteCarloPricer(std::shared_ptr<ComplexFunction2D> sigmaLoc,
		std::shared_ptr<ComplexFunction> payoff,
		Complex z0,
		double r,
		double T,
		int steps,
		int paths);

	Complex Price();
};

#endif
#pragma once
