#ifndef CALL_SURFACE_FUNCTION_HPP
#define CALL_SURFACE_FUNCTION_HPP

#include "VolSurfaceInterpolator.h"
#include "ComplexFramework.h"

class RealFunction2D : public ComplexFunction2D {
public:
	virtual double Eval(double strike, double t) const = 0;

	Complex Evaluate(const Complex& strike, const Complex& t) const override {
		return Complex(Evaluate(strike.Re(), t.Re()));
	}

	Complex Evaluate(const Complex& strike, double t) const override {
		return Complex(Evaluate(strike.Re(), t));
	}
	Complex Evaluate(double strike, double t) const override {
		return Complex(Eval(strike, t));
	}

	virtual ~RealFunction2D() = default;
};

class CallSurfaceFunction : public RealFunction2D {
private:
	VolSurfaceInterpolator volSurface;
	double spot;
	double rate;

public:
	CallSurfaceFunction(const VolSurfaceInterpolator& surface, double spot, double rate);
	double Eval(double strike, double t) const override;
};

#endif
