// RealFunction.h
#ifndef REALFUNCTION_H
#define REALFUNCTION_H

#include <memory>

class RealFunction {
public:
	virtual double Evaluate(double x) const = 0;
	virtual ~RealFunction() {}
};

class SumFunction : public RealFunction {
private:
	std::unique_ptr<RealFunction> f;
	std::unique_ptr<RealFunction> g;

public:
	SumFunction(std::unique_ptr<RealFunction> f, std::unique_ptr<RealFunction> g);
	double Evaluate(double x) const override;
};


class ProductFunction : public RealFunction {
private:
	std::unique_ptr<RealFunction> f;
	std::unique_ptr<RealFunction> g;

public:
	ProductFunction(std::unique_ptr<RealFunction> f, std::unique_ptr<RealFunction> g);
	double Evaluate(double x) const override;
};

class DivFunction : public RealFunction {
private:
	std::unique_ptr<RealFunction> f;
	std::unique_ptr<RealFunction> g;

public:
	DivFunction(std::unique_ptr<RealFunction> f, std::unique_ptr<RealFunction> g);
	double Evaluate(double x) const override;
};

class ConstFunction : public RealFunction {
private:
	double c;

public:
	ConstFunction(double c);
	double Evaluate(double x) const override;
};



std::unique_ptr<RealFunction> operator+(std::unique_ptr<RealFunction> f, std::unique_ptr<RealFunction> g);
std::unique_ptr<RealFunction> operator+(std::unique_ptr<RealFunction> f, double a);
std::unique_ptr<RealFunction> operator+(double a, std::unique_ptr<RealFunction> f);

std::unique_ptr<RealFunction> operator*(std::unique_ptr<RealFunction> f, std::unique_ptr<RealFunction> g);
std::unique_ptr<RealFunction> operator*(std::unique_ptr<RealFunction> f, double a);
std::unique_ptr<RealFunction> operator*(double a, std::unique_ptr<RealFunction> f);

std::unique_ptr<RealFunction> operator/(std::unique_ptr<RealFunction> f, std::unique_ptr<RealFunction> g);
std::unique_ptr<RealFunction> operator/(std::unique_ptr<RealFunction> f, double a);
std::unique_ptr<RealFunction> operator/(double a, std::unique_ptr<RealFunction> f);

#endif // REALFUNCTION_H
