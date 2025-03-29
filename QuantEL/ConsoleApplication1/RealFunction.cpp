#include "RealFunction.h"

SumFunction::SumFunction(std::unique_ptr<RealFunction> f, std::unique_ptr<RealFunction> g)
	: f(std::move(f)), g(std::move(g)) {}

double SumFunction::Evaluate(double x) const {
	return f->Evaluate(x) + g->Evaluate(x);
}

ProductFunction::ProductFunction(std::unique_ptr<RealFunction> f, std::unique_ptr<RealFunction> g)
	: f(std::move(f)), g(std::move(g)) {}

double ProductFunction::Evaluate(double x) const {
	return f->Evaluate(x) * g->Evaluate(x);
}

DivFunction::DivFunction(std::unique_ptr<RealFunction> f, std::unique_ptr<RealFunction> g)
	: f(std::move(f)), g(std::move(g)) {}

double DivFunction::Evaluate(double x) const {
	double denominator = g->Evaluate(x);
	if (denominator == 0) {
		throw std::exception("Division by zero");
	}
	return f->Evaluate(x) / denominator;
}

ConstFunction::ConstFunction(double c) : c(c) {}

double ConstFunction::Evaluate(double x) const {
	return c;
}

std::unique_ptr<RealFunction> operator+(std::unique_ptr<RealFunction> f, std::unique_ptr<RealFunction> g) {
	return std::make_unique<SumFunction>(std::move(f), std::move(g));
}

std::unique_ptr<RealFunction> operator+(std::unique_ptr<RealFunction> f, double a) {
	return std::make_unique<SumFunction>(std::move(f), std::make_unique<ConstFunction>(a));
}

std::unique_ptr<RealFunction> operator+(double a, std::unique_ptr<RealFunction> f) {
	return std::make_unique<SumFunction>(std::make_unique<ConstFunction>(a), std::move(f));
}

std::unique_ptr<RealFunction> operator*(std::unique_ptr<RealFunction> f, std::unique_ptr<RealFunction> g) {
	return std::make_unique<ProductFunction>(std::move(f), std::move(g));
}

std::unique_ptr<RealFunction> operator*(std::unique_ptr<RealFunction> f, double a) {
	return std::make_unique<ProductFunction>(std::move(f), std::make_unique<ConstFunction>(a));
}

std::unique_ptr<RealFunction> operator*(double a, std::unique_ptr<RealFunction> f) {
	return std::make_unique<ProductFunction>(std::make_unique<ConstFunction>(a), std::move(f));
}

std::unique_ptr<RealFunction> operator/(std::unique_ptr<RealFunction> f, std::unique_ptr<RealFunction> g) {
	return std::make_unique<DivFunction>(std::move(f), std::move(g));
}

std::unique_ptr<RealFunction> operator/(std::unique_ptr<RealFunction> f, double a) {
	return std::make_unique<DivFunction>(std::move(f), std::make_unique<ConstFunction>(a));
}

std::unique_ptr<RealFunction> operator/(double a, std::unique_ptr<RealFunction> f) {
	return std::make_unique<DivFunction>(std::make_unique<ConstFunction>(a), std::move(f));
}
