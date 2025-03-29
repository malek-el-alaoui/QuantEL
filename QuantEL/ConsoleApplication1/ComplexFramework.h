#ifndef COMPLEX_FRAMEWORK_HPP
#define COMPLEX_FRAMEWORK_HPP

#include <iostream>
#include <cmath>
#include <memory>

// ------------------ Complex ------------------
class Complex {
private:
	double real, imag;

public:
	Complex(double r = 0.0, double i = 0.0);
	double Re() const;
	double Im() const;

	Complex operator+(const Complex& other) const;
	Complex operator-(const Complex& other) const;
	Complex operator*(const Complex& other) const;
	Complex operator/(const Complex& other) const;

	Complex conj() const;
	double abs() const;
	double arg() const;
	Complex sqrt() const;
	Complex exp() const;
	Complex log() const;
	Complex sin() const;
	Complex cos() const;
	static Complex erf(const Complex& z);
	static Complex norm_cdf(const Complex& z);
	friend std::ostream& operator<<(std::ostream& os, const Complex& c);
};

// ------------------ Function 1D ------------------
class ComplexFunction {
public:
	virtual Complex Evaluate(const Complex& x) const = 0;
	virtual ~ComplexFunction() = default;
};

using FunctionPtr = std::shared_ptr<ComplexFunction>;

class ConstFunction : public ComplexFunction {
	Complex c;
public:
	ConstFunction(Complex c);
	Complex Evaluate(const Complex& x) const override;
};

class IdentityFunction : public ComplexFunction {
public:
	Complex Evaluate(const Complex& x) const override;
};

class ExpFunction : public ComplexFunction {
public:
	Complex Evaluate(const Complex& x) const override;
};

class SinFunction : public ComplexFunction {
public:
	Complex Evaluate(const Complex& x) const override;
};

class AddFunction : public ComplexFunction {
	FunctionPtr f, g;
public:
	AddFunction(FunctionPtr f, FunctionPtr g);
	Complex Evaluate(const Complex& x) const override;
};

class MulFunction : public ComplexFunction {
	FunctionPtr f, g;
public:
	MulFunction(FunctionPtr f, FunctionPtr g);
	Complex Evaluate(const Complex& x) const override;
};

FunctionPtr operator+(FunctionPtr f, FunctionPtr g);
FunctionPtr operator*(FunctionPtr f, FunctionPtr g);

// ------------------ Function 2D ------------------
class ComplexFunction2D {
public:
	virtual Complex Evaluate(const Complex& z1, const Complex& z2) const = 0;

	// Real-valued second argument
	virtual Complex Evaluate(const Complex& z1, double t) const {
		return Evaluate(z1, Complex(t, 0.0));
	}

	// Real-valued arguments
	virtual Complex Evaluate(double x, double t) const {
		return Evaluate(Complex(x, 0.0), Complex(t, 0.0));
	}

	virtual ~ComplexFunction2D() = default;
};

using Function2DPtr = std::shared_ptr<ComplexFunction2D>;

class Add2DFunction : public ComplexFunction2D {
public:
	Complex Evaluate(const Complex& z1, const Complex& z2) const override;
};

class ExpMulFunction : public ComplexFunction2D {
public:
	Complex Evaluate(const Complex& z1, const Complex& z2) const override;
};

class RealSecondArgFunction : public ComplexFunction2D {
	Function2DPtr inner;
public:
	RealSecondArgFunction(Function2DPtr f);
	Complex Evaluate(const Complex& z1, const Complex& z2) const override;
};

// Composite 2D function operators
class AddFunction2D : public ComplexFunction2D {
	Function2DPtr f, g;
public:
	AddFunction2D(Function2DPtr f, Function2DPtr g);
	Complex Evaluate(const Complex& z1, const Complex& z2) const override;
};

class SubFunction2D : public ComplexFunction2D {
	Function2DPtr f, g;
public:
	SubFunction2D(Function2DPtr f, Function2DPtr g);
	Complex Evaluate(const Complex& z1, const Complex& z2) const override;
};

class MulFunction2D : public ComplexFunction2D {
	Function2DPtr f, g;
public:
	MulFunction2D(Function2DPtr f, Function2DPtr g);
	Complex Evaluate(const Complex& z1, const Complex& z2) const override;
};

class DivFunction2D : public ComplexFunction2D {
	Function2DPtr f, g;
public:
	DivFunction2D(Function2DPtr f, Function2DPtr g);
	Complex Evaluate(const Complex& z1, const Complex& z2) const override;
};

Function2DPtr operator+(Function2DPtr f, Function2DPtr g);
Function2DPtr operator-(Function2DPtr f, Function2DPtr g);
Function2DPtr operator*(Function2DPtr f, Function2DPtr g);
Function2DPtr operator/(Function2DPtr f, Function2DPtr g);

#endif
