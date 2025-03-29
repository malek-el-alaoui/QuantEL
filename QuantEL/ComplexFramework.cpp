#include "ComplexFramework.h"

// ------------------ Complex ------------------
Complex::Complex(double r, double i) : real(r), imag(i) {}
double Complex::Re() const { return real; }
double Complex::Im() const { return imag; }

Complex Complex::operator+(const Complex& o) const { return Complex(real + o.real, imag + o.imag); }
Complex Complex::operator-(const Complex& o) const { return Complex(real - o.real, imag - o.imag); }
Complex Complex::operator*(const Complex& o) const {
	return Complex(real * o.real - imag * o.imag, real * o.imag + imag * o.real);
}
Complex Complex::operator/(const Complex& o) const {
	double denom = o.real * o.real + o.imag * o.imag;
	return Complex((real * o.real + imag * o.imag) / denom, (imag * o.real - real * o.imag) / denom);
}

Complex Complex::conj() const { return Complex(real, -imag); }
double Complex::abs() const { return std::sqrt(real * real + imag * imag); }
double Complex::arg() const { return std::atan2(imag, real); }
Complex Complex::sqrt() const {
	double r = this->abs();
	double a = this->Re();
	double b = this->Im();

	double realPart = std::sqrt((r + a) / 2.0);
	double imagPart = std::sqrt((r - a) / 2.0);
	if (b < 0) imagPart = -imagPart;

	return Complex(realPart, imagPart);
}
Complex Complex::exp() const {
	double e = std::exp(real);
	return Complex(e * std::cos(imag), e * std::sin(imag));
}
Complex Complex::log() const {
	return Complex(std::log(this->abs()), this->arg());
}
Complex Complex::sin() const {
	return Complex(std::sin(real) * std::cosh(imag), std::cos(real) * std::sinh(imag));
}
Complex Complex::cos() const {
	return Complex(std::cos(real) * std::cosh(imag), -std::sin(real) * std::sinh(imag));
}
static Complex erf(const Complex& z) {
	const int N = 30;
	const double pi = 3.14159265358979323846;
	Complex sum(0.0);
	Complex term = z;

	for (int n = 0; n < N; ++n) {
		double coeff = std::pow(-1, n) / (std::tgamma(n + 1) * (2 * n + 1));
		sum = sum + term * coeff;
		term = term * z * z;
	}

	return sum * (2.0 / std::sqrt(pi));
}

static Complex norm_cdf(const Complex& z) {
	return Complex(0.5) * (Complex(1.0) + erf(z / std::sqrt(2.0)));
}
std::ostream& operator<<(std::ostream& os, const Complex& c) {
	os << c.real;
	if (c.imag >= 0) os << " + " << c.imag << "i";
	else os << " - " << -c.imag << "i";
	return os;
}

// ------------------ Function 1D ------------------
ConstFunction::ConstFunction(Complex c_) : c(c_) {}
Complex ConstFunction::Evaluate(const Complex&) const { return c; }

Complex IdentityFunction::Evaluate(const Complex& x) const { return x; }
Complex ExpFunction::Evaluate(const Complex& x) const { return x.exp(); }
Complex SinFunction::Evaluate(const Complex& x) const { return x.sin(); }

AddFunction::AddFunction(FunctionPtr f_, FunctionPtr g_) : f(f_), g(g_) {}
Complex AddFunction::Evaluate(const Complex& x) const { return f->Evaluate(x) + g->Evaluate(x); }

MulFunction::MulFunction(FunctionPtr f_, FunctionPtr g_) : f(f_), g(g_) {}
Complex MulFunction::Evaluate(const Complex& x) const { return f->Evaluate(x) * g->Evaluate(x); }

FunctionPtr operator+(FunctionPtr f, FunctionPtr g) { return std::make_shared<AddFunction>(f, g); }
FunctionPtr operator*(FunctionPtr f, FunctionPtr g) { return std::make_shared<MulFunction>(f, g); }

// ------------------ Function 2D ------------------
Complex Add2DFunction::Evaluate(const Complex& z1, const Complex& z2) const {
	return z1 + z2;
}

Complex ExpMulFunction::Evaluate(const Complex& z1, const Complex& z2) const {
	return z1 * z2.exp();
}

RealSecondArgFunction::RealSecondArgFunction(Function2DPtr f) : inner(f) {}
Complex RealSecondArgFunction::Evaluate(const Complex& z1, const Complex& z2) const {
	return inner->Evaluate(z1, Complex(z2.Re(), 0.0));
}

// --- Composite ---
AddFunction2D::AddFunction2D(Function2DPtr f_, Function2DPtr g_) : f(f_), g(g_) {}
Complex AddFunction2D::Evaluate(const Complex& z1, const Complex& z2) const {
	return f->Evaluate(z1, z2) + g->Evaluate(z1, z2);
}

SubFunction2D::SubFunction2D(Function2DPtr f_, Function2DPtr g_) : f(f_), g(g_) {}
Complex SubFunction2D::Evaluate(const Complex& z1, const Complex& z2) const {
	return f->Evaluate(z1, z2) - g->Evaluate(z1, z2);
}

MulFunction2D::MulFunction2D(Function2DPtr f_, Function2DPtr g_) : f(f_), g(g_) {}
Complex MulFunction2D::Evaluate(const Complex& z1, const Complex& z2) const {
	return f->Evaluate(z1, z2) * g->Evaluate(z1, z2);
}

DivFunction2D::DivFunction2D(Function2DPtr f_, Function2DPtr g_) : f(f_), g(g_) {}
Complex DivFunction2D::Evaluate(const Complex& z1, const Complex& z2) const {
	return f->Evaluate(z1, z2) / g->Evaluate(z1, z2);
}

Function2DPtr operator+(Function2DPtr f, Function2DPtr g) {
	return std::make_shared<AddFunction2D>(f, g);
}
Function2DPtr operator-(Function2DPtr f, Function2DPtr g) {
	return std::make_shared<SubFunction2D>(f, g);
}
Function2DPtr operator*(Function2DPtr f, Function2DPtr g) {
	return std::make_shared<MulFunction2D>(f, g);
}
Function2DPtr operator/(Function2DPtr f, Function2DPtr g) {
	return std::make_shared<DivFunction2D>(f, g);
}
