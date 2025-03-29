#include "EquitySmileQuote.h"
#include <iostream>
#include <cmath> // For std::isnan, std::isinf, and std::isfinite

EquitySmileQuote::EquitySmileQuote(double k, double v) : strike(k), impliedVol(v){}

EquitySmileQuote::EquitySmileQuote(double strike, double impliedVol): strike(strike), impliedVol(impliedVol) {}

EquitySmileQuote::EquitySmileQuote(const EquitySmileQuote &c) : strike(c.strike), impliedVol(c.impliedVol) {}

EquitySmileQuote&::EquitySmileQuote::operator=(const EquitySmileQuote &c) {
	if (this == &c) {
		return *this;
	}
	strike = c.strike;
	impliedVol = c.impliedVol;

	return *this;
}

double ::EquitySmileQuote::value() const {
	return impliedVol;
}
bool ::EquitySmileQuote::isValid() const
{
	if (std::isnan(impliedVol) || std::isnan(strike))
	{
		return false;
	}
	return true;
}