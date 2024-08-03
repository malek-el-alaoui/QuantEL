#pragma once
#include "MarketQuote.h"
#include "EquityId.h"
#include <iostream>

class EquitySmileQuote : MarketQuote {
public:
	EquitySmileQuote(double strike, double impliedVol);
	EquitySmileQuote(const std ::string& id, double strike, double impliedVol);
	EquitySmileQuote(const EquitySmileQuote &c);
	EquitySmileQuote& operator=(const EquitySmileQuote &c);
	double value() const override;
	bool isValid() const override;
private:
	double strike;
	double impliedVol;

};

