#pragma once
#include "MarketQuote.h"
#include "EquityId.h"

class ForwardQuote : MarketQuote
{
public:
	ForwardQuote(EquityId& id, double value);
	ForwardQuote(const ForwardQuote &c);
	ForwardQuote& operator=(const ForwardQuote &c);
	double value() const override;
	bool isValid() const override;
private:
	double _value;
};

