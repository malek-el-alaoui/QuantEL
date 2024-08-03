#pragma once
class MarketQuote
{
public:
	~MarketQuote() = default;
	//! returns the current value
	virtual double value() const = 0;
	//! returns true if the Quote holds a valid value
	virtual bool isValid() const = 0;
};


