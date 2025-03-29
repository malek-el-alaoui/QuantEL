#pragma once
#include <iostream>
#include <vector>
#include "EquitySmileQuote.h"
class EquitySmileData
{
public:
	EquitySmileData(EquityId& id);
	void addQuote(const EquitySmileQuote& quote);
	const std::vector<EquitySmileQuote>& getQuotes() const;
	const EquityId& getEquityId() const;

private:
	EquityId& equityId;
	std::vector<EquitySmileQuote> quotes;
};

