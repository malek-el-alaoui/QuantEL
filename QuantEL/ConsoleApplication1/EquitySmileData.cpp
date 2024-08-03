#include "EquitySmileData.h"

EquitySmileData::EquitySmileData(EquityId& id) : equityId(id) {

}
void EquitySmileData :: addQuote(const EquitySmileQuote& quote) {
	quotes.push_back(quote);
}

const std::vector<EquitySmileQuote>& EquitySmileData::getQuotes() const {
	return quotes;
}

const EquityId& EquitySmileData :: getEquityId() const {
	return equityId;
}
