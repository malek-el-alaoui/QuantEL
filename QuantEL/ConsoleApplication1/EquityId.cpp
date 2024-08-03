#include "EquityId.h"

EquityId::EquityId(const std::string id) : identifier(id) {};
EquityId::EquityId(const EquityId& id) : identifier(id.identifier) {};
