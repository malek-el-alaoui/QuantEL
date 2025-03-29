#pragma once
#include <iostream>

class EquityId
{
public :
	EquityId(const std ::string id);
	EquityId(const EquityId& c);
private:
	std::string identifier;
};

