// Copyright 2018 Romanov Alexander

#ifndef ____EUROPEAN_OPTION____
#define ____EUROPEAN_OPTION____

#include <time.h>
#include <fstream>
#include <algorithm>
#include <string>

#include "../../../include/B00_BlackScholes.h"

#define			K	 100.0f // strike price -- price fixed in option
#define	 INVSQRT2	 0.707106781f
#define			N	1000000 //amount of options 

class EuropeanOption : public BlackScholes {
public:
	// TODO: make abstract 
	//virtual float GetMCPrice();
};

#endif // !____EUROPEAN_OPTION____