// Copyright 2018 Romanov Alexander

#ifndef ____BLACK_SHOLES_MARKET_MODEL____
#define ____BLACK_SHOLES_MARKET_MODEL____

#include "mkl.h"
#include <math.h>

#if defined(_WIN64) || defined(_WIN32) 
	#pragma warning(disable : 4005) // fopen warning disable
	#pragma warning(disable : 4996) // REDEFINITION WARNING DISABLED!
#endif


#define TIME		3.0f	// option execute time (years)
#define SIG			0.2f	// volatility; percent per year 0.2 -> 20%
#define R			0.05f	// the interest rate; percent per year 0.05 -> 5%		
#define S0			100.0f	// option price at t == 0

#define __SEED__	20000000ll


class BlackScholes {
	// todo: add ReadMe
public:
	BlackScholes() {};
	~BlackScholes() {};
	float GetStockPrice(float wiener_diff, float time);
	virtual void Execute() = 0;
	virtual void TGetOptionPrice();
};

#endif // !____BLACK_SHOLES_MARKET_MODEL____
