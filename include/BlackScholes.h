// Copyright 2018 Romanov Alexander

#ifndef ____BLACK_SHOLES_MARKET_MODEL____
#define ____BLACK_SHOLES_MARKET_MODEL____

#include "mkl.h"
#include <math.h>

#pragma warning(disable : 4005) // fopen warning disable
#pragma warning(disable : 4996) // REDEFINITION WARNING DISABLED!

#define TIME		3.0		// option execute time (years)
#define SIG			0.2		// volatility; percent per year 0.2 -> 20%
#define R			0.05	// the interest rate; percent per year 0.05 -> 5%		
#define S0			100.0	// option price at t == 0

#define __SEED__	20000000	



class BlackSñholes {
	// todo: add ReadMe
public:
	//BlackSholes();
	//~BlackSholes();
	double GetStockPrice(double wiener_diff, double time);
	virtual void Execute() = 0;
	//virtual void GetOptionPrice();
};

#endif // !____BLACK_SHOLES_MARKET_MODEL____
