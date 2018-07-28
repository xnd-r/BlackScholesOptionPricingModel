// Copyright 2018 Romanov Alexander

#ifndef ____BLACK_SHOLES_MARKET_MODEL____
#define ____BLACK_SHOLES_MARKET_MODEL____

#include "mkl.h"

// TODO: check which defines should be moved
// TODO: add paths to the properties

#define TIME		3.0		// option execute time (years)
#define NPATHS		50	// amo of trajectories
#define PART_PATH	20 // amo of paths to build avg trajectory of stock priece

#define SIG			0.2		// volatility; percent per year 0.2 -> 20%
#define R			0.05	// the interest rate; percent per year 0.05 -> 5%		
#define S0			100.0	// option price at t == 0

#define NSTEPS		300 
#define __SEED__	20000000	


class BlackSholes {
	// todo: add ReadMe
	// Abstract/virtual class?
	// constructor/destr ?
public:
	//BlackSholes();
	//~BlackSholes();
	virtual void GetStockPriece();
	virtual void GetOptPriece();
};

#endif // !____BLACK_SHOLES_MARKET_MODEL____