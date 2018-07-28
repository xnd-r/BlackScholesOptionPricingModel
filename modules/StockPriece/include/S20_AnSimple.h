#ifndef ____SDE_ANALITYCAL_EQUATION_SIMPLE____
#define ____SDE_ANALITYCAL_EQUATION_SIMPLE____

#include "S10_Wiener.h"

class AnSimple : public StockPriece {
#define NPATHS 10000
public:
	virtual void WriteToCsv(double *buffer, int nRows, int nColumns);
	double GetStockPrice(double z, double t);
	double SimulateStockPrices(int nPaths, double Time, double *sBuffer);
};

#endif // !____SDE_ANALITYCAL_EQUATION_SIMPLE____