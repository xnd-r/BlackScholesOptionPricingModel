#ifndef ____SDE_ANALITYCAL_EQUATION_SIMPLE____
#define ____SDE_ANALITYCAL_EQUATION_SIMPLE____

#include "S10_Wiener.h"

class AnSimple : public StockPrice {
#define NPATHS 10000
public:
	virtual void WriteToCsv(float *buffer, int nRows, int nColumns, float avg);
	__declspec(noinline) float SimulateStockPrices(int nPaths, float Time, float *sBuffer);
	void Execute();
};

#endif // !____SDE_ANALITYCAL_EQUATION_SIMPLE____