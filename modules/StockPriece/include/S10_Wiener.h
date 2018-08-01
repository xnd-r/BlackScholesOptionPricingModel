// Copyright 2018 Romanov Alexander

#ifndef ____STOCK_PTIECE____
#define ____STOCK_PTIECE____

#include <time.h>
#include <fstream>
#include <algorithm>
#include <string>

#include "../../../include/BlackScholes.h"
#include "../../../etc/RowTable/RowTable.h"

#define NPATHS		50	// amo of trajectories
#define NSTEPS		300 

class StockPrice : public BlackSñholes {

public:

	VSLStreamStatePtr InitGen();
	void FreeGen(VSLStreamStatePtr stream);
	void GenerateGauss(double expect_val, double deviation, int amou, VSLStreamStatePtr stream, double *dest_array);
	
	virtual void SimulateWienerProcess(int nPaths, int nSteps, double Time, double **buffer);
	virtual void WriteToCsv(double **buffer, int nRows, int nColumns);
	void Execute();
};

#endif // !____STOCK_PTIECE____