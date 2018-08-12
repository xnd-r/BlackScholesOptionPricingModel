// Copyright 2018 Romanov Alexander

#ifndef ____STOCK_PTIECE____
#define ____STOCK_PTIECE____

#include <time.h>
#include <fstream>
#include <algorithm>
#include <string>
#include <vector>

#include "../../../include/B00_BlackScholes.h"
#include "../../../etc/RowTable/E_RowTable.h"

#define NPATHS		50	// amo of trajectories
#define NSTEPS		300 

class StockPrice : public BlackScholes {
	
public:
	std::vector<float> WienerTraject;
	VSLStreamStatePtr InitGen();
	void FreeGen(VSLStreamStatePtr stream);
	__declspec(noinline) void GenerateGauss(float expect_val, float deviation, int amou, VSLStreamStatePtr stream, float *dest_array);
	__declspec(noinline) virtual void SimulateWienerProcess(int nPaths, int nSteps, float Time, float **buffer);
	virtual void WriteToCsv(float **buffer, int nRows, int nColumns);
	void Execute();
};

#endif // !____STOCK_PTIECE____