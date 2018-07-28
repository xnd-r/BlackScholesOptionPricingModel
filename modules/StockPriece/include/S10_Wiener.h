// Copyright 2018 Romanov Alexander

#ifndef ____STOCK_PTIECE____
#define ____STOCK_PTIECE____

#include <math.h>
#include <time.h>
#include <fstream>
#include <algorithm>
#include <string>
#include "mkl.h"

#include "../../../include/BlackScholes.h"
#include "../../../etc/RowTable/RowTable.h"

#define NPATHS		50	// amo of trajectories
#define PART_PATH	20 // amo of paths to build avg trajectory of stock priece
#define NSTEPS		300 
#define __SEED__	20000000	

class StockPriece : public BlackSholes {

public:

	virtual VSLStreamStatePtr InitGen();
	virtual void FreeGen(VSLStreamStatePtr stream);
	virtual void GenerateGauss(double expect_val, double deviation, int amou, VSLStreamStatePtr stream, double *dest_array);
	virtual void SimulateWienerProcess(int nPaths, int nSteps, double Time, double **buffer);
	virtual void WriteToCsv(double **buffer, int nRows, int nColumns);
	virtual void PrintToFile(double **buffer, int nRows, int nColumns, char * FileName);
};

#endif // !____STOCK_PTIECE____