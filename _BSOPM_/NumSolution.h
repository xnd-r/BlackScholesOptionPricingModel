// Copyright 2018 Romanov Alexander
#include "BSOPM.h"
#include "E_RowTable.h"
#ifndef ____BLACK_SHOLES_NUMERICAL_SOLUTION____
#define ____BLACK_SHOLES_NUMERICAL_SOLUTION____


class NumSolution : public BSOPM {
public:
	NumSolution() {};
	~NumSolution() {};
	bool IsConvergence(int StepIndex, VSLStreamStatePtr stream, int npaths, int nsteps, float pS0, float pR, float pSig, float time, float *error);
	float SimulateStockPrices(int StepIndex, int npaths, int nsteps, float pS0, float pR, float pSig, float time);
	void WriteToCsv(float* Errors, int nsteps, int nrows, float time, int scale, int stepIndex);
	void Execute(int StepIndex, int npaths, int nsteps, float pS0, float pR, float pSig, float time);
	float SimulateStockPricesVol(int StepIndex, int npaths, int nsteps, float pS0, float* pR, float* pSig, float time);
protected:
	float EulMarStep(float S, float dt, float dw, float r, float sig);
	float MilsteinStep(float S, float dt, float dw, float r, float sig);
	float RK1Step(float S, float dt, float dw, float r, float sig);
	float BurragePlatenStep(float S, float dt, float dw, float dz, float r, float sig);
	void wAndZProcesses(VSLStreamStatePtr stream, int nsteps, float time, float *buffer);

};

#endif // !____BLACK_SHOLES_NUMERICAL_SOLUTION____
