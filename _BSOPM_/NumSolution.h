// Copyright 2018 Romanov Alexander
#include "BSOPM.h"
#include "E_RowTable.h"
#ifndef ____BLACK_SHOLES_NUMERICAL_SOLUTION____
#define ____BLACK_SHOLES_NUMERICAL_SOLUTION____

#define M_PIF 3.14159265358979323846f

class NumSolution : public BSOPM {
public:
	NumSolution() {};
	~NumSolution() {};
	bool checkConvergence(int StepIndex, VSLStreamStatePtr stream, int npaths, int nsteps, float pS0, float pR, float pSig, float time, float *error, unsigned int seed, int indexGen);
	void writeMethodConvergence(int StepIndex, int indexGen, int npaths, int nsteps, float pS0, float pR, float pSig, float time, unsigned int seed);
	void getErrors(int nsteps, int indexGen, int N, unsigned int seed, float K, float R, float Time, float SIG, float pS0, int sampleStep, float fair);
	void MCParExecute(int StepIndex, int nsteps, int indexGen, int N, unsigned int seed, float K, float R, float Time, float SIG, float pS0);
	float getMCPricePar(int NumThreads, int StepIndex, int nsteps, int indexGen, int N, unsigned int seed, float K, float R, float Time, float SIG, float S0, double& workTime);
	float SimulateStockPrices(int StepIndex, int indexGen, int npaths, int nsteps, float pS0, float pR, float pSig, float time, unsigned int seed);
	float SimulateStockPricesVol(int StepIndex, int npaths, int nsteps, float pS0, float* pR, float* pSig, float time, unsigned int seed, int indexGen);
	float getMCPrice(int StepIndex, int nsteps, int indexGen, int N, unsigned int seed, float K, float R, float Time, float SIG, float S0);
	void WriteMethodErrors(float* Errors, int nsteps, int nrows, float Time, int scale, int stepIndex);
	
	void GetRPricePar	(float a, float b, int scale, int NumThreads, int N, float r, float sig, float* pT, float* pK, float* pS0, float* pC);
	void GetTPricePar	(float a, float b, int scale, int NumThreads, int N, float r, float sig, float* pT, float* pK, float* pS0, float* pC);
	void GetSPricePar	(float a, float b, int scale, int NumThreads, int N, float r, float sig, float* pT, float* pK, float* pS0, float* pC);
	void Get3_8PricePar	(float a, float b, int scale, int NumThreads, int N, float r, float sig, float* pT, float* pK, float* pS0, float* pC);
	void GetSPrice(float a, float b, int scale, int N, float r, float sig, float* pT, float* pK, float* pS0, float* pC);

protected:
	clock_t start, finish;
	double t;
	void wienerArrayLarge(VSLStreamStatePtr stream, int nsteps, int nsamples, float time, float *w);
	void wAndZArrayLarge(VSLStreamStatePtr stream, int nsteps, int nsamples, float time, float *w);
	void wAndZProcesses(VSLStreamStatePtr stream, int nsteps, float time, float *buffer);
	float* memoryWrajectAlloc(VSLStreamStatePtr stream, int StepIndex, int npaths, int nsteps, float time, float* wtraject);
	float Integrand(float z, float pS0, float K, float tmp1, float tmp2, int scale);
	float EulMarStep(float S, float dt, float dw, float r, float sig);
	float MilsteinStep(float S, float dt, float dw, float r, float sig);
	float RK1Step(float S, float dt, float dw, float r, float sig);
	float BurragePlatenStep(float S, float dt, float dw, float dz, float r, float sig);
	float stockPricesIntegrator(VSLStreamStatePtr stream, float* wtraject, int StepIndex, int nsteps, float pS0, float pR, float pSig, float time);
	float stockPricesIntegratorVol(VSLStreamStatePtr stream, float* wtraject, int StepIndex, int nsteps, float pS0, float* pR, float* pSig, float time);
};

#endif // !____BLACK_SHOLES_NUMERICAL_SOLUTION____
