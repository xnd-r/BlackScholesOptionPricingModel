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
	void WriteMethodErrors(float* Errors, int nsteps, int nrows, float time, int scale, int stepIndex);
	void Execute(int StepIndex, int indexGen, int npaths, int nsteps, float pS0, float pR, float pSig, float time, unsigned int seed);
	void getErrors(int nsteps, int indexGen, int N, unsigned int seed, float K, float R, float Time, float SIG, float S0,int sampleStep, float fair);
	void MCParExecute(int StepIndex, int nsteps, int indexGen, int N, unsigned int seed, float K, float R, float Time, float SIG, float pS0);
	void Execute();
	float getMCPricePar(int NumThreads, int StepIndex, int nsteps, int indexGen, int N, unsigned int seed, float K, float R, float Time, float SIG, float S0, double& workTime);
	float SimulateStockPrices(int StepIndex, int indexGen, int npaths, int nsteps, float pS0, float pR, float pSig, float time, unsigned int seed);
	float SimulateStockPricesVol(int StepIndex, int npaths, int nsteps, float pS0, float* pR, float* pSig, float time, unsigned int seed, int indexGen);
	
	float getMCPrice(int StepIndex, int nsteps, int indexGen, int N, unsigned int seed, float K, float R, float Time, float SIG, float S0);
	float GetRPrice		(float a, float b);
	float GetTPrice		(float a, float b);
	float GetSPrice		(float a, float b);
	float Get3_8Price	(float a, float b);
	float GetRPrice		(float a, float b, int NumThread/*, float* s_array, float* expf_array*/);
	float GetTPrice		(float a, float b, int NumThread);
	float GetSPrice		(float a, float b, int NumThread);
	float Get3_8Price	(float a, float b, int NumThread);
	float GetRPrice(float a, float b, int StepIndex, int npaths, int nsteps, float pS0, float pR, float pSig, float time, float K);
protected:
	clock_t start, finish;
	double t;
	void SetS(int amo);
	void wienerProcessLarge(VSLStreamStatePtr stream, int nsteps, int nsamples, float time, float *w);
	void wAndZProcessesLarge(VSLStreamStatePtr stream, int nsteps, int nsamples, float time, float *w);
	void wAndZProcesses(VSLStreamStatePtr stream, int nsteps, float time, float *buffer);
	float* memoryWrajectAlloc(VSLStreamStatePtr stream, int StepIndex, int npaths, int nsteps, float time, float* wtraject);
	float* s_array;
	float* exp_array;
	float Integrand(float x);
	float EulMarStep(float S, float dt, float dw, float r, float sig);
	float MilsteinStep(float S, float dt, float dw, float r, float sig);
	float RK1Step(float S, float dt, float dw, float r, float sig);
	float BurragePlatenStep(float S, float dt, float dw, float dz, float r, float sig);
	float stockPricesIntegrator(VSLStreamStatePtr stream, float* wtraject, int StepIndex, int nsteps, float pS0, float pR, float pSig, float time);
	float stockPricesIntegratorVol(VSLStreamStatePtr stream, float* wtraject, int StepIndex, int nsteps, float pS0, float* pR, float* pSig, float time);
};

#endif // !____BLACK_SHOLES_NUMERICAL_SOLUTION____
