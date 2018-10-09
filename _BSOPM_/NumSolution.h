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
	bool IsConvergence(int StepIndex, VSLStreamStatePtr stream, int npaths, int nsteps, float pS0, float pR, float pSig, float time, float *error);
	float stockPricesIntegrator(VSLStreamStatePtr stream, float* wtraject, int StepIndex, int nsteps, float pS0, float pR, float pSig, float time);
	void WriteToCsv(float* Errors, int nsteps, int nrows, float time, int scale, int stepIndex);
	void Execute(int StepIndex, int npaths, int nsteps, float pS0, float pR, float pSig, float time);
	float SimulateStockPrices(int StepIndex, int npaths, int nsteps, float pS0, float pR, float pSig, float time);
	float SimulateStockPricesVol(int StepIndex, int npaths, int nsteps, float pS0, float* pR, float* pSig, float time);
	float getMCPrice(int StepIndex, int nsteps, int indexGen, int N, unsigned int seed, float K, float R, float Time, float SIG, float S0);
	float getMCPricePar(int StepIndex, int nsteps, int indexGen, int N, unsigned int seed, float K, float R, float Time, float SIG, float S0);
	void  getErrors(int nsteps, int indexGen, int N, unsigned int seed, float K, float R, float Time, float SIG, float S0,int sampleStep);
	float GetRPrice(float a, float b, int StepIndex, int npaths, int nsteps, float pS0, float pR, float pSig, float time, float K);
	float Integrand(float x);

	float GetRPrice(float a, float b);
	float GetTPrice(float a, float b);
	float GetSPrice(float a, float b);
	float Get3_8Price(float a, float b);
	void SetS(int amo);
	void Execute();

	float* s_array;
	float* exp_array;

	float GetRPrice(float a, float b, int NumThread/*, float* s_array, float* expf_array*/);
	float GetTPrice(float a, float b, int NumThread);
	float GetSPrice(float a, float b, int NumThread);
	float Get3_8Price(float a, float b, int NumThread);
protected:
	clock_t start, finish;
	double t;
	float EulMarStep(float S, float dt, float dw, float r, float sig);
	float MilsteinStep(float S, float dt, float dw, float r, float sig);
	float RK1Step(float S, float dt, float dw, float r, float sig);
	float BurragePlatenStep(float S, float dt, float dw, float dz, float r, float sig);
	void wAndZProcesses(VSLStreamStatePtr stream, int nsteps, float time, float *buffer);

};

#endif // !____BLACK_SHOLES_NUMERICAL_SOLUTION____
