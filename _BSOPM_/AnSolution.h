// Copyright 2018 Romanov Alexander

#include "BSOPM.h"
#ifndef ____BLACK_SHOLES_ANALYTICAL_SOLUTION____
#define ____BLACK_SHOLES_ANALYTICAL_SOLUTION____

class AnSolution : public BSOPM {
public:
	AnSolution() {};
	~AnSolution() {};

	float simulateStockPrice(int npaths, float s0, float r, float sig, float time, int nsamples, float* sbuffer);

	void baseVer		(float* pT, float* pK, float* pS0, float* pC, int nsamples, float r, float sig);
	void eqBaseVer		(float* pT, float* pK, float* pS0, float* pC, int nsamples, float r, float sig);
	void updEqBaseVer	(float* pT, float* pK, float* pS0, float* pC, int nsamples, float r, float sig);
	void fpVer			(float* pT, float* pK, float* pS0, float* pC, int nsamples, float r, float sig);
	void erfVer			(float* pT, float* pK, float* pS0, float* pC, int nsamples, float r, float sig);
	void fpVerOld		(float* pT, float* pK, float* pS0, float* pC, int nsamples, float r, float sig);
	void erfVerOld		(float* pT, float* pK, float* pS0, float* pC, int nsamples, float r, float sig);
	void vectVer		(float* pT, float* pK, float* pS0, float* pC, int nsamples, float r, float sig);
	void vectErfVer		(float* pT, float* pK, float* pS0, float* pC, int nsamples, float r, float sig);
	void vectInvSqrtVer	(float* pT, float* pK, float* pS0, float* pC, int nsamples, float r, float sig);
	void parVer			(float* pT, float* pK, float* pS0, float* pC, int nsamples, float r, float sig);
	void parVerOld		(float* pT, float* pK, float* pS0, float* pC, int nsamples, float r, float sig);
	void parNontempVer	(float* pT, float* pK, float* pS0, float* pC, int nsamples, float r, float sig);
	void restrictVer	(float* __restrict pT, float* __restrict pK, float* __restrict pS0, float* __restrict pC, int nsamples, float r, float sig);
	void stdVectorVer	(std::vector<float>& T, std::vector<float>& K, std::vector<float>& S0, std::vector<float>& C, int nsamples, float r, float sig);

	void writeToFile	(float* buffer, int nrows, float avg);
	void writeAllToFile(int num_Threads, int N, float* pT, float* pK, float* pS0, float* pC, float r, float sig, std::vector<float>& T, std::vector<float>& K, std::vector<float>& S0, std::vector<float>& C);
	float writeToFile	(int numVer, int num_Threads, int N, float* pT, float* pK, float* pS0, float* pC, float r, float sig, std::vector<float>& T, std::vector<float>& K, std::vector<float>& S0, std::vector<float>& C);
	void execute		(int nsteps, int npaths, float time, float s0, float r, float sig, int nsamples, float* sbuffer);


};

#endif // !____BLACK_SHOLES_ANALYTICAL_SOLUTION____
