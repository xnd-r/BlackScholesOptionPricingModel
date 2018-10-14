#include <iostream>
#include "AnSolution.h"
#include "NumSolution.h"
#include "stdafx.h"
#define TIME		1.f	// option execute time (years)
#define SIG			0.2f	// volatility; percent per year 0.2 -> 20%
#define R			0.05f	// the interest rate; percent per year 0.05 -> 5%		
#define S0			100.0f	// option price at t == 0
#define NSTEPS		1000
#define NPATHS		2000
#define NSAMPLES	15000
#define NUMTHREAD	4
#define INDEXGEN	0
#define K			100.0f // strike price -- price fixed in option


void memoryInit(){
}
int main(int argc, char* argv[]) {
	//int numVer	 = std::stoi(argv[1]);
	//int NSAMPLES 	 = std::stoi(argv[2]);
	//int numThread	 = std::stoi(argv[3]);
	
	float* pR	= new float[2 * NSTEPS];
	float* pSig = pR + NSTEPS;
	float* pT	= new float[4 * NSAMPLES];
	float* pK	= pT + NSAMPLES;
	float* pS0	= pT + 2 * NSAMPLES;
	float* pC	= pT + 3 * NSAMPLES;

	#if defined(__INTEL_COMPILER) 
	#pragma ivdep
	#pragma vector always	
	#endif
	for (int i = 0; i < NSTEPS; i++) {
		pR[i] = R + i * 0.00001;
		pSig[i] = SIG + i * 0.00001;
	}
	for (int i = 0; i < NSAMPLES; i++) {
		pT[i]	= TIME;
		pS0[i]	= S0;
		pK[i]	= K;
	}
	float *sbuffer = new float[NPATHS];
	std::vector<float> vT (NSAMPLES, TIME);
	std::vector<float> vS0(NSAMPLES, S0);
	std::vector<float> vK (NSAMPLES, K);
	std::vector<float> vC (NSAMPLES, .0f);

	AnSolution as;
	//as.executeStockPrice(NSTEPS, NPATHS, TIME, S0, R, SIG, sbuffer, __SEED__, INDEXGEN);
	//as.writeAllFairs(NUMTHREAD, NSAMPLES, pT, pK, pS0, pC,R, SIG, vT, vK, vS0, vC);
	//as.writeOneVersion(2, NUMTHREAD, NSAMPLES, pT, pK, pS0, pC, R, SIG, vT, vK, vS0, vC);

	NumSolution ns;
	//
	ns.Execute(1, INDEXGEN, NPATHS, NSTEPS, S0, R, SIG, TIME, __SEED__);
	//std::cout << "Analytical: \t\t";
	//std::cout << as.simulateStockPriceAn(NPATHS, S0, R, SIG, TIME, sbuffer, __SEED__, INDEXGEN) << "\n";
	//std::cout << "Numerical: \n";
	//std::cout << "Euler: \t\t\t"			<< ns.SimulateStockPrices(0, INDEXGEN, NPATHS, NSTEPS, S0, R, SIG, TIME, __SEED__) << "\n";
	//std::cout << "Milstein: \t\t"			<< ns.SimulateStockPrices(1, INDEXGEN, NPATHS, NSTEPS, S0, R, SIG, TIME, __SEED__) << "\n";
	//std::cout << "RK1: \t\t\t"				<< ns.SimulateStockPrices(2, INDEXGEN, NPATHS, NSTEPS, S0, R, SIG, TIME, __SEED__) << "\n";
	//std::cout << "Burrage-Platen: \t"		<< ns.SimulateStockPrices(3, INDEXGEN, NPATHS, NSTEPS, S0, R, SIG, TIME, __SEED__) << "\n";
	//std::cout << "Milstein Vol: \t\t" << ns.SimulateStockPricesVol(1, NPATHS, NSTEPS, S0, pR, pSig, TIME, __SEED__, INDEXGEN) << "\n";
	//std::cout << "MC Price via RK1: \t" << ns.getMCPrice(2, NSTEPS, INDEXGEN, NSAMPLES, __SEED__, K, R, TIME, SIG, S0) << "\n";



	//double t1, t2, t = 0.0;
	//as.fpVer(pT, pK, pS0, pC, 1, R, SIG);
	//std::cout << "BS: \t\t" << pC[0] << "\n";
	//t1 = omp_get_wtime();
	//float Call = ns.getMCPrice(0, NSTEPS, 0, NSAMPLES, __SEED__, K, R, TIME, SIG, S0);
	//t2 = omp_get_wtime();
	//std::cout << "Call\t\t" << Call << "\t" << t2 - t1 << "\n";

	//t1 = omp_get_wtime();
	//Call = ns.getMCPricePar(NUMTHREAD, 0, NSTEPS, 0, NSAMPLES, __SEED__, K, R, TIME, SIG, S0, t);
	//t2 = omp_get_wtime();
	//std::cout << "Par Call\t" << Call << "\t" << t2 - t1 << "\n";

	//ns.MCExecute(2, NSTEPS, INDEXGEN, NSAMPLES, __SEED__, K, R, TIME, SIG, S0);

	//std::cout << "Call via Rectangle: " << ns.GetRPrice(-5.15f, 6.f, 2, NPATHS, NSTEPS, S0, R, SIG, TIME, K) << "\n";

	//ns.getErrors(NSTEPS, 0, 100000, __SEED__, K, R, TIME, SIG, S0, 1000);
	delete[] sbuffer;
	delete[] pR;
	delete[] pT;
	//system("pause");
	return 0;
}