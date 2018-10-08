#include <iostream>
#include "AnSolution.h"
#include "NumSolution.h"
#include "stdafx.h"
#define TIME		3.0f	// option execute time (years)
#define SIG			0.2f	// volatility; percent per year 0.2 -> 20%
#define R			0.05f	// the interest rate; percent per year 0.05 -> 5%		
#define S0			100.0f	// option price at t == 0
#define NSTEPS		512
#define NPATHS		1000
#define NSAMPLES	1000
#define NUMTHREAD	4
#define K			100.0f // strike price -- price fixed in option

int main(int argc, char* argv[]) {
	//int numVer	 = std::stoi(argv[1]);
	//int NSAMPLES 	 = std::stoi(argv[1]);
	//int numThread	 = std::stoi(argv[2]);
	
	AnSolution as;
	float* pR	= new float[2 * NSTEPS];
	float* pSig = pR + NSTEPS;
	float* pT	= new float[4 * NSAMPLES];
	float* pK	= pT + NSAMPLES;
	float* pS0	= pT + 2 * NSAMPLES;
	float* pC	= pT + 3 * NSAMPLES;
	//float* pS	= pT + 4 * NSAMPLES;

	#if defined(__INTEL_COMPILER) 
	#pragma simd
	#pragma vector always	
	#endif
	for (int i = 0; i < NSTEPS; i++) {
		pR[i]	= R;
		pSig[i] = SIG;
	}
	for (int i = 0; i < NSAMPLES; i++) {
		pT[i]	= TIME;
		pS0[i]	= S0;
		pK[i]	= K;
	}
	float *sbuffer = new float[NSAMPLES];
	std::vector<float> vT (NSAMPLES, TIME);
	std::vector<float> vS0(NSAMPLES, S0);
	std::vector<float> vK (NSAMPLES, K);
	std::vector<float> vC (NSAMPLES, .0f);

	//as.execute(NSTEPS, NPATHS, TIME, S0, R, SIG, NSAMPLES, sbuffer);
	//clock_t start, finish;
	//start = clock();
	//as.baseVer(pT, pK, pS0, pC, NSAMPLES, R, SIG);
	//finish = clock();
	//printf("%lf;\n", (float)(finish - start) / CLOCKS_PER_SEC);

	//as.writeAllToFile(numThread, NSAMPLES, pT, pK, pS0, pC,R, SIG, vT, vK, vS0, vC);

	//start = clock();
	//as._V0(pT, pK, pS0, pC, NSAMPLES, R, SIG);
	//finish = clock();

	//printf("time: %lf;\n", (float)(finish - start) / CLOCKS_PER_SEC);
	//std::cout << "Fair price via BS-formula: " << pC[0] << "\n";
	NumSolution ns;
	//ns.Execute(2, NPATHS, NSTEPS, S0, R, SIG, TIME);
	//std::cout << "Analytical: \n";
	//std::cout << as.simulateStockPrice(NSTEPS, S0, R, SIG, TIME, NSAMPLES, sbuffer) << "\n";
	//std::cout << "Numerical: \n";
	//std::cout << "Euler: "			<< ns.SimulateStockPrices(0, NPATHS, NSTEPS, S0, R, SIG, TIME) << "\n";
	//std::cout << "Milstein: "		<< ns.SimulateStockPrices(1, NPATHS, NSTEPS, S0, R, SIG, TIME) << "\n";
	//std::cout << "RK1: "			<< ns.SimulateStockPrices(2, NPATHS, NSTEPS, S0, R, SIG, TIME) << "\n";
	//std::cout << "Burrage-Platen: " << ns.SimulateStockPrices(3, NPATHS, NSTEPS, S0, R, SIG, TIME) << "\n";
	//std::cout << ns.SimulateStockPricesVol(1, NPATHS, NSTEPS, S0, pR, pSig, TIME) << "\n";

	double t1, t2;
	as.fpVer(pT, pK, pS0, pC, 1, R, SIG);
	std::cout << "BS: \t\t" << pC[0] << "\n";
	t1 = omp_get_wtime();
	float Call = ns.getMCPrice(2, NSTEPS, 0, NSAMPLES, __SEED__, K, R, TIME, SIG, S0);
	t2 = omp_get_wtime();
	std::cout << "Call\t\t" << Call << "\t" << t2 - t1 << "\n";

	//t1 = omp_get_wtime();
	//Call = ns.getMCPricePar(2, NSTEPS, 0, NSAMPLES, __SEED__, K, R, TIME, SIG, S0);
	//t2 = omp_get_wtime();
	//std::cout << "Par Call\t" << Call << "\t" << t2 - t1 << "\n";

	std::cout << "Call via Rectangle: " << ns.GetRPrice(-5.15f, 6.f, 2, NPATHS, NSTEPS, S0, R, SIG, TIME, K) << "\n";

	//ns.getErrors(NSTEPS, 0, 100000, __SEED__, K, R, TIME, SIG, S0, 1000);
	//ns.SimulateStockPrices(3, NPATHS, NSTEPS, S0, R, SIG, TIME) << "\n";
	//std::cout << 
	//delete[] sbuffer;
	delete[] pR;
	delete[] pT;
	//system("pause");
	return 0;
}