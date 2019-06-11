#include <iostream>
#include "AnSolution.h"
#include "NumSolution.h"
#include "config.h"
#include "MDim.h"
#include <vector>
#include <algorithm>

int main(int argc, char* argv[], char* envp[]) {

	Config config("BSOPM.config", envp);

	if (config.pBool("USE_BSOPM_DEFINES")) {
#define TIME		config.pFloat("TIME")	// option execute time (years)
#define SIG			config.pFloat("SIG")	// volatility; percent per year 0.2 -> 20%
#define R			config.pFloat("R")		// the interest rate; percent per year 0.05 -> 5%
#define S0			config.pFloat("S0")		// option price at t == 0
#define NSTEPS		config.pInt("NSTEPS")
#define NPATHS		config.pInt("NPATHS")
#define NSAMPLES	config.pInt("NSAMPLES")
#define NUMTHREAD	config.pInt("NUMTHREAD")
#define INDEXGEN	config.pInt("INDEXGEN")
#define K			config.pFloat("K")		// strike price -- price fixed in option
#define __SEED__	config.pInt("__SEED__")
#define STEP_INDEX  config.pInt("STEP_INDEX")
#define VERSION_INDEX config.pInt("VERSION_INDEX")
	}

	float* pR = new float[2 * NSTEPS];

	if (config.pBool("USE_VOLATILE_PARAMS")) {

#define DELTA_R			config.pFloat("DELTA_R")
#define DELTA_SIG		config.pFloat("DELTA_SIG")


		float* pSig = pR + NSTEPS;

#if defined(__INTEL_COMPILER) 
#pragma ivdep
#pragma vector always	
#endif
		for (int i = 0; i < NSTEPS; i++) {
			pR[i] = R + i * DELTA_R;
			pSig[i] = SIG + i * DELTA_SIG;
		}
	}

	float* pT = new float[4 * NSAMPLES];
	float* pK = pT + NSAMPLES;
	float* pS0 = pT + 2 * NSAMPLES;
	float* pC = pT + 3 * NSAMPLES;

	for (int i = 0; i < NSAMPLES; i++) {
		pT[i] = TIME;
		pS0[i] = S0;
		pK[i] = K;
	}
	float *sbuffer = new float[NPATHS];
	std::vector<float> vT(NSAMPLES, TIME);
	std::vector<float> vS0(NSAMPLES, S0);
	std::vector<float> vK(NSAMPLES, K);
	std::vector<float> vC(NSAMPLES, .0f);

	double workTime;
	bool is_optimized;

	config.pBool("IS_OPTIMIZED") ? is_optimized = true : is_optimized = false;

	if (config.pBool("GET_AN_SOLUTION")) {
		AnSolution as;
		as.executeStockPrice(NSTEPS, NPATHS, TIME, S0, R, SIG, sbuffer, __SEED__, INDEXGEN);	// Gives Stock An data in csv
	}
	if (config.pBool("AN_SOLUTION_FAIR")) {
		AnSolution as;

		//as.writeOneVersion(VERSION_INDEX, NUMTHREAD, NSAMPLES, pT, pK, pS0, pC, R, SIG, vT, vK, vS0, vC);	// Gives concrete Fair An data in csv
	}
	if (config.pBool("WRITE_ALL_FAIR")) {
		AnSolution as;
		as.writeAllFairs(NUMTHREAD, NSAMPLES, pT, pK, pS0, pC, R, SIG, vT, vK, vS0, vC);			//  Gives All Fair An data in csv
	}
	if (config.pBool("AN_SOLUTION_PAR")) {
		AnSolution as;
		as.simulateStockPriceAnOptedPar(NPATHS, S0, R, SIG, TIME, __SEED__, INDEXGEN, is_optimized, NUMTHREAD, workTime);
		std::cout << workTime << "\n";
	}

	NumSolution ns;
	AnSolution as;
	double t1, t2;
	t1 = omp_get_wtime();
	as.baseVer(pT, pK, pS0, pC, 1, R, SIG);
	t2 = omp_get_wtime();
	std::cout << "BS time " << t2 - t1 << std::endl;


	delete[] sbuffer;
	delete[] pR;
	delete[] pT;
	//system("pause");
	#define DIM	config.pInt("DIM")
	float* SP = new float[DIM];
	std::cout.precision(7);
	float eps = 1e-3f;

	float I_1_n, I_1_2n;
	//for (int ind = 0; ind < 4; ++ind) {
	//	int n0 = 50;
	//	float delta = 1000.f;
	//	while (delta > eps) {
	//		I_1_n = ns.SimulateStockPrices(ind, INDEXGEN, n0, NSTEPS, S0, R, SIG, TIME, __SEED__);
	//		I_1_2n = ns.SimulateStockPrices(ind, INDEXGEN, 2 * n0, NSTEPS, S0, R, SIG, TIME, __SEED__);
	//		delta = abs(1 / 3.f * (I_1_n - I_1_2n));
	//		std::cout << n0 << "\t" << delta << "\t" << I_1_n << "\t" << I_1_2n << std::endl;
	//		n0 *= 2;
	//	}
	//}
	//std::cout << getRRNG(DIM, n0, -15.f, 6.f, 0.05f, 3.f, SP, 100.f, 100.f) << std::endl;

	//int n_steps = 100;
	//for (int ind = 0; ind < 4; ++ind) {
	//	
	//}

	//std::vector<double> times;
	//for (int ind = 0; ind < 3; ++ind) {
	//	for (int i = 0; i < 11; ++i) {
	//		times.push_back(ns.SimulateStockPrices(ind, INDEXGEN, 150000, NSTEPS, S0, R, SIG, TIME, __SEED__));
	//	}
	//	std::sort(times.begin(), times.end());
	//	std::cout << "Avg time is " << times[5] << std::endl;
	//	times.clear();
	//}
	//float SPAn = 116.4776;
	// euler, RK1, mil -- 150000 0.005
	// BP -- 80000
	// 116.177
	t1 = omp_get_wtime();
	ns.getMCPrice(0, 512, 0, 1, __SEED__, K, R, TIME, SIG, S0);
	t2 = omp_get_wtime();
	std::cout << "MC time " << t2 - t1 << std::endl;


	
	//Npaths for BP == 290k, enother -- 640k
	//as.simulateStockPriceAn(NSAMPLES, S0, R, SIG, TIME, __SEED__, INDEXGEN, workTime);
	//for (int i = 0; i < 100; ++i) {
	//	//std::cout /*<< i * 1000 << "\t" */<< as.simulateStockPriceAn(i * 1000, 100.0, 0.05, 0.2, 3.0, __SEED__, 0, workTime) << std::endl;
	//	//std::cout << 100 + i * 10 << "\n";
	//	std::cout /*<< 100 + i * 10 << "\t"*/ << ns.SimulateStockPrices(3, 0, 10000 * i, 290, 100.f, 0.05f, 0.2f, 3.f, __SEED__) - 116.177 << std::endl;
	//}

	//NSTEPS for BP == 290
	//std::cout << "\t" << as.simulateStockPriceAn(300000, S0, R, SIG, TIME, __SEED__, 0, workTime);
	//ns.SimulateStockPrices(0, INDEXGEN, 1000000, 520, S0, R, SIG, TIME, __SEED__);
	//ns.SimulateStockPrices(1, INDEXGEN, 1000000, 520, S0, R, SIG, TIME, __SEED__);
	//ns.SimulateStockPrices(2, INDEXGEN, 1000000, 520, S0, R, SIG, TIME, __SEED__);
	//ns.SimulateStockPrices(3, INDEXGEN, 1000000, 290, S0, R, SIG, TIME, __SEED__);
	/*ns.SimulateStockPrices(3, INDEXGEN, 80000, NSTEPS, S0, R, SIG, TIME, __SEED__);*/
	//getMCPriceParNew(NUMTHREAD, 0, NSTEPS, INDEXGEN, NSAMPLES, 5000, __SEED__, K, R, TIME, SIG, S0, DIM);

	//ns.MCParExecute(STEP_INDEX, NSTEPS, INDEXGEN, NSAMPLES, __SEED__, K, R, TIME, SIG, S0);
	//ns.getMCPrice(STEP_INDEX, NSTEPS, INDEXGEN, NSAMPLES, __SEED__, K, R, TIME, SIG, S0);
	delete[] SP;
	system("pause");
	return 0;
}
