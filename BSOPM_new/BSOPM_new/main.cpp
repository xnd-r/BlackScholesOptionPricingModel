#include <iostream>
#include "AnSolution.h"
#include "NumSolution.h"
#include "config.h"
#include "MDim.h"

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
		as.writeOneVersion(VERSION_INDEX, NUMTHREAD, NSAMPLES, pT, pK, pS0, pC, R, SIG, vT, vK, vS0, vC);	// Gives concrete Fair An data in csv
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
	delete[] sbuffer;
	delete[] pR;
	delete[] pT;
	//system("pause");
	#define DIM	config.pInt("DIM")
	float* SP = new float[DIM];
	std::cout.precision(7);
	float eps = 1e-3f;
	int n0 = 50;
	float delta = 1000.f;
	float I_1_n, I_1_2n;
	//while (delta > eps) {
	//	I_1_n = getRRNG(DIM, n0, -15.f, 6.f, 0.05f, 3.f, SP, 100.f, 100.f);
	//	I_1_2n = getRRNG(DIM, 2 * n0, -15.f, 6.f, 0.05f, 3.f, SP, 100.f, 100.f);
	//	delta = abs(1 / 3.f * (I_1_n - I_1_2n));
	//	std::cout << n0 << "\t" << delta << "\t" << I_1_n << "\t" << I_1_2n << std::endl;
	//	n0 *= 2;
	//}
	//std::cout << getRRNG(DIM, n0, -15.f, 6.f, 0.05f, 3.f, SP, 100.f, 100.f) << std::endl;

	
	n0 = 50; 
	delta = 1000.f;
	while (delta > eps) {
		I_1_n = getMCPriceParNew(NUMTHREAD, 0, NSTEPS, INDEXGEN, n0, 5000, __SEED__, K, R, TIME, SIG, S0, DIM);
		I_1_2n = getMCPriceParNew(NUMTHREAD, 0, NSTEPS, INDEXGEN, n0 * 2, 5000, __SEED__, K, R, TIME, SIG, S0, DIM);
		delta = abs(1 / 3.f * (I_1_n - I_1_2n));
		std::cout << n0 << "\t" << delta << "\t" << I_1_n << "\t" << I_1_2n << std::endl;
		n0 *= 2;
	}

	delete[] SP;
	system("pause");
	return 0;
}
