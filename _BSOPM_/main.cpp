#include <iostream>
#include "AnSolution.h"
#include "NumSolution.h"
#include "stdafx.h"
#include "../3rdparty/ConfigParser/config.h"
#include "../3rdparty/OptionParser/optionparser.h"

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

	if (config.pBool("USE_VOLATILE_PARAMS")) {
	#define VOLATILE_PARAMS
	#define DELTA_R			config.pFloat("DELTA_R")
	#define DELTA_SIG		config.pFloat("DELTA_SIG")
	}

	//int numVer	 = std::stoi(argv[1]);
	//int NSAMPLES 	 = std::stoi(argv[2]);
	//int numThread	 = std::stoi(argv[3]);
	
	float* pR	= new float[2 * NSTEPS];
	float* pSig = pR + NSTEPS;
	float* pT	= new float[4 * NSAMPLES];
	float* pK	= pT + NSAMPLES;
	float* pS0	= pT + 2 * NSAMPLES;
	float* pC	= pT + 3 * NSAMPLES;

	#if defined(VOLATILE_PARAMS)
	#if defined(__INTEL_COMPILER) 
	#pragma ivdep
	#pragma vector always	
	#endif
	for (int i = 0; i < NSTEPS; i++) {
		pR[i] = R + i * DELTA_R;
		pSig[i] = SIG + i * DELTA_SIG;
	}
	#endif
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
	if (config.pBool("GET_AN_SOLUTION")) {
		as.executeStockPrice(NSTEPS, NPATHS, TIME, S0, R, SIG, sbuffer, __SEED__, INDEXGEN);	// Gives Stock An data in csv
	}
	if (config.pBool("AN_SOLUTION_FAIR")) {
		as.writeOneVersion(VERSION_INDEX, NUMTHREAD, NSAMPLES, pT, pK, pS0, pC, R, SIG, vT, vK, vS0, vC);	// Gives concrete Fair An data in csv
	}
	if (config.pBool("WRITE_ALL_FAIRS")) {
		as.writeAllFairs(NUMTHREAD, NSAMPLES, pT, pK, pS0, pC, R, SIG, vT, vK, vS0, vC);			//  Gives All Fair An data in csv
	}

	NumSolution ns;
	//
	//ns.writeMethodConvergence(2, INDEXGEN, NPATHS, NSTEPS, S0, R, SIG, TIME, __SEED__); // Burrage-Platen is incorrect
	std::cout << "Analytical: \t\t";
	//std::cout << as.simulateStockPriceAn(NPATHS, S0, R, SIG, TIME, sbuffer, __SEED__, INDEXGEN) << "\n";
	//std::cout << "Numerical: \n";
	//std::cout << "Euler: \t\t\t"			<< ns.SimulateStockPrices(0, INDEXGEN, NPATHS, NSTEPS, S0, R, SIG, TIME, __SEED__) << "\n";
	//std::cout << "Milstein: \t\t"			<< ns.SimulateStockPrices(1, INDEXGEN, NPATHS, NSTEPS, S0, R, SIG, TIME, __SEED__) << "\n";
	//std::cout << "RK1: \t\t\t"				<< ns.SimulateStockPrices(2, INDEXGEN, NPATHS, NSTEPS, S0, R, SIG, TIME, __SEED__) << "\n";
	//std::cout << "Burrage-Platen: \t"		<< ns.SimulateStockPrices(3, INDEXGEN, NPATHS, NSTEPS, S0, R, SIG, TIME, __SEED__) << "\n";
	//std::cout << "Milstein Vol: \t\t"		<< ns.SimulateStockPricesVol(1, NPATHS, NSTEPS, S0, pR, pSig, TIME, __SEED__, INDEXGEN) << "\n";
	//std::cout << "MC Fair Price via RK1: \t"<< ns.getMCPrice(2, NSTEPS, INDEXGEN, NSAMPLES, __SEED__, K, R, TIME, SIG, S0) << "\n";

	//as.fpVer(pT, pK, pS0, pC, 1, R, SIG);
	//std::cout << "BS: \t\t\t" << pC[0] << "\n";

	//double t1, t2, t = 0.0;

	//t1 = omp_get_wtime();
	//ns.GetRPricePar(-5.15f, 5.5f, 2000, 1, NSAMPLES, R, SIG, pT, pK, pS0, pC);
	//t2 = omp_get_wtime();
	//std::cout << "Par Call\t" << pC[0] << "\t" << t2 - t1 << "\n";

	//t1 = omp_get_wtime();
	//ns.GetTPricePar(-5.15f, 5.5f, 2000, 1, NSAMPLES, R, SIG, pT, pK, pS0, pC);
	//t2 = omp_get_wtime();
	//std::cout << "Par Call\t" << pC[0] << "\t" << t2 - t1 << "\n";
	//
	//t1 = omp_get_wtime();
	//ns.GetSPricePar(-5.15f, 5.5f, 2000, 1, NSAMPLES, R, SIG, pT, pK, pS0, pC);
	//t2 = omp_get_wtime();
	//std::cout << "Par Call\t" << pC[0] << "\t" << t2 - t1 << "\n";

	//t1 = omp_get_wtime();
	//ns.Get3_8PricePar(-5.15f, 5.5f, 2000, 1, NSAMPLES, R, SIG, pT, pK, pS0, pC);
	//t2 = omp_get_wtime();
	//std::cout << "Par Call\t" << pC[0] << "\t" << t2 - t1 << "\n";



	//t1 = omp_get_wtime();
	//ns.Get3_8PricePar(-5.15f, 6.f, 2000, NUMTHREAD, NSAMPLES, R, SIG, pT, pK, pS0, pC);
	//t2 = omp_get_wtime();
	//std::cout << "MC Simp\t\t" << pC[0] << "\t" << t2 - t1 << "\n";


	//t1 = omp_get_wtime();
	//float Call = ns.getMCPrice(3, NSTEPS, INDEXGEN, NSAMPLES, __SEED__, K, R, TIME, SIG, S0);
	//t2 = omp_get_wtime();
	//std::cout << "Call\t\t" << Call << "\t" << t2 - t1 << "\n";

	//t1 = omp_get_wtime();
	//Call = ns.getMCPricePar(NUMTHREAD, 0, NSTEPS, 0, NSAMPLES, __SEED__, K, R, TIME, SIG, S0, t);
	//t2 = omp_get_wtime();
	//std::cout << "Par Call\t" << Call << "\t" << t2 - t1 << "\n";

	//ns.MCParExecute(3, NSTEPS, INDEXGEN, NSAMPLES, __SEED__, K, R, TIME, SIG, S0);

	//std::cout << "Call via Rectangle: " << ns.GetRPrice(-5.15f, 6.f, 2, NPATHS, NSTEPS, S0, R, SIG, TIME, K) << "\n";

	//ns.getErrors(NSTEPS, 0, 100000, __SEED__, K, R, TIME, SIG, S0, 1000);
	delete[] sbuffer;
	delete[] pR;
	delete[] pT;
	//system("pause");
	return 0;
}