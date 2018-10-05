#include "AnSolution.h"
#include <iostream>

float AnSolution::simulateStockPrice(
	int npaths, float s0, float r, float sig, float time, int nsamples, float* sbuffer) {

	VSLStreamStatePtr stream = initGen();
	float* dw = new float[npaths]; // random values with N(0, time) buffer
	float stockPrice = .0f;

	normalGenerator(0.f, sqrtf(time), npaths, stream, dw);
#if defined(__INTEL_COMPILER) 
		#pragma simd
		#pragma vector always	
	#endif
	for (int i = 0; i < npaths; i++) {
		sbuffer[i] = getStockPrice(s0, r, sig, dw[i], time);
		stockPrice += sbuffer[i];
	}

	freeGen(stream);
	delete[] dw;
	return stockPrice / npaths;
}

void AnSolution::writeToFile(float* buffer, int nrows, float avg) {

	time_t timestamp;
	time(&timestamp);
	std::string date = asctime(localtime(&timestamp));
	date.pop_back();
	date.append("_GBM_Analitycal_Simple.csv");

	for (std::string::iterator it = date.begin(); it<date.end(); ++it) {
		if (*it == ':') {
			date.erase(it);
		}
		std::replace(date.begin(), date.end(), ' ', '_');
	}

	FILE *f = fopen(date.c_str(), "w");
	for (int j = 0; j < nrows; j++) {
		std::string tmp_cell = std::to_string(buffer[j]);

		for (std::string::iterator it = tmp_cell.begin(); it<tmp_cell.end(); ++it) {
			std::replace(tmp_cell.begin(), tmp_cell.end(), '.', ',');
		}
		fprintf(f, /*"%lf;\n"*/"%s;", /*buffer[j]*/tmp_cell.c_str());
	}
	fprintf(f, "\n");
	fprintf(f, "Average price = %lf;\n", avg);
	fclose(f);
}

void AnSolution::execute(int nsteps, int npaths, float time, float pS0, float pR, float pSig, int nsamples, float* sbuffer ) {

	float avg = simulateStockPrice(npaths, pS0, pR, pSig, time, nsamples, sbuffer);
	writeToFile(sbuffer, npaths, avg);
	printf("Average price = %lf\n", avg);
}

// preliminary
void AnSolution::baseVer(float* pT, float* pK, float* pS0, float* pC, int nsamples, float r, float sig)
{
	float d1, d2, p1, p2;
	for (int i = 0; i < nsamples; i++)
	{
		d1 = (log(pS0[i] / pK[i]) + (r + sig * sig * 0.5) *
			pT[i]) / (sig * sqrt(pT[i]));
		d2 = (log(pS0[i] / pK[i]) + (r - sig * sig * 0.5) *
			pT[i]) / (sig * sqrt(pT[i]));
		vsCdfNorm(1, &d1, &p1);
		vsCdfNorm(1, &d2, &p2);
		pC[i] = pS0[i] * p1 - pK[i] * exp((-1.0) * r * pT[i]) * p2;
	}
}

// equivalent transformation of d[1]

void AnSolution::eqBaseVer(float* pT, float* pK, float* pS0, float* pC, int nsamples, float r, float sig)
{
	float d1, d2, p1, p2;
	for (int i = 0; i < nsamples; i++)
	{
		d1 = (log(pS0[i] / pK[i]) + (r + sig * sig * 0.5) *
			pT[i]) / (sig * sqrt(pT[i]));
		d2 = d1 - sig * sqrt(pT[i]);
		vsCdfNorm(1, &d1, &p1);
		vsCdfNorm(1, &d2, &p2);
		pC[i] = pS0[i] * p1 - pK[i] * exp((-1.0) * r * pT[i]) * p2;
	}
}

// upd equivalent transformation of d[1]
void AnSolution::updEqBaseVer(float* pT, float* pK, float* pS0, float* pC, int nsamples, float r, float sig) {
	float d[2], p[2]; 
	for (int i = 0; i < nsamples; i++) {
		d[0] = (log(pS0[i] / pK[i]) + (r + sig * sig * 0.5) * pT[i]) / (sig * sqrt(pT[i]));
		d[1] = d[0] - sig * sqrt(pT[i]);
		vsCdfNorm(2, d, p);
		pC[i] = pS0[i] * p[0] - pK[i] * exp((-1.0) * r * pT[i]) * p[1];
	}
}

// preliminary via std::vector
void AnSolution::stdVectorVer(std::vector<float>& T, std::vector<float>& K, std::vector<float>& S0, std::vector<float>& C, int nsamples, float r, float sig)
{
	float d[2], p[2];
	for (int i = 0; i < nsamples; i++)
	{
		d[0] = (logf(S0[i] / K[i]) + (r + sig * sig * 0.5f) *
			T[i]) / (sig * sqrtf(T[i]));
		d[1] = d[0] - sig * sqrtf(T[i]);

		vsCdfNorm(2, d, p);

		C[i] = S0[i] * p[0] - K[i] * expf((-1.0f) * r * T[i]) * p[1];
	}

}

// fp
void AnSolution::fpVer(float* pT, float* pK, float* pS0, float* pC, int nsamples, float r, float sig)
{
	float d[2], p[2];
	for (int i = 0; i < nsamples; i++)
	{
		d[0] = (logf(pS0[i] / pK[i]) + (r + sig * sig * 0.5f) *
			pT[i]) / (sig * sqrtf(pT[i]));
		d[1] = d[0] - sig * sqrtf(pT[i]);

		vsCdfNorm(2, d, p);

		pC[i] = pS0[i] * p[0] - pK[i] * expf((-1.0f) * r * pT[i]) * p[1];
	}
}

// fpOld
void AnSolution::fpVerOld(float* pT, float* pK, float* pS0, float* pC, int nsamples, float r, float sig)
{
	float d1, d2, p1, p2;
	for (int i = 0; i < nsamples; i++)
	{
		d1 = (logf(pS0[i] / pK[i]) + (r + sig * sig * 0.5f) *
			pT[i]) / (sig * sqrtf(pT[i]));
		d2 = (logf(pS0[i] / pK[i]) + (r - sig * sig * 0.5f) *
			pT[i]) / (sig * sqrtf(pT[i]));

		vsCdfNorm(1, &d1, &p1);
		vsCdfNorm(1, &d2, &p2);

		pC[i] = pS0[i] * p1 - pK[i] * expf((-1.0f) * r * pT[i]) * p2;
	}
}

// erf
void AnSolution::erfVer(float* pT, float* pK, float* pS0, float* pC, int nsamples, float r, float sig)
{
	float d[2], cdf1, cdf2;
	for (int i = 0; i < nsamples; i++)
	{
		d[0] = (logf(pS0[i] / pK[i]) + (r + sig * sig * 0.5f) * pT[i]) / (sig * sqrtf(pT[i]));
		d[1] = d[0] - sig * sqrtf(pT[i]);

		cdf1 = 0.5f + 0.5f * erff(d[0] / sqrtf(2.0f));
		cdf2 = 0.5f + 0.5f * erff(d[1] / sqrtf(2.0f));

		pC[i] = pS0[i] * cdf1 - pK[i] * expf((-1.0f) * r * pT[i]) * cdf2;
	}
}

// erfVerOld
void AnSolution::erfVerOld(float* pT, float* pK, float* pS0, float* pC, int nsamples, float r, float sig)
{
	float d1, d2, erf1, erf2;
	for (int i = 0; i < nsamples; i++)
	{
		d1 = (logf(pS0[i] / pK[i]) + (r + sig * sig * 0.5f) *
			pT[i]) / (sig * sqrtf(pT[i]));
		d2 = (logf(pS0[i] / pK[i]) + (r - sig * sig * 0.5f) *
			pT[i]) / (sig * sqrtf(pT[i]));

		erf1 = 0.5f + 0.5f * erff(d1 / sqrtf(2.0f));
		erf2 = 0.5f + 0.5f * erff(d2 / sqrtf(2.0f));

		pC[i] = pS0[i] * erf1 - pK[i] * expf((-1.0f) * r * pT[i]) * erf2;
	}
}

// restrict
void AnSolution::restrictVer(float* __restrict pT, float* __restrict pK, float* __restrict pS0, float* __restrict pC, int nsamples, float r, float sig)
{
	float d[2], cdf1, cdf2;
	for (int i = 0; i < nsamples; i++)
	{
		d[0] = (logf(pS0[i] / pK[i]) + (r + sig * sig * 0.5f) *
			pT[i]) / (sig * sqrtf(pT[i]));
		d[1] = (logf(pS0[i] / pK[i]) + (r - sig * sig * 0.5f) *
			pT[i]) / (sig * sqrtf(pT[i]));

		cdf1 = 0.5f + 0.5f * erff(d[0] / sqrtf(2.0f));
		cdf2 = 0.5f + 0.5f * erff(d[1] / sqrtf(2.0f));

		pC[i] = pS0[i] * cdf1 - pK[i] * expf((-1.0f) * r * pT[i]) * cdf2;
	}
}

// fp_simd_vector_always
void AnSolution::vectVer(float* pT, float* pK, float* pS0, float* pC, int nsamples, float r, float sig)
{
	float d[2], p[2];
#if defined (__INTEL_COMPILER)
#pragma simd
#pragma ivdep
#pragma vector always	
#endif
	for (int i = 0; i < nsamples; i++)
	{
		d[0] = (logf(pS0[i] / pK[i]) + (r + sig * sig * 0.5f) *
			pT[i]) / (sig * sqrtf(pT[i]));
		d[1] = (logf(pS0[i] / pK[i]) + (r - sig * sig * 0.5f) *
			pT[i]) / (sig * sqrtf(pT[i]));

		vsCdfNorm(2, d, p);
		pC[i] = pS0[i] * p[0] - pK[i] * expf((-1.0f) * r * pT[i]) * p[1];
	}
}

// fp_erf_simd_vector_always
void AnSolution::vectErfVer(float* pT, float* pK, float* pS0, float* pC, int nsamples, float r, float sig)
{
	float d[2], cdf1, cdf2;
#if defined (__INTEL_COMPILER)
#pragma simd
#pragma vector always
#endif
	for (int i = 0; i < nsamples; i++)
	{
		d[0] = (logf(pS0[i] / pK[i]) + (r + sig * sig * 0.5f) *
			pT[i]) / (sig * sqrtf(pT[i]));
		d[1] = (logf(pS0[i] / pK[i]) + (r - sig * sig * 0.5f) *
			pT[i]) / (sig * sqrtf(pT[i]));
		cdf1 = 0.5f + 0.5f * erff(d[0] / sqrtf(2.0f));
		cdf2 = 0.5f + 0.5f * erff(d[1] / sqrtf(2.0f));

		pC[i] = pS0[i] * cdf1 - pK[i] * expf((-1.0f) * r * pT[i]) * cdf2;
	}
}

// using_defined_INVSQRT2
void AnSolution::vectInvSqrtVer(float* pT, float* pK, float* pS0, float* pC, int nsamples, float r, float sig)
{
	float d[2], cdf1, cdf2;
	float invsqrt2 = 0.707106781f;
#if defined (__INTEL_COMPILER)
#pragma simd
#endif
	for (int i = 0; i < nsamples; i++)
	{
		d[0] = (logf(pS0[i] / pK[i]) + (r + sig * sig * 0.5f) * pT[i]) / (sig * sqrtf(pT[i]));
		d[1] = (logf(pS0[i] / pK[i]) + (r - sig * sig * 0.5f) * pT[i]) / (sig * sqrtf(pT[i]));

		cdf1 = 0.5f + 0.5f * erff(d[0] * invsqrt2);
		cdf2 = 0.5f + 0.5f * erff(d[1] * invsqrt2);

		pC[i] = pS0[i] * cdf1 - pK[i] * expf((-1.0f) * r * pT[i]) * cdf2;
	}
}

//omp
void AnSolution::parVer(float* pT, float* pK, float* pS0, float* pC, int nsamples, float r, float sig)
{
	float d[2], cdf1, cdf2;
#if defined (__INTEL_COMPILER)
#pragma simd
#endif
	#pragma omp parallel for private(d, cdf1, cdf2)
		for (int i = 0; i < nsamples; i++)
		{
			d[0] = (logf(pS0[i] / pK[i]) + (r + sig * sig * 0.5f) *
				pT[i]) / (sig * sqrtf(pT[i]));
			d[1] = (logf(pS0[i] / pK[i]) + (r - sig * sig * 0.5f) *
				pT[i]) / (sig * sqrtf(pT[i]));
	
			cdf1 = 0.5f + 0.5f * erff(d[0] / sqrtf(2.0f));
			cdf2 = 0.5f + 0.5f * erff(d[1] / sqrtf(2.0f));
	
			pC[i] = pS0[i] * cdf1 - pK[i] * expf((-1.0f) * r * pT[i]) * cdf2;
		}
}

void AnSolution::parVerOld(float* pT, float* pK, float* pS0, float* pC, int N, float r, float sig)
{
	int		i;
	float	d1, d2, erf1, erf2;

#pragma simd
#pragma omp parallel for private(d1, d2, erf1, erf2)
	for (i = 0; i < N; i++)
	{
		d1 = (logf(pS0[i] / pK[i]) + (r + sig * sig * 0.5f) *
			pT[i]) / (sig * sqrtf(pT[i]));
		d2 = (logf(pS0[i] / pK[i]) + (r - sig * sig * 0.5f) *
			pT[i]) / (sig * sqrtf(pT[i]));

		erf1 = 0.5f + 0.5f * erff(d1 / sqrtf(2.0f));
		erf2 = 0.5f + 0.5f * erff(d2 / sqrtf(2.0f));

		pC[i] = pS0[i] * erf1 - pK[i] * expf((-1.0f) * r * pT[i]) * erf2;
	}
}

//omp_nontemporal
void AnSolution::parNontempVer(float* pT, float* pK, float* pS0, float* pC, int nsamples, float r, float sig)
{
	float d[2], cdf1, cdf2;
#if defined (__INTEL_COMPILER)
#pragma simd 
#pragma vector nontemporal
#endif
	#pragma omp parallel for private(d, cdf1, cdf2)
		for (int i = 0; i < nsamples; i++)
		{
			d[0] = (logf(pS0[i] / pK[i]) + (r + sig * sig * 0.5f) *
				pT[i]) / (sig * sqrtf(pT[i]));
			d[1] = (logf(pS0[i] / pK[i]) + (r - sig * sig * 0.5f) *
				pT[i]) / (sig * sqrtf(pT[i]));
	
			cdf1 = 0.5f + 0.5f * erff(d[0] / sqrtf(2.0f));
			cdf2 = 0.5f + 0.5f * erff(d[1] / sqrtf(2.0f));
	
			pC[i] = pS0[i] * cdf1 - pK[i] * expf((-1.0f) * r * pT[i]) * cdf2;
		}
}

void AnSolution::writeAllToFile(int num_Threads, int N, float* pT, float* pK, float* pS0, float* pC, float r, float sig, std::vector<float>& T, std::vector<float>& K, std::vector<float>& S0, std::vector<float>& C) {
	omp_set_num_threads(num_Threads);

	time_t timestamp;
	time(&timestamp);
	std::string date = asctime(localtime(&timestamp));
	date.pop_back();
	date.append("_Call_Option_Price.csv");
	clock_t start, finish;

	for (std::string::iterator it = date.begin(); it<date.end(); ++it) {
		if (*it == ':') {
			date.erase(it);
		}
		std::replace(date.begin(), date.end(), ' ', '_');
	}

	FILE *f = fopen(date.c_str(), "w");


	fprintf(f, "%s;%s;%s;%s;%s;%s;\n", "S0", "TIME", "K", "r", "sig", "N");
	fprintf(f, "%lf;%lf;%lf;%lf;%lf;%i\n\n", pS0[0], pT[0], pK[0], r, sig, N);
	fprintf(f, "%s;%s;%s\n", "Index", "Option Value", "time (sec)");

		start = clock();
		baseVer(pT, pK, pS0, pC, N, r, sig);
		finish = clock();
		fprintf(f, "%i;%lf;%lf;\n", 0, pC[0], (float)(finish - start) / CLOCKS_PER_SEC);
		printf("Base\t\t%lf\t%lf;\n",pC[0], (float)(finish - start) / CLOCKS_PER_SEC);


		start = clock();
		eqBaseVer(pT, pK, pS0, pC, N, r, sig);
		finish = clock();
		fprintf(f, "%i;%lf;%lf;\n", 1, pC[0], (float)(finish - start) / CLOCKS_PER_SEC);
		printf("EqBase\t\t%lf\t%lf\n", pC[0], (float)(finish - start) / CLOCKS_PER_SEC);

		start = clock();
		updEqBaseVer(pT, pK, pS0, pC, N, r, sig);
		finish = clock();
		fprintf(f, "%i;%lf;%lf;\n", 2, pC[0], (float)(finish - start) / CLOCKS_PER_SEC);
		printf("d[2]EqBase\t%lf\t%lf\n", pC[0], (float)(finish - start) / CLOCKS_PER_SEC);

		start = clock();
		fpVer(pT, pK, pS0, pC, N, r, sig);
		finish = clock();
		fprintf(f, "%i;%lf;%lf;\n", 3, pC[0], (float)(finish - start) / CLOCKS_PER_SEC);
		printf("fp\t\t%lf\t%lf\n", pC[0], (float)(finish - start) / CLOCKS_PER_SEC);

		start = clock();
		fpVerOld(pT, pK, pS0, pC, N, r, sig);
		finish = clock();
		fprintf(f, "%i;%lf;%lf;\n", 3, pC[0], (float)(finish - start) / CLOCKS_PER_SEC);
		printf("fpOld\t\t%lf\t%lf\n", pC[0], (float)(finish - start) / CLOCKS_PER_SEC);

		start = clock();
		erfVer(pT, pK, pS0, pC, N, r, sig);
		finish = clock();
		fprintf(f, "%i;%lf;%lf;\n", 4, pC[0], (float)(finish - start) / CLOCKS_PER_SEC);
		printf("erf\t\t%lf\t%lf\n", pC[0], (float)(finish - start) / CLOCKS_PER_SEC);

		start = clock();
		erfVerOld(pT, pK, pS0, pC, N, r, sig);
		finish = clock();
		fprintf(f, "%i;%lf;%lf;\n", 4, pC[0], (float)(finish - start) / CLOCKS_PER_SEC);
		printf("erfOld\t\t%lf\t%lf\n", pC[0], (float)(finish - start) / CLOCKS_PER_SEC);

		start = clock();
		restrictVer(pT, pK, pS0, pC, N, r, sig);
		finish = clock();
		fprintf(f, "%i;%lf;%lf;\n", 5, pC[0], (float)(finish - start) / CLOCKS_PER_SEC);
		printf("restrict\t%lf\t%lf\n", pC[0], (float)(finish - start) / CLOCKS_PER_SEC);

		start = clock();
		vectVer(pT, pK, pS0, pC, N, r, sig);
		finish = clock();
		fprintf(f, "%i;%lf;%lf;\n", 6, pC[0], (float)(finish - start) / CLOCKS_PER_SEC);
		printf("vect\t\t%lf\t%lf\n", pC[0], (float)(finish - start) / CLOCKS_PER_SEC);

		start = clock();
		vectErfVer(pT, pK, pS0, pC, N, r, sig);
		finish = clock();
		fprintf(f, "%i;%lf;%lf;\n", 7, pC[0], (float)(finish - start) / CLOCKS_PER_SEC);
		printf("vectErfVer\t%lf\t%lf\n", pC[0], (float)(finish - start) / CLOCKS_PER_SEC);

		start = clock();
		vectInvSqrtVer(pT, pK, pS0, pC, N, r, sig);
		finish = clock();
		fprintf(f, "%i;%lf;%lf;\n", 8, pC[0], (float)(finish - start) / CLOCKS_PER_SEC);
		printf("vectInvSqrtVer\t%lf\t%lf\n", pC[0], (float)(finish - start) / CLOCKS_PER_SEC);

		start = clock();
		parVer(pT, pK, pS0, pC, N, r, sig);
		finish = clock();
		fprintf(f, "%i;%lf;%lf;\n", 9, pC[0], (float)(finish - start) / CLOCKS_PER_SEC);
		printf("parVer\t\t%lf\t%lf\n", pC[0], (float)(finish - start) / CLOCKS_PER_SEC);

		start = clock();
		parVerOld(pT, pK, pS0, pC, N, r, sig);
		finish = clock();
		fprintf(f, "%i;%lf;%lf;\n", 9, pC[0], (float)(finish - start) / CLOCKS_PER_SEC);
		printf("parVerOld\t%lf\t%lf\n", pC[0], (float)(finish - start) / CLOCKS_PER_SEC);

		start = clock();
		parNontempVer(pT, pK, pS0, pC, N, r, sig);
		finish = clock();
		fprintf(f, "%i;%lf;%lf;\n", 10, pC[0], (float)(finish - start) / CLOCKS_PER_SEC);
		printf("parNontempVer\t%lf\t%lf\n", pC[0], (float)(finish - start) / CLOCKS_PER_SEC);

		start = clock();
		stdVectorVer(T, K, S0, C, N, r, sig);
		finish = clock();
		fprintf(f, "%i;%lf;%lf;\n", 11, C[0], (float)(finish - start) / CLOCKS_PER_SEC);
		printf("stdVectorVer\t%lf\t%lf\n", pC[0], (float)(finish - start) / CLOCKS_PER_SEC);

	fprintf(f, "\n");
	fclose(f);

}

float AnSolution::writeToFile(int numVer, int num_Threads, int N, float* pT, float* pK, float* pS0, float* pC, float r, float sig, std::vector<float>& T, std::vector<float>& K, std::vector<float>& S0, std::vector<float>& C) {

	//times = new float[VersionCount];
	//OptionValue = new float[VersionCount];
	omp_set_num_threads(num_Threads);

	time_t timestamp;
	time(&timestamp);
	std::string date = asctime(localtime(&timestamp));
	date.pop_back();
	date.append("_Call_Option_Price.csv");
	clock_t start, finish;

	for (std::string::iterator it = date.begin(); it<date.end(); ++it) {
		if (*it == ':') {
			date.erase(it);
		}
		std::replace(date.begin(), date.end(), ' ', '_');
	}

	FILE *f = fopen(date.c_str(), "w");


	fprintf(f, "%s;%s;%s;%s;%s;%s;\n", "S0", "TIME", "K", "r", "sig", "N");
	fprintf(f, "%lf;%lf;%lf;%lf;%lf;%i\n\n", pS0[0], pT[0], pK[0], r, sig, N);
	fprintf(f, "%s;%s;%s\n", "Index", "Option Value", "time (sec)");
	switch (numVer) {
	case 0:
		start = clock();
		baseVer(pT, pK, pS0, pC, N, r, sig);
		finish = clock();
		fprintf(f, "%i;%lf;%lf;\n", 0, pC[0], (float)(finish - start) / CLOCKS_PER_SEC);
		return (float)(finish - start) / CLOCKS_PER_SEC;
		break;
	case 1:
		start = clock();
		fpVer(pT, pK, pS0, pC, N, r, sig);
		finish = clock();
		fprintf(f, "%i;%lf;%lf;\n", 1, pC[0], (float)(finish - start) / CLOCKS_PER_SEC);
		return (float)(finish - start) / CLOCKS_PER_SEC;
		break;
	case 2:
		start = clock();
		erfVer(pT, pK, pS0, pC, N, r, sig);
		finish = clock();
		fprintf(f, "%i;%lf;%lf;\n", 2, pC[0], (float)(finish - start) / CLOCKS_PER_SEC);
		return (float)(finish - start) / CLOCKS_PER_SEC;
		break;
	case 3:
		start = clock();
		restrictVer(pT, pK, pS0, pC, N, r, sig);
		finish = clock();
		fprintf(f, "%i;%lf;%lf;\n", 3, pC[0], (float)(finish - start) / CLOCKS_PER_SEC);
		return (float)(finish - start) / CLOCKS_PER_SEC;
		break;
	case 4:
		start = clock();
		vectVer(pT, pK, pS0, pC, N, r, sig);
		finish = clock();
		fprintf(f, "%i;%lf;%lf;\n", 4, pC[0], (float)(finish - start) / CLOCKS_PER_SEC);
		return (float)(finish - start) / CLOCKS_PER_SEC;
		break;
	case 5:
		start = clock();
		vectErfVer(pT, pK, pS0, pC, N, r, sig);
		finish = clock();
		fprintf(f, "%i;%lf;%lf;\n", 5, pC[0], (float)(finish - start) / CLOCKS_PER_SEC);
		return (float)(finish - start) / CLOCKS_PER_SEC;
		break;
	case 6:
		start = clock();
		vectInvSqrtVer(pT, pK, pS0, pC, N, r, sig);
		finish = clock();
		fprintf(f, "%i;%lf;%lf;\n", 6, pC[0], (float)(finish - start) / CLOCKS_PER_SEC);
		return (float)(finish - start) / CLOCKS_PER_SEC;
		break;
	case 7:
		start = clock();
		parVer(pT, pK, pS0, pC, N, r, sig);
		finish = clock();
		fprintf(f, "%i;%lf;%lf;\n", 7, pC[0], (float)(finish - start) / CLOCKS_PER_SEC);
		return (float)(finish - start) / CLOCKS_PER_SEC;
		break;
	case 8:
		start = clock();
		parNontempVer(pT, pK, pS0, pC, N, r, sig);
		finish = clock();
		fprintf(f, "%i;%lf;%lf;\n", 8, pC[0], (float)(finish - start) / CLOCKS_PER_SEC);
		return (float)(finish - start) / CLOCKS_PER_SEC;
		break;
	case 9:
		start = clock();
		stdVectorVer(T, K, S0, C, N, r, sig);
		finish = clock();
		fprintf(f, "%i;%lf;%lf;\n", 9, C[0], (float)(finish - start) / CLOCKS_PER_SEC);
		return (float)(finish - start) / CLOCKS_PER_SEC;
		break;
	}
	//start = finish = 0.0;
	//int sum = 0;
	//int a[1024] = { 0 };
	////#pragma omp parallel for shared(a) reduction (+: sum) 
	//{
	//	//# pragma omp for
	//	for (int i = 0; i < 4096; ++i)
	//		sum += a[i];
	//}
	//(this->*version_array[VersionCount - 1])(pT, pK, pS0, pC);
	//(this->*version_array[VersionCount - 1])(pT, pK, pS0, pC);
	//start = clock();
	//(this->*version_array[VersionCount - 1])(pT, pK, pS0, pC);
	//finish = clock();
	//fprintf(f, "%i;%lf;%lf;\n", VersionCount - 1, pC[0], (float)(finish - start) / CLOCKS_PER_SEC);

	fprintf(f, "\n");
	fclose(f);
}
