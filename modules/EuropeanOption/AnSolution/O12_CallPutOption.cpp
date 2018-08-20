#include "O12_CallPutOption.h"

//base_version

void CallPutOption::_V0(float *pT, float *pK, float *pS0, float *pC, float *pP)
{

	for (int i = 0; i < N; i++)
	{
		d1 = (log(pS0[i] / pK[i]) + (R + SIG * SIG * 0.5) *
			pT[i]) / (SIG * sqrt(pT[i]));
		d2 = (log(pS0[i] / pK[i]) + (R - SIG * SIG * 0.5) *
			pT[i]) / (SIG * sqrt(pT[i]));
		vsCdfNorm(1, &d1, &p1);
		vsCdfNorm(1, &d2, &p2);

		pC[i] = pS0[i] * p1 - pK[i] * exp((-1.0) * R * pT[i]) * p2;
		pP[i] = pC[i] - pS0[i] + pK[i] * exp((-1.0) * R * pT[i]);
	}
}

//fp

void CallPutOption::_V1(float *pT, float *pK, float *pS0, float *pC, float *pP)
{

	for (int i = 0; i < N; i++)
	{
		d1 = (logf(pS0[i] / pK[i]) + (R + SIG * SIG * 0.5f) *
			pT[i]) / (SIG * sqrtf(pT[i]));
		d2 = (logf(pS0[i] / pK[i]) + (R - SIG * SIG * 0.5f) *
			pT[i]) / (SIG * sqrtf(pT[i]));

		vsCdfNorm(1, &d1, &p1);
		vsCdfNorm(1, &d2, &p2);

		pC[i] = pS0[i] * p1 - pK[i] * expf((-1.0f) * R * pT[i]) * p2;
		pP[i] = pC[i] - pS0[i] + pK[i] * expf((-1.0f) * R * pT[i]);
	}
}

//erf

void CallPutOption::_V2(float *pT, float *pK, float *pS0, float *pC, float *pP)
{

	for (int i = 0; i < N; i++)
	{
		d1 = (logf(pS0[i] / pK[i]) + (R + SIG * SIG * 0.5f) *
			pT[i]) / (SIG * sqrtf(pT[i]));
		d2 = (logf(pS0[i] / pK[i]) + (R - SIG * SIG * 0.5f) *
			pT[i]) / (SIG * sqrtf(pT[i]));

		erf1 = 0.5f + 0.5f * erff(d1 / sqrtf(2.0f));
		erf2 = 0.5f + 0.5f * erff(d2 / sqrtf(2.0f));

		pC[i] = pS0[i] * erf1 - pK[i] * expf((-1.0f) * R * pT[i]) * erf2;
		pP[i] = pC[i] - pS0[i] + pK[i] * expf((-1.0f) * R * pT[i]);
	}
}

//restrict
#if defined (__INTEL_COMPILER)
void CallPutOption::_V3(float* restrict pT, float* restrict pK, float* restrict pS0, float* restrict pC, float* restrict pP)
{

	for (int i = 0; i < N; i++)
	{
		d1 = (logf(pS0[i] / pK[i]) + (r + sig * sig * 0.5f) *
			pT[i]) / (sig * sqrtf(pT[i]));
		d2 = (logf(pS0[i] / pK[i]) + (r - sig * sig * 0.5f) *
			pT[i]) / (sig * sqrtf(pT[i]));

		erf1 = 0.5f + 0.5f * erff(d1 / sqrtf(2.0f));
		erf2 = 0.5f + 0.5f * erff(d2 / sqrtf(2.0f));

		pC[i] = pS0[i] * erf1 - pK[i] * expf((-1.0f) * r * pT[i]) * erf2;
		pP[i] = pC[i] - S0[i] + pK[i] * expf((-1.0f) * r * pT[i]);
	}
}
#endif
//fp_simd_vector_always

void CallPutOption::_V4(float *pT, float *pK, float *pS0, float *pC, float *pP)
{
#if defined (__INTEL_COMPILER)
#pragma simd
#pragma vector always	
#endif
	for (int i = 0; i < N; i++)
	{
		d1 = (logf(pS0[i] / pK[i]) + (R + SIG * SIG * 0.5f) *
			pT[i]) / (SIG * sqrtf(pT[i]));
		d2 = (logf(pS0[i] / pK[i]) + (R - SIG * SIG * 0.5f) *
			pT[i]) / (SIG * sqrtf(pT[i]));

		vsCdfNorm(1, &d1, &p1);
		vsCdfNorm(1, &d2, &p2);
		pC[i] = pS0[i] * p1 - pK[i] * expf((-1.0f) * R * pT[i]) * p2;
		pP[i] = pC[i] - pS0[i] + pK[i] * expf((-1.0f) * R * pT[i]);
	}
}

//fp_erf_simd_vector_always

void CallPutOption::_V5(float* pT, float* pK, float* pS0, float* pC, float *pP)
{
#if defined (__INTEL_COMPILER)
#pragma simd
#pragma vector always
#endif
	for (int i = 0; i < N; i++)
	{
		d1 = (logf(pS0[i] / pK[i]) + (R + SIG * SIG * 0.5f) *
			pT[i]) / (SIG * sqrtf(pT[i]));
		d2 = (logf(pS0[i] / pK[i]) + (R - SIG * SIG * 0.5f) *
			pT[i]) / (SIG * sqrtf(pT[i]));
		erf1 = 0.5f + 0.5f * erff(d1 / sqrtf(2.0f));
		erf2 = 0.5f + 0.5f * erff(d2 / sqrtf(2.0f));

		pC[i] = pS0[i] * erf1 - pK[i] * expf((-1.0f) * R * pT[i]) * erf2;
		pP[i] = pC[i] - pS0[i] + pK[i] * expf((-1.0f) * R * pT[i]);
	}
}

//using_defined_invsqrt2

void CallPutOption::_V6(float* pT, float* pK, float* pS0, float* pC, float *pP)
{

#if defined (__INTEL_COMPILER)
#pragma simd
#endif
	for (int i = 0; i < N; i++)
	{
		d1 = (logf(pS0[i] / pK[i]) + (R + SIG * SIG * 0.5f) * pT[i]) / (SIG * sqrtf(pT[i]));
		d2 = (logf(pS0[i] / pK[i]) + (R - SIG * SIG * 0.5f) * pT[i]) / (SIG * sqrtf(pT[i]));

		erf1 = 0.5f + 0.5f * erff(d1 * INVSQRT2);
		erf2 = 0.5f + 0.5f * erff(d2 * INVSQRT2);

		pC[i] = pS0[i] * erf1 - pK[i] * expf((-1.0f) * R * pT[i]) * erf2;
		pP[i] = pC[i] - pS0[i] + pK[i] * expf((-1.0f) * R * pT[i]);
	}
}

//omp

void CallPutOption::_V7(float *pT, float *pK, float *pS0, float *pC, float *pP)
{

#if defined (__INTEL_COMPILER)
#pragma simd
#endif
#pragma omp parallel for private(invf, d1, d2, erf1, erf2)
	for (int i = 0; i < N; i++)
	{
		d1 = (logf(pS0[i] / pK[i]) + (R + SIG * SIG * 0.5f) *
			pT[i]) / (SIG * sqrtf(pT[i]));
		d2 = (logf(pS0[i] / pK[i]) + (R - SIG * SIG * 0.5f) *
			pT[i]) / (SIG * sqrtf(pT[i]));

		erf1 = 0.5f + 0.5f * erff(d1 / sqrtf(2.0f));
		erf2 = 0.5f + 0.5f * erff(d2 / sqrtf(2.0f));

		pC[i] = pS0[i] * erf1 - pK[i] * expf((-1.0f) * R * pT[i]) * erf2;
		pP[i] = pC[i] - pS0[i] + pK[i] * expf((-1.0f) * R * pT[i]);
	}
}

//omp_nontemporal

void CallPutOption::_V8(float *pT, float *pK, float *pS0, float *pC, float *pP)
{

#if defined (__INTEL_COMPILER)
#pragma simd
#pragma vector nontemporal
#endif

#pragma omp parallel for private(d1, d2, erf1, erf2)
	for (int i = 0; i < N; i++)
	{
		d1 = (logf(pS0[i] / pK[i]) + (R + SIG * SIG * 0.5f) *
			pT[i]) / (SIG * sqrtf(pT[i]));
		d2 = (logf(pS0[i] / pK[i]) + (R - SIG * SIG * 0.5f) *
			pT[i]) / (SIG * sqrtf(pT[i]));

		erf1 = 0.5f + 0.5f * erff(d1 / sqrtf(2.0f));
		erf2 = 0.5f + 0.5f * erff(d2 / sqrtf(2.0f));

		pC[i] = pS0[i] * erf1 - pK[i] * expf((-1.0f) * R * pT[i]) * erf2;
		pP[i] = pC[i] - pS0[i] + pK[i] * expf((-1.0f) * R * pT[i]);
	}
}

void CallPutOption::WriteToCsv(int threads) {
	int VersionCount;
#if defined (__INTEL_COMPILER)
	VersionCount = 9;
#else 
	VersionCount = 8;
#endif
	//times = new double[VersionCount];
	//OptionValue = new double[VersionCount];

	num_Threads = threads;
	omp_set_num_threads(num_Threads);

	float* pT = new float[5 * N];
	float* pK = pT + N;
	float* pS0 = pT + 2 * N;
	float* pC = pT + 3 * N;
	float* pP = pT + 4 * N;

#if defined(__INTEL_COMPILER) 
#pragma simd
#pragma vector always	
#endif
	for (int i = 0; i < N; i++) {
		pT[i] = TIME;
		pS0[i] = S0;
		pK[i] = K;
	}

	time_t rawtime;
	time(&rawtime);
	std::string date = asctime(localtime(&rawtime));
	date.pop_back();
	date.append("_Call_And_Put_Option_Price.csv");

	for (std::string::iterator it = date.begin(); it<date.end(); ++it) {
		if (*it == ':') {
			date.erase(it);
		}
		std::replace(date.begin(), date.end(), ' ', '_');
	}

	FILE *f = fopen(date.c_str(), "w");
	fprintf(f, "%s;%s;%s;%s;%s;%s;\n", "S0", "TIME", "K", "R", "SIG", "N");
	fprintf(f, "%lf;%lf;%lf;%lf;%lf;%i\n\n", S0, TIME, K, R, SIG, N);
	fprintf(f, "%s;%s;%s%s;\n", "Index", "Call Option Value", "Put Option Value", "Time (sec)");
	for (int i = 0; i < VersionCount - 1; i++) {
		start = clock();
		(this->*version_array[i])(pT, pK, pS0, pC, pP);
		finish = clock();

		fprintf(f, "%i;%lf;%lf;%lf\n", i, pC[0], pP[0], (double)(finish - start) / CLOCKS_PER_SEC);
	}
	start = finish = 0.0;
	int sum = 0;
	int a[1024] = { 0 };
#pragma omp parallel for shared(a) reduction (+: sum) 
	{
# pragma omp for
		for (int i = 0; i < 4096; ++i)
			sum += a[i];
	}
	(this->*version_array[VersionCount - 1])(pT, pK, pS0, pC, pP);
	(this->*version_array[VersionCount - 1])(pT, pK, pS0, pC, pP);
	start = clock();
	(this->*version_array[VersionCount - 1])(pT, pK, pS0, pC, pP);
	finish = clock();
	fprintf(f, "%i;%lf;%lf;%lf\n", VersionCount - 1, pC[0], pP[0], (double)(finish - start) / CLOCKS_PER_SEC);

	fprintf(f, "\n");

	fclose(f);
}
