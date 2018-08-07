#include "O20_CallOption.h"

//base_version

void CallOption::_V0(float *pT, float *pK, float *pS0, float *pC)
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
	}
}

//fp

void CallOption::_V1(float *pT, float *pK, float *pS0, float *pC)
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
	}
}

//erf

void CallOption::_V2(float *pT, float *pK, float *pS0, float *pC)
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
	}
}

//restrict
#if defined (__INTEL_COMPILER)
void CallOption::_V3(float* restrict pT, float* restrict pK, float* restrict pS0, float* restrict pC)
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
	}
}
#endif
//fp_simd_vector_always

void CallOption::_V4(float *pT, float *pK, float *pS0, float *pC)
{

#pragma simd
#pragma vector always	
	for (int i = 0; i < N; i++)
	{
		d1 = (logf(pS0[i] / pK[i]) + (R + SIG * SIG * 0.5f) *
			pT[i]) / (SIG * sqrtf(pT[i]));
		d2 = (logf(pS0[i] / pK[i]) + (R - SIG * SIG * 0.5f) *
			pT[i]) / (SIG * sqrtf(pT[i]));

		vsCdfNorm(1, &d1, &p1);
		vsCdfNorm(1, &d2, &p2);
		pC[i] = pS0[i] * p1 - pK[i] * expf((-1.0f) * R * pT[i]) * p2;
	}
}

//fp_erf_simd_vector_always

void CallOption::_V5(float* pT, float* pK, float* pS0, float* pC)
{

#pragma simd
#pragma vector always
	for (int i = 0; i < N; i++)
	{
		d1 = (logf(pS0[i] / pK[i]) + (R + SIG * SIG * 0.5f) *
			pT[i]) / (SIG * sqrtf(pT[i]));
		d2 = (logf(pS0[i] / pK[i]) + (R - SIG * SIG * 0.5f) *
			pT[i]) / (SIG * sqrtf(pT[i]));
		erf1 = 0.5f + 0.5f * erff(d1 / sqrtf(2.0f));
		erf2 = 0.5f + 0.5f * erff(d2 / sqrtf(2.0f));

		pC[i] = pS0[i] * erf1 - pK[i] * expf((-1.0f) * R * pT[i]) * erf2;
	}
}

//using_defined_INVSQRT2

void CallOption::_V6(float* pT, float* pK, float* pS0, float* pC)
{

#pragma simd
	for (int i = 0; i < N; i++)
	{
		d1 = (logf(pS0[i] / pK[i]) + (R + SIG * SIG * 0.5f) * pT[i]) / (SIG * sqrtf(pT[i]));
		d2 = (logf(pS0[i] / pK[i]) + (R - SIG * SIG * 0.5f) * pT[i]) / (SIG * sqrtf(pT[i]));

		erf1 = 0.5f + 0.5f * erff(d1 * INVSQRT2);
		erf2 = 0.5f + 0.5f * erff(d2 * INVSQRT2);

		pC[i] = pS0[i] * erf1 - pK[i] * expf((-1.0f) * R * pT[i]) * erf2;
	}
}

//omp

void CallOption::_V7(float *pT, float *pK, float *pS0, float *pC)
{

#pragma simd
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
	}
}

//omp_nontemporal

void CallOption::_V8(float *pT, float *pK, float *pS0, float *pC)
{

#pragma simd 
#pragma vector nontemporal
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
	}
}