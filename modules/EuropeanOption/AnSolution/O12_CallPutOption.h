#ifndef ____CALL_AND_PUT_OPTION____
#define ____CALL_AND_PUT_OPTION____

#include "O10_AnSolution.h"

class  CallPutOption : public AnSolutionOption {
public:
	typedef void(CallPutOption::*TGetOptionPrice)(float *pT, float *pK, float *pS0, float *pC, float* pP);

	void _V0(float *pT, float *pK, float *pS0, float *pC, float* pP);
	void _V1(float *pT, float *pK, float *pS0, float *pC, float* pP);
	void _V2(float *pT, float *pK, float *pS0, float *pC, float* pP);

#if defined (__INTEL_COMPILER)
	void _V3(float* __restrict pT, float* __restrict pK, float* __restrict pS0, float* __restrict pC, float* pP);
#endif

	void _V4(float *pT, float *pK, float *pS0, float *pC, float* pP);
	void _V5(float *pT, float *pK, float *pS0, float *pC, float* pP);
	void _V6(float *pT, float *pK, float *pS0, float *pC, float* pP);
	void _V7(float *pT, float *pK, float *pS0, float *pC, float* pP);
	void _V8(float *pT, float *pK, float *pS0, float *pC, float* pP);
	void WriteToCsv(int threads);

#if !defined (__INTEL_COMPILER)
	CallPutOption::TGetOptionPrice version_array[8] =
	{
		&CallPutOption::_V0, // preference 1
		&CallPutOption::_V1, // preference 2
		&CallPutOption::_V2, // erf
		&CallPutOption::_V4, // vect without restrict
		&CallPutOption::_V5, // #pragma simd #pragma vector always
		&CallPutOption::_V6, // #pragma simd invsqrt2_1
		&CallPutOption::_V7, // #pragma simd #pragma omp parallel for private
		&CallPutOption::_V8  // _V7 + #pragma vector nontemporal
	};
#endif

#if defined (__INTEL_COMPILER)
	CallPutOption::TGetOptionPrice version_array[9] =
	{
		&CallPutOption::_V0,
		&CallPutOption::_V1,
		&CallPutOption::_V2,
		&CallPutOption::_V3, // restrict
		&CallPutOption::_V4,
		&CallPutOption::_V5,
		&CallPutOption::_V6,
		&CallPutOption::_V7,
		&CallPutOption::_V8
	};
#endif
};

#endif // !____CALL_AND_PUT_OPTION____
