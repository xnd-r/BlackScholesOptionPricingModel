#ifndef ____CALL_OPTION____
#define ____CALL_OPTION____

#include "O10_AnSolution.h"

class  CallOption : public AnSolutionOption {
public:
	typedef void(CallOption::*TGetOptionPrice)(float *pT, float *pK, float *pS0, float *pC);

	void _V0(float *pT, float *pK, float *pS0, float *pC);
	void _V1(float *pT, float *pK, float *pS0, float *pC);
	void _V2(float *pT, float *pK, float *pS0, float *pC);

#if defined (__INTEL_COMPILER)
	void _V3(float* restrict pT, float* restrict pK, flot pS0, float* restrict pC);
#endif

	void _V4(float *pT, float *pK, float *pS0, float *pC);
	void _V5(float *pT, float *pK, float *pS0, float *pC);
	void _V6(float *pT, float *pK, float *pS0, float *pC);
	void _V7(float *pT, float *pK, float *pS0, float *pC);
	void _V8(float *pT, float *pK, float *pS0, float *pC);
	void WriteToCsv(int threads);

	CallOption::TGetOptionPrice version_array[8] = 
	{ 
		&CallOption::_V0, // preference 1
		&CallOption::_V1, // preference 2
		&CallOption::_V2, // erf
		&CallOption::_V4, // vect without restrict
		&CallOption::_V5, // #pragma simd #pragma vector always
		&CallOption::_V6, // #pragma simd invsqrt2_1
		&CallOption::_V7, // #pragma simd #pragma omp parallel for private
		&CallOption::_V8  // _V7 + #pragma vector nontemporal
	};

#if defined (__INTEL_COMPILER)
	CallOption::TGetOptionPrice version_array[9] = 
	{ 
		&CallOption::_V0,
		&CallOption::_V1,
		&CallOption::_V2,
		&CallOption::_V3, // restrict
		&CallOption::_V4,
		&CallOption::_V5,
		&CallOption::_V6, 
		&CallOption::_V7,
		&CallOption::_V8 
	};
#endif
};

#endif // !____CALL_OPTION____
