#ifndef ____ANALITYCAL_SOLUTION_OPTION____
#define ____ANALITYCAL_SOLUTION_OPTION____
#include "../include/O00_EuropeanOption.h"


class AnSolutionOption : public EuropeanOption {
public:

	int		num_Threads;
	int 	version;
	double	_time;
	double	start, finish;

	float	d1, d2, erf1, erf2, erf3, erf4;
	float	sig2;
	float	invf;
	float	p1, p2;

	// TODO: try to make abstract
	//TGetOptionPrice option_array[9];
	virtual void _V0() {};
	virtual void _V1() {};
	virtual void _V2() {};
	virtual void _V3() {};
	virtual void _V4() {};
	virtual void _V5() {};
	virtual void _V6() {};
	virtual void _V7() {};
	virtual void _V8() {};
};
#endif // !____ANALITYCAL_SOLUTION_OPTION____