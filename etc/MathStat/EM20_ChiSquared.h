#ifndef ____CHI_SQUARED_DISTRIBUTION____
#define ____CHI_SQUARED_DISTRIBUTION____

#include "EM10_MathStat.h"

class ChiSquared : public RVCharacteristics {
	int		k;		// amount of sections to calculate R0
	float	R0;		// statustics for ChiSquared distribution
	float	alpha;	// significance level 
	float*	arrZ;
	int*	arrN;
	float*	arrQ;
	float	FDash;

public:
	ChiSquared(float* wd, int _len, float _h, int _k, float al);
	float	CDF(float mean, float variance, float x); // cummulative distribution function for normal distribution
	float	ChiSquaredDencity(float x);
	void	ChiSquaredDistribute();
	void	SetIntervals();
	void	SetNj();
	void	SetQj();
	void	SetR0();
	void	WriteToCsv();
	void	Execute();
	char*	IsHypoAccepted();


};

#endif // !____CHI_SQUARED_DISTRIBUTION____