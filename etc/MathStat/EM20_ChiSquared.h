#ifndef ____CHI_SQUARED_DISTRIBUTION____
#define ____CHI_SQUARED_DISTRIBUTION____

#include "EM10_MathStat.h"

class ChiSquared : public RVCharacteristics {
	int k; // amount of sections to calculate R0
	float R0; // statustics for ChiSquared distribution
	float alpha; // significance level 
	float* arrZ;
	float* arrN;
	float* arrQ;
public:
	ChiSquared(float* wd, int _len, float _h, int _k, float al);
	float CDF(float mean, float variance, float x); // cummulative distribution function for normal distribution
	float ChiSquaredDencity(int r, float x);
	float ChiSquared::ChiSquaredDistribute(int r);
	void SetNj();
	void SetQj();
	void SetR0();
	void SetIntervals();
	char* IsHypoAccepted();
	void WriteToCsv();
	void Execute();

};

#endif // !____CHI_SQUARED_DISTRIBUTION____