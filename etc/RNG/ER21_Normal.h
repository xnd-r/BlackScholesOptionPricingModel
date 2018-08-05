#ifndef ____NORMAL_GENERATOR____
#define ____NORMAL_GENERATOR____

#include <math.h>
#include "ER10_RNG.h"

class  NormalGen : public MCG59 {

	const float p0 = 0.322232431088f;     const float q0 = 0.099348462606f;
	const float p1 = 1.0f;                const float q1 = 0.588581570495f;
	const float p2 = 0.342242088547f;     const float q2 = 0.531103462366f;
	const float p3 = 0.204231210245e-1f;  const float q3 = 0.103537752850f;
	const float p4 = 0.453642210148e-4f;  const float q4 = 0.385607006340e-2f;

public:
	NormalGen(const float _mean, const float _var);
	float mean;
	float variance;
	float u, t, p, q, z;

	void RandomArray(float* dest_arr, int len);
	float GetFloat();
};

#endif // !____NORMAL_GENERATOR____