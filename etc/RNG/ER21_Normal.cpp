#include "ER21_Normal.h"

NormalGen::NormalGen(const float _m, const float _v) {
	mean = _m;
	variance = _v;
}

float NormalGen::GetFloat() {

	u = GetFloat();
	if (u < 0.5f)
		t = sqrtf(-2.0f * logf(u));
	else
		t = sqrtf(-2.0f * logf(1.0f - u));
	p = p0 + t * (p1 + t * (p2 + t * (p3 + t * p4)));
	q = q0 + t * (q1 + t * (q2 + t * (q3 + t * q4)));
	if (u < 0.5f)
		z = (p / q) - t;
	else
		z = t - (p / q);
	return (mean + sqrtf(variance) * z);
}

void NormalGen::RandomArray(float* dest_arr, int len) {
	for (int i = 0; i < len; ++i)
		dest_arr[i] = GetFloat();
}
