#include "ER21_Normal.h"

NormalGen::NormalGen(const float _m, const float _v, long int _seed) : MCG59(_seed){
	mean = _m;
	variance = _v;
}

float NormalGen::GetFloat() {

	u = MCG59::GetFloat();
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
	return 0;
}

void NormalGen::RandomArray(/*float* input_arr,*/ float* dest_arr, int len) {
	for (int i = 0; i < len; ++i) {
		dest_arr[i] = GetFloat();
		//if (input_arr[i] < 0.5f)
		//	t = sqrtf(-2.0f * logf(input_arr[i]));
		//else
		//	t = sqrtf(-2.0f * logf(1.0f - input_arr[i]));
		//p = p0 + t * (p1 + t * (p2 + t * (p3 + t * p4)));
		//q = q0 + t * (q1 + t * (q2 + t * (q3 + t * q4)));
		//if (input_arr[i] < 0.5f)
		//	z = (p / q) - t;
		//else
		//	z = t - (p / q);
		//dest_arr[i] = (mean + sqrtf(variance) * z);
	}
}
