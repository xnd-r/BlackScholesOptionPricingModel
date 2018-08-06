#include "EM20_ChiSquared.h"

ChiSquared::ChiSquared(float* wd, int _len, float _h,int _k, float al) 
	: RVCharacteristics(wd, _len, _h) {
	k = _k;
	alpha = al;
}

void ChiSquared::SetIntervals() {
	if (k >= 15) {
		arrZ = new float[k - 1];
		arrZ[k - 2] = wd_sorted[len - 1];
		double d = (arrZ[k - 2] - wd_sorted[5]) / (k - 2);
		arrZ[0] = wd_sorted[5];
		for (int i = 1; i < k - 2; i++)
			arrZ[i] = arrZ[0] + d * i;
	}
}

void ChiSquared::SetNj() {
	int counter = 0;
	arrN = new float[len];
	for (int i = 0; i < len; ++i)
	{
		if (wd_sorted[i] <= arrZ[0])
		{
			counter++;
		}
	}
	arrN[0] = counter;
	counter = 0;

	for (int j = 1; j < k - 1; ++j) {
		counter = 0;
		for (int i = 0; i < len; ++i)
			if (arrZ[j - 1] < wd_sorted[i] && wd_sorted[i] <= arrZ[j])
				counter++;
		arrN[j] = counter;
		counter = 0;
	}

	for (int i = 0; i < len; ++i)
		if (arrZ[k - 2] < wd_sorted[i])
			counter++;

	arrN[k - 1] = counter;
}
float ChiSquared::CDF(float mean, float variance, float x) {
	return 0.5f * (1.0f + erff(-((x - mean) * (x - mean)) / (2.0f * variance)));
}

void ChiSquared::SetQj()
{
	arrQ[0] = CDF(MeanCh, VarCh, arrZ[0]);
	for (int i = 1; i < k - 1; ++i)
	{
		arrQ[i] = CDF(MeanCh, VarCh, arrZ[i]) - CDF(MeanCh, VarCh, arrZ[i - 1]);;
	}
	arrQ[k - 1] = 1.0f - CDF(MeanCh, VarCh, arrZ[k - 2]);
}

void ChiSquared::SetR0() {
	float res = 0.0f;
	for (int i = 0; i < k; ++i)
		res += powf(arrN[i] / (float)len - arrQ[i], 2) / arrQ[i];
	R0 = len * res;
}

float ChiSquared::ChiSquaredDencity(int r, float x) {
	return x > 0.0f ? powf(2.0f, -r / 2.0f) * powf(tgammaf(r), -1) : 0.0f;
}

float ChiSquared::ChiSquaredDistribute(int r) // integration of ChiSquaredDencity from 0.0 to R0
{
	double res = 0.0;
	int n = 1000; // amount of sections in integral calculating
	for (int k = 1; k <= n; ++k)
	{
		res += (ChiSquaredDencity(R0 * (k - 1) / n, r) + ChiSquaredDencity(R0 * k / n, r)) * (R0 / (2.0f * n));
	}
	return res;
}

char* ChiSquared::IsHypoAccepted() {
	return 1.0f - ChiSquaredDistribute(k - 1) < alpha ? "Hypothsesis accepted" : "Hypothsesis rejected";
}