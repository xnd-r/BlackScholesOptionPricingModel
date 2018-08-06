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
	arrN = new float[k - 1];
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
	arrQ = new float[k - 1];

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

void ChiSquared::Execute() {
	SetIntervals();
	SetNj();
	SetQj();
	SetR0();
	ChiSquaredDistribute(k - 1);
}

void ChiSquared::WriteToCsv() {
		time_t rawtime;
		time(&rawtime);
		std::string date = asctime(localtime(&rawtime));
		date.pop_back();
		date.append("_GBM_Analitycal_Simple.csv");

		for (std::string::iterator it = date.begin(); it<date.end(); ++it) {
			if (*it == ':') {
				date.erase(it);
			}
			std::replace(date.begin(), date.end(), ' ', '_');
		}

		FILE *f = fopen(date.c_str(), "w");
		for (int j = 0; j < k; j++) {
			std::string tmp_cell = std::to_string(arrZ[j]);
			std::string tmp_cell2 = std::to_string(arrN[j]);
 			std::string tmp_cell3 = std::to_string(arrQ[j]);

			for (std::string::iterator it = tmp_cell.begin(); it<tmp_cell.end(); ++it)
				std::replace(tmp_cell.begin(), tmp_cell.end(), '.', ',');
			for (std::string::iterator it = tmp_cell2.begin(); it<tmp_cell2.end(); ++it)
				std::replace(tmp_cell2.begin(), tmp_cell2.end(), '.', ',');
			for (std::string::iterator it = tmp_cell3.begin(); it<tmp_cell3.end(); ++it)
				std::replace(tmp_cell3.begin(), tmp_cell3.end(), '.', ',');

			fprintf(f, /*"%lf;\n"*/"%s;%s;%s;", /*buffer[j]*/tmp_cell.c_str(), tmp_cell2.c_str(), tmp_cell3.c_str());
		}
		fprintf(f, "\n");
		fclose(f);
}