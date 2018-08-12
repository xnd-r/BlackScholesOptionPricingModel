#include "EM20_ChiSquared.h"
#define M_PIF 3.14159265358979323846f

ChiSquared::ChiSquared(float* wd, int _len, float _h,int _k, float al) 
	: RVCharacteristics(wd, _len, _h) {
	k = _k;
	alpha = al;
}

void ChiSquared::SetIntervals() {
	arrZ = new float[k - 1];
	arrZ[k - 2] = wd_sorted[len - 1];

	if (k >= 15) {
		float d = (arrZ[k - 2] - wd_sorted[5]) / (k - 2);
		arrZ[0] = wd_sorted[5];
		for (int i = 1; i < k - 2; i++)
			arrZ[i] = arrZ[0] + d * i;
	}
	else {
		float d = (arrZ[k - 2] - wd_sorted[1]) / (k - 2);
		arrZ[0] = wd_sorted[1];
		for (int i = 1; i < k - 2; i++)
			arrZ[i] = arrZ[0] + d * i;
	}
}

void ChiSquared::SetNj() {
	int counter = 0;
	arrN = new int[k - 1];
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
	return 0.5f * (1.0f + erff((x - mean) / sqrtf(2.0f * variance)));
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

float Gamma(int x)
{
	if (x == 2)	
		return 1.0f;
	else if (x == 1) 
		return sqrtf(M_PIF);
	else 
		return (x / 2.0f - 1.0f) * Gamma(x - 2);
}


float ChiSquared::ChiSquaredDencity(float x) {
	return x > 0.0f ? powf(2.0f, -(k - 1) / 2.0f) * powf(tgammaf((k - 1) / 2.0f), -1) * powf(x, (float)(k - 1) / 2.0f - 1.0f) * expf(-x / 2.0f) : 0.0f;
}

void ChiSquared::ChiSquaredDistribute() // integration of ChiSquaredDencity from 0.0 to R0
{
	float res = 0.0;
	int n = 1000; // amount of sections in integral calculating
	for (int i = 1; i <= n; ++i)
		res += (ChiSquaredDencity(R0 * (i - 1) / (float)n) + ChiSquaredDencity(R0 * i / (float)n)) * (R0 / (2.0f * (float)n));
	FDash = res;
}

char* ChiSquared::IsHypoAccepted() {
	return 1.0f - FDash < alpha ? "Hypothsesis accepted" : "Hypothsesis rejected";
}

void ChiSquared::Calculate() {
	SetIntervals();
	SetNj();
	SetQj();
	SetR0();
	ChiSquaredDistribute();
}

void ChiSquared::Execute() {
	Calculate();
	WriteToCsv();
}

void ChiSquared::WriteToCsv() {
		time_t rawtime;
		time(&rawtime);
		std::string date = asctime(localtime(&rawtime));
		date.pop_back();
		date.append("_Chi_Squared_Results.csv");

		for (std::string::iterator it = date.begin(); it<date.end(); ++it) {
			if (*it == ':') {
				date.erase(it);
			}
			std::replace(date.begin(), date.end(), ' ', '_');
		}

		FILE *f = fopen(date.c_str(), "w");
		fprintf(f, "%s\n", "index(j);Z[j];N[j];Q[j]");
		fprintf(f, "%s%lf%s;%i\n", "0;[-inf, ", arrZ[0], "]", arrN[0]);
		for (int j = 1; j < k - 1; j++) {
			std::string tmp_cell0 = std::to_string(arrZ[j - 1]);
			std::string tmp_cell1 = std::to_string(arrZ[j]);
			std::string tmp_cell2 = std::to_string(arrN[j]);
 			std::string tmp_cell3 = std::to_string(arrQ[j]);
			for (std::string::iterator it = tmp_cell0.begin(); it<tmp_cell0.end(); ++it)
				std::replace(tmp_cell0.begin(), tmp_cell0.end(), '.', ',');
			for (std::string::iterator it = tmp_cell1.begin(); it<tmp_cell1.end(); ++it)
				std::replace(tmp_cell1.begin(), tmp_cell1.end(), '.', ',');
			for (std::string::iterator it = tmp_cell2.begin(); it<tmp_cell2.end(); ++it)
				std::replace(tmp_cell2.begin(), tmp_cell2.end(), '.', ',');
			for (std::string::iterator it = tmp_cell3.begin(); it<tmp_cell3.end(); ++it)
				std::replace(tmp_cell3.begin(), tmp_cell3.end(), '.', ',');

			fprintf(f, "%i;%s%s%s%s%s;%s;%s;\n", j, "[", tmp_cell0.c_str(), ", ", tmp_cell1.c_str(), "]", tmp_cell2.c_str(), tmp_cell3.c_str());
		}
		fprintf(f, "%i;%s%lf%s;%i\n\n", k - 1, "[", arrZ[k - 2], ", +inf]", arrN[k - 1]);
		fprintf(f, "%s%lf;\n", "R0: ",R0);
		fprintf(f, "%s%lf;\n", "Integrated Chi Squared Dencity from 0.0 to R0: ", FDash);
		fprintf(f, "%s%lf;\n", "Significance level: ", alpha);
		fprintf(f, "%s\n", IsHypoAccepted());

		fclose(f);
}