#include "../include/S30_NumMethods.h"

void NumMethods::WriteToCsv(double* Errors, int nSteps, int nRows, double Time, int scale, char* FileName) {
	
	row_table rt;
	time_t rawtime;
	time(&rawtime);
	std::string date = asctime(localtime(&rawtime));
	date.pop_back();
	date.append(FileName);

	for (std::string::iterator it = date.begin(); it<date.end(); ++it) {
		if (*it == ':') {
			date.erase(it);
		}
		std::replace(date.begin(), date.end(), ' ', '_');
	}
	

	FILE *f = fopen(date.c_str(), "w");
	fprintf(f, "h;e;log(h);log(e);\n");
	for (int i = 0; i < nRows; ++i) {
		rt.SetValues(Errors[i], nSteps, (int)(scale / pow(2,i)), nRows, Time);
		std::string tmp_string = rt.CreateString();
		for (std::string::iterator it = tmp_string.begin(); it<tmp_string.end(); ++it) {
			std::replace(tmp_string.begin(), tmp_string.end(), '.', ',');
		}
		fprintf(f, "%s\n", tmp_string.c_str());
	}
	fclose(f);
}