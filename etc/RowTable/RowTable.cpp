#include "RowTable.h"

void row_table::SetValues(double e, int nst, int sc, int nr, double t) {
	error = e; 
	nSteps = nst; 
	scale = sc; 
	nRows = nr; 
	time = t;
}

std::string row_table::CreateString() {
	std::string res;

	res.append(std::to_string(time/(double)(nSteps / scale)));
	res.append(";");
	res.append(std::to_string(error));
	res.append(";");
	res.append(std::to_string(log10(time / (double)(nSteps / scale))));
	res.append(";");
	res.append(std::to_string(log10(error)));
	res.append(";");
	
	return res;
}