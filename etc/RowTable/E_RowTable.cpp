#include "E_RowTable.h"

void row_table::SetValues(float e, int nst, int sc, int nr, float t) {
	error = e; 
	nSteps = nst; 
	scale = sc; 
	nRows = nr; 
	time = t;
}

std::string row_table::CreateString() {
	std::string res;

	res.append(std::to_string(time/(float)(nSteps / scale)));
	res.append(";");
	res.append(std::to_string(error));
	res.append(";");
	res.append(std::to_string(log10(time / (float)(nSteps / scale))));
	res.append(";");
	res.append(std::to_string(log10(error)));
	res.append(";");
	
	return res;
}