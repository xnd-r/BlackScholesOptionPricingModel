#include "E_RowTable.h"
#include <vector>
void row_table::setValues(double e, int nst, int sc, int nr, float t) {
	Error = e; 
	nsteps = nst; 
	scale = sc; 
	nrows = nr; 
	time = t;
}

std::string row_table::createString() {
	std::string res;

	res.append(std::to_string(time / static_cast<float>(nsteps / scale)));
	res.append(";");
	res.append(std::to_string(Error));
	res.append(";");
	res.append(std::to_string(log10(time / static_cast<float>(nsteps / scale))));
	res.append(";");
	res.append(std::to_string(log10(Error)));
	res.append(";");
	
	return res;
}