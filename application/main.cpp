#include "../modules/StockPriece/include/S10_Wiener.h"
#include "../modules/StockPriece/include/S20_AnSimple.h"
#include "../modules/StockPriece/include/S21_AnExtended.h"

#include <iostream>

int main() { 
	StockPrice sp;
	sp.Execute();
	std::cout << "Wiener Success" << std::endl;

	AnSimple as;
	as.Execute();
	std::cout << "AnSimple Success" << std::endl;
	
	AnExtended ae;
	ae.Execute();
	std::cout << "AnExtended Success" << std::endl;

	system("pause");
	return 0; 
}