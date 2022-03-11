#include "solution.h"
#include <iostream>



int main() {
	Functions f(0.001);
	f.divide();
	f.chord();
	f.secants();
	f.Newton();

	return 0;
}
