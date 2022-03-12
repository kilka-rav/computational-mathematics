#include "interpol.hpp"



int main() {
	Interpol ex(1, 18.5, 8, 2.5);
	ex.Newton();
	Interpol ex2(1, 16, 7, 2.5);
	ex2.Lagrange();
	Interpol ex3(1, 16, 7, 2.5);
	ex3.cube();
	system("python3 graph.py");
	system("python3 graph2.py");
	system("python3 graph3.py");
	return 0;
}
	
