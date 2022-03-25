#include "diff.h"



int main() {
	double x0 = 0;
	double y0 = 1;
	double step = 0.005;
	Function ex(x0, y0, step);
	ex.explicit_euler(true);
	system("python3 graph.py");
	Function ex2(x0, y0, step);
	ex2.implicit_euler(true);
	system("python3 graph2.py");
	Function ex3(x0, y0, step);
	ex3.recalc_euler(true);
	system("python3 graph3.py");
	Function ex4(x0, y0, step);
	ex4.runge_kutta(true);
	system("python3 graph4.py");
	count_error(x0, y0);
	system("python3 graph5.py");
	return 0;
}
