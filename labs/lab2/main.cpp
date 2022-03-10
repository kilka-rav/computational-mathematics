#include "matrix.h"



int main() {
	int size;
	std::cout << "Input dimension: ";
	std::cin >> size;
	double** arr = nullptr;
    	arr = init(arr, size, size + 1);
	std::cout << "Input system of linear equations\n";
	double x;
	for(int i = 0; i < size; ++i) {
		for(int j = 0; j < size + 1; ++j) {
			std::cin >> x;
			arr[i][j] = x;
		}
	}
	Matrix A(size, size + 1, arr);
	Matrix B(size, size + 1, arr);
	Matrix C(size, size + 1, arr);
	Matrix D(size, size + 1, arr);
	Matrix E(size, size + 1, arr);
	Matrix F(size, size + 1, arr);
	Matrix G(size, size + 1, arr);
	std::cout << "count of inconsistencies = " << E.count_incon() << std::endl;
	//E.count_incon();
	//std::cout << "mul = " << E.scalar_multiply(E) << std::endl;
	A.gauss();
	std::cout << "Gauss:\n";
	A.print_reply();
	B.gauss_main_elem();
	std::cout << "Gauss + choice\n";
	B.print_reply();
	bool flag = C.shuttle();
	if ( flag == false ) {
		std::cout << "ERROR IN DIAG\n";
	}
	std::cout << "Yacobi\n";		//add_checking
	D.yacobi(0.0001);
	std::cout << "Zeydel\n";		//add_checking
	E.zeydel(0.0001);
	std::cout << "Gradient down\n";		//add checking
	G.gradient_down(0.0001);
	std::cout << "the smallest discrepancy\n";
	F.discrepancy(0.0001);			//add check
	return 0;
}
