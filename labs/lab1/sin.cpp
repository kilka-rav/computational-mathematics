#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

class Function {
	private:
		std::string name;
		int number;
		int length;
		std::vector<double> data_x;
		std::vector<double> data;
		std::vector<double> data1;
		std::vector<double> data1_1;
		std::vector<double> data2;
		std::vector<double> data2_1;
	public:
		Function(std::string name, int N, double length);
		void record(std::string filename, std::string filename1, std::string filename2, std::string filename3);
		void count();
		void count2(double& e1, double& e2, double& e1_1, double& e2_1);
};

Function::Function(std::string _name, int _number, double _length ) : name(_name), number(_number), length(_length) {}

void _record(std::string filename, std::vector<double>& data) {
	std::ofstream out;
	out.open(filename + ".txt");
        if ( out.is_open() ) {
                for( auto&& t : data) {
                        out << t << std::endl;
                }
        }
        else {
                std::cout << "ERROR IN RECORD 1" << std::endl;
        }
}


void Function::record(std::string filename1, std::string filename2, std::string filename3, std::string filename4) {
	_record(filename1, data);
	_record(filename2, data1);
	_record(filename3, data2);
	_record(filename4, data_x);
}

void Function::count() {
	double pointer = -1 * length / 2.0;
	double end = -1 * pointer;
	double dx = length / (static_cast<double>(number));
	for(int i = 0; i < number; ++i) {
		data.push_back(sin(pointer));
		data_x.push_back(pointer);
		pointer += dx;
	}
	data1.push_back((4 * data[1] - 3 * data[0] - data[2]) / (2 * dx));
	data1_1.push_back((data[1] - data[0]) / dx);
	data2.push_back((2 * data[0] - 5 * data[1] + 4 * data[2] - data[3]) / (dx * dx));
	data2_1.push_back(data[0] - 2 * data[1] + data[2] / (dx * dx));
	pointer = -1 * length / 2.0;
	for(int i = 1; i < number - 1; ++i) {
		data1.push_back((data[i+1] - data[i-1]) / (2*dx));
		data2.push_back((data[i+1] + data[i-1] - 2 * data[i]) / (dx * dx));
		data1_1.push_back((data[i] - data[i - 1] ) / dx);
		data2_1.push_back((data[i - 1] - 2 * data[i] + data[i + 1] ) / (dx * dx));
		pointer += dx;
	}
	data1.push_back((data[number - 3] - 4 * data[number - 2] + 3 * data[number - 1]) / (2 * dx)); 
        data2.push_back((-1 * data[number - 4] + 4 * data[number - 3] - 5 * data[number - 2] + 2 *data[number - 1]) / (dx * dx));
	data1_1.push_back((data[number - 1] - data[number - 2]) / dx);
	data2_1.push_back((data[number - 3] - 2 * data[number - 2] + data[number - 1] ) / (dx * dx));
}

double modul(double a) {
	if ( a > 0.0 ) {
		return a;
	}
	return -1 * a;
}

void Function::count2(double& e1, double& e2, double& e1_1, double& e2_1) {
	double pointer = -1 * length / 2.0;
	double end = -1 * pointer;
	double dx = length / (static_cast<double>(number));
	for(int i = 0; i < number; ++i) {
		if ( modul(cos(pointer) - data1[i]) > e1 ) {
			e1 = modul(cos(pointer) - data1[i]);
			e1_1 = modul(cos(pointer) - data1_1[i]);
		}
		if ( modul(-1 * sin(pointer) - data2[i]) > e2 ) {
			e2 = modul(-1 * sin(pointer) - data2[i]);
			e2_1 = modul(-1 * sin(pointer) - data2_1[i]);
			//std::cout << "e2 = " << e2 << std::endl;
		}
		pointer += dx;
	}
}

void stepping() {
	std::vector<double> dx;
	std::vector<double> max_error;
	std::vector<double> max_error_1;
	std::vector<double> max_error_two;
	std::vector<double> max_error_two_1;
	int leng = 5;
	double e1, e2, e1_1, e2_1;
	for(int N = 33; N < 10000; N *= 2) {
		e1 = 0.0;
		e2 = 0.0;
		e1_1 = 0.0;
		e2_1 = 0.0;
		dx.push_back(static_cast<double>(leng) / static_cast<double>(N));
		Function s("sin", N, leng);
		s.count();
		s.count2(e1, e2, e1_1, e2_1);
		max_error.push_back(e1);
		max_error_1.push_back(e1_1);
		max_error_two.push_back(e2);
		max_error_two_1.push_back(e2_1);
	}
	_record("dx", dx);
	_record("error_1", max_error);
	_record("error_1_1", max_error_1);
	_record("error_2", max_error_two);
	_record("error_2_2", max_error_two_1);
}

int main() {
	Function sinus("sin", 20, 5);
	sinus.count();
	sinus.record("data1", "data2", "data3", "data_x");
	system("python3 graph.py");
	system("python3 graph2.py");
	stepping();
	system("python3 graph3.py");
	system("python3 graph4.py");
	return 0;
}
