#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

class Interpol {
	private:
		double left;
		double right;
		double len;
		double step;
		int N;
	public:
		std::vector<double> X;
		std::vector<double> Y;
		std::vector<double> D;
	public:
		Interpol(double l, double r, int _N, double s);
		double F(double x) const;
		void Newton();
		double polynom_newton(double x);
		void record_newton();
		void Lagrange();
		double l(double x, int num);
		double polynom_lagrange(double x);
		void cube();
		double cub_spline(double x);
};

void _record(std::string filename, std::vector<double>& data);

double Interpol::F(double x) const {
	return log(x) + sqrt(1 + x);
}

template <typename T>
void print(T x) {
	std::cout << x << std::endl;
}

template <typename T>
void print_vector(std::vector<T> vec) {
	for( T x : vec ) {
		std::cout << x << " ";
	}
	std::cout << std::endl;
}

Interpol::Interpol(double l, double r, int _N, double s) {
	left = l;
	right = r;
	N = _N;
	step = s;
	for(double cur = l; cur <= right; cur += step) {
		X.push_back(cur);
		Y.push_back(F(cur));
	}
}

void upload(std::vector<double>& dest, std::vector<double>& source) {
	for(int i = 1; i < source.size(); ++i) {
		dest.push_back(source[i] - source[i - 1]);
	}
}

long factorial(int n) {
	if ( n == 0 ) {
		return 1;
	}
	return n * factorial(n - 1);
}

double stepen(double q, int n) {	///CHEck
	double reply = q;
	//std::cout << " q = " << q << std::endl;
	for(int i = 0; i < n - 1; ++i) {
		q--;
		reply *= q;
	}
	//std::cout << "stepen_reply = " << reply << std::endl;
	return reply;
}

void Interpol::Newton() {
	D.push_back(Y[0]);
	std::vector<double> d_cur;
	std::vector<double> d_prev;
	d_prev = Y;
	for(int i = 1; i < N; ++i) {
		upload(d_cur, d_prev);
		D.push_back(d_prev[1] - d_prev[0]);
		d_prev.clear();
		d_prev = d_cur;
		d_cur.clear();
	}
	record_newton();
}

void Interpol::Lagrange() {
	std::vector<double> x_lag;
	std::vector<double> y_lag;
	for(double cur = left; cur <= right; cur += 0.1) {
                x_lag.push_back(cur);
                y_lag.push_back(polynom_lagrange(cur));
        }
        _record("x_lag", x_lag);
        _record("y_lag", y_lag);

}

void Interpol::cube() {
	std::vector<double> x_cub;
	std::vector<double> y_cub;
	for(double cur = left; cur <= right; cur += step) {
                x_cub.push_back(cur);
                y_cub.push_back(F(cur));
        }
        _record("x_cub", x_cub);
        _record("y_cub", y_cub);
}


double Interpol::cub_spline(double x) {
	return 0;

}

double Interpol::l(double x, int num) {
	double reply = 1;
	for(int i = 0; i < N; ++i) {
		if ( i != num ) {
			reply *= ( x - X[i] ) / ( X[num] - X[i] );
		}
	}
	return reply;
}

double Interpol::polynom_newton(double x) {
	double q = ( x - left ) / step;
	double reply = D[0];
	for(int i = 1; i < N; ++i) {
		reply += D[i] * stepen(q, i) / factorial(i);
	}
	return reply;
}

double Interpol::polynom_lagrange(double x) {
	double reply = 0;
	for(int i = 0; i < N; ++i) {
		reply = ( Y[i] * l(x, i) ) + reply;
	}
	return reply;
}

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

void Interpol::record_newton() {
	std::vector<double> x_new;
	std::vector<double> y_new;
	//print(polynom_newton(3.5));
	for(double cur = left; cur <= right; cur += 0.1) {
		x_new.push_back(cur);
		y_new.push_back(polynom_newton(cur));
	}
	_record("x_newton", x_new);
	_record("y_newton", y_new);
}






