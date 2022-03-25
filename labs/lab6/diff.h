#pragma once

#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

class Function {
	private:
		double y0;
		double x0;
		std::vector<double> X;
		std::vector<double> Y;
		double step;
		double eps;
	public:
		Function(double _x0, double _y0, double s);
		double F(double x, double y) const;
		double solution(double x) const;
		void record(std::string name);
		void explicit_euler(bool flag);
		void implicit_euler(bool flag);
		void recalc_euler(bool flag);
		void runge_kutta(bool flag);
		double get_eps() const;
};

double Function::F(double x, double y) const {
	return 2 * x * y;
}

double Function::solution(double x) const {
	return exp(x * x);
}

Function::Function(double _x0, double _y0, double s) {
	x0 = _x0;
	y0 = _y0;
	step = s;
	eps = 0;
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

void Function::record(std::string name) {
        _record("x_" + name, X);
        _record("y_" + name, Y);
}



void Function::explicit_euler(bool flag) {
	Y.push_back(y0);
	X.push_back(x0);
	double cur = x0 + step;
	for(int i = 1; cur < 1; i++) {
		X.push_back(cur);
		Y.push_back(Y[i - 1] + step * F(X[i - 1], Y[i - 1]));
		cur += step;
		if ( fabs(solution(X[i]) - Y[i]) > eps ) {
                        eps = fabs(solution(X[i]) - Y[i]);
                }

	}
	if ( flag == true ) {
		record("explicit_euler");
	}
}

double newton_mini(double x_next, double y_prev, double h) {
	return (y_prev / ( 1 - 2 * h * x_next));
}

void Function::implicit_euler(bool flag) {
	X.push_back(x0);
	Y.push_back(y0);
	double cur = x0 + step;
	for(int i = 1; cur < 1; ++i) {
		X.push_back(cur);
		double res = newton_mini(X[i], Y[i - 1], step);
		Y.push_back(res);
		cur += step;
		if ( fabs(solution(X[i]) - Y[i]) > eps ) {
			eps = fabs(solution(X[i]) - Y[i]);
		}
	}
	if ( flag == true ) {
		record("implicit_euler");
	}
}

void Function::recalc_euler(bool flag) {
	double y_mini;
	X.push_back(x0);
	Y.push_back(y0);
	y_mini = y0;
	double cur = x0 + step;
	for(int i = 1; cur < 1; ++i) {
		y_mini = Y[i - 1] + step * F(X[i - 1], Y[i - 1]);
		X.push_back(cur);
		Y.push_back(Y[i - 1] + step / 2 * (F(X[i - 1], Y[i - 1]) + F(X[i], y_mini)));
		cur += step;
		if ( fabs(solution(X[i]) - Y[i]) > eps ) {
                        eps = fabs(solution(X[i]) - Y[i]);
                }

	}
	if ( flag == true ) {
		record("recalc_euler");
	}
}

void Function::runge_kutta(bool flag) {
	double k1, k2, k3, k4;
	X.push_back(x0);
	Y.push_back(y0);
	double cur = x0 + step;
	for(int i = 1; cur < 1; ++i) {
		X.push_back(cur);
		k1 = F(X[i - 1], Y[i - 1]);
		k2 = F(X[i - 1] + step / 2, Y[i - 1] + step / 2 * k1);
		k3 = F(X[i - 1] + step / 2, Y[i - 1] + step / 2 * k2);
		k4 = F(X[i - 1] + step, Y[i - 1] + step * k3);
		Y.push_back(Y[i - 1] + step * (k1 + 2 * k2 + 2 * k3 + k4) / 6);
		cur += step;
		if ( fabs(solution(X[i]) - Y[i]) > eps ) {
                        eps = fabs(solution(X[i]) - Y[i]);
                }

	}
	if ( flag == true ) {
		record("runge_kutta");
	}
}

double Function::get_eps() const {
	return eps;
}

void count_error(double x0, double y0) {
	double h0 = 0.01;
	std::vector<double> h;
	std::vector<double> ex_err;
	std::vector<double> im_err;
	std::vector<double> reca_err;
	std::vector<double> runge_err;
	for(int i = 0; i < 3; ++i) {
		h.push_back(h0);
		h0 /= 2;
		Function ex1(x0, y0, h0);
		Function ex2(x0, y0, h0);
		Function ex3(x0, y0, h0);
		Function ex4(x0, y0, h0);
		ex1.explicit_euler(false);
		ex2.implicit_euler(false);
		ex3.recalc_euler(false);
		ex4.runge_kutta(false);
		ex_err.push_back(ex1.get_eps());
		im_err.push_back(ex2.get_eps());
		reca_err.push_back(ex3.get_eps());
		runge_err.push_back(ex4.get_eps());
	}
	_record("step", h);
	_record("ex_err", ex_err);
	_record("im_err", im_err);
	_record("recalc_err", reca_err);
	_record("runge_err", runge_err);
}
