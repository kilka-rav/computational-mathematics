#pragma once


#include <iostream>
#include <cmath>
#include <vector>

class Functions {
	private:
		double left;
		double right;
		double step;
		std::vector<double> X;
		std::vector<double> Y;
	public:
		Functions(double l, double r, double s);
		Functions(double l, double s, double accuracy, bool flag);
		double F(double x) const;
		double F2(double x) const;
		double rectangle() const;
		double rectangle(double accuracy);
		double trapezoid() const;
		double trapezoid(double accuracy);
		double sympson() const;
		double sympson(double accuracy);
		double runge_modif;
		double unown();
};

Functions::Functions(double l, double r, double s) {
	left = l;
	right = r;
	step = s;
	for(double cur = left; cur < right; cur += step) {
		X.push_back(cur);
		Y.push_back(F(cur));
	}
}

Functions::Functions(double l, double s, double accuracy, bool flag) {
	if ( flag == false ) {
		std::cout << "DANGEROUS" << std::endl;
		exit(1);
	}
	left = l;
	step = s;
	right = sqrt(fabs(log(accuracy)));
	//std::cout << "right = " << right << std::endl;
	for(double cur = left; cur < right; cur += step) {
		X.push_back(cur);
		Y.push_back(F2(cur));
	}
}

double Functions::F(double x) const {
	return 2 * sqrt(1 - x * x);
}

double Functions::F2(double x) const {
	return 2 * exp(-1 * x * x);
}

double Functions::unown() {
	double reply = 0;
	double accuracy = 0.01;
	double was = step / 2;
	double reply_min = this->sympson();
	//std::cout << std::endl;
	while(1) {
		Functions three(left, was, accuracy, true);
		reply = three.sympson();
		//std::cout << "rep = " << reply << std::endl;
		if ( fabs(reply - reply_min) < accuracy ) {
                        break;
                }
		was /= 2;
                reply_min = reply;
        }
        return reply;
}

double Functions::rectangle() const {
	double reply = 0;
	double width = X[1] - X[0];
	for( double y : Y ) {
		reply += y * width;
	}
	return reply;
}

double Functions::rectangle(double accuracy) {
	double reply = 0;
	double reply_min = this->rectangle();
	while(1) {
		Functions in(left, right, step / 2);
		reply = in.rectangle();
		if ( fabs(reply - reply_min) < accuracy ) {
			runge_modif = fabs(reply - reply_min);
			//std::cout << "\nRunge modification rect = " << fabs(reply - reply_min) / 3 << std::endl;
			break;
		}
		reply_min = reply;
	}
	return reply;
}

double Functions::trapezoid() const {
	double reply = 0;
	double width = X[1] - X[0];
	for(int i = 0; i < Y.size() - 1; ++i) {
		reply += ( Y[i] + Y[i + 1] ) * width / 2;
	}
	return reply;
}

double Functions::trapezoid(double accuracy) {
	double reply = 0;
	double reply_min = this->trapezoid();
        while(1) {
                Functions in(left, right, step / 2); 
                reply = in.trapezoid();        
                if ( fabs(reply - reply_min) < accuracy ) {
			runge_modif = fabs(reply - reply_min);
			//std::cout << "\nRunge modification trap = " << fabs(reply - reply_min) / 3 << std::endl;
                        break;
                }
		reply_min = reply;
        }       
        return reply;
}

double Functions::sympson() const {
	double reply = 0;
	double width = X[1] - X[0];
	int n = X.size() / 2;
	for(int i = 0; i < n - 1; ++i) {
		reply += width / 3 * ( Y[2 * i] + 4 * Y[2 * i + 1] + Y[2 * i + 2] );
	}
	return reply;
}

double Functions::sympson(double accuracy) {
	double reply = 0;
	double reply_min = this->sympson();
        while(1) {
                Functions in(left, right, step / 2); 
                reply = in.sympson();        
                if ( ( fabs(reply - reply_min) / 15 ) < accuracy ) {
			runge_modif = fabs(reply - reply_min);
			//std::cout << "\nRunge modification Symp = " << fabs(reply - reply_min) / 15 << std::endl;
			break;
                }
		reply_min = reply;
        }       
        return reply;
}

