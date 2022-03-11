#pragma once
#include <iostream>
#include <cmath>
#include <vector>


class Functions {
	private:
		double accuracy;
		double F(double x) const;
	public:
		Functions(double e);
		void divide() const;
		void chord() const;
		void secants() const;
		void Newton() const;
		double search_tangent(double x) const;
};

Functions::Functions(double e) {
	accuracy = e;
}

double Functions::F(double x) const {
	return ( x - 1 ) * ( x - 1 ) * ( x - 1 );
}

void Functions::divide() const {
	double a_k = -10;
	double b_k = 10;
	int count = 0;
	double middle = ( a_k + b_k ) / 2;
	double F_left = F(a_k);
	double F_right = F(b_k);
	double F_middle = F(middle);
	while( fabs(F_middle) > accuracy ) {
		count++;
		if ( F_middle * F_right < 0 ) {
			a_k = middle;
			middle = ( a_k + b_k ) / 2;
			F_left = F(a_k);
			F_middle = F(middle);
		}
		else if ( F_middle * F_left < 0 ) {
			b_k = middle;
			middle = ( a_k + b_k ) / 2;
			F_middle = F(middle);
			F_right = F(b_k);
		}
	}
	std::cout << "solution_divide = " << middle << " count = " << count << std::endl;
}

double search_zero(double a, double b, double fa, double fb) {
	double k = ( fb - fa ) / ( b - a );
	double c = fa - k * a;
	return -1 * c / k;
} 

void Functions::chord() const {
	double a_k = -10;
	double b_k = 10;
	int count = 0;
	double x_sol = search_zero(a_k, b_k, F(a_k), F(b_k));
	while ( fabs(F(x_sol)) > accuracy ) {
		x_sol = search_zero(x_sol, b_k, F(x_sol), F(b_k));
		count++;
	}
	std::cout << "solution of chord = " << x_sol << " count = " << count << std::endl;
}

void Functions::secants() const {
	std::vector<double> x;
	x.push_back(-10);
	x.push_back(10);
	int count = 0;
	while ( fabs( x[count + 1] - x[count] ) > accuracy ) {
		double x_next = x[count] - ( F(x[count]) * ( x[count + 1] - x[count] ) / ( F(x[count + 1] ) - F(x[count]) ) );
		x.push_back(x_next);
		count++;
	}
	std::cout << "solution of secants = " << x[count + 1] << " count = " << count << std::endl;
}

double Functions::search_tangent(double x_prev) const {
	double deriavate = 3 * ( x_prev - 1 ) * ( x_prev - 1); 
	return x_prev - F(x_prev) / deriavate;
}


void Functions::Newton() const {
	double x_prev = 10;
	double x_next = this->search_tangent(x_prev);
	double delta = x_next - x_prev;
	int count = 0;
	while ( fabs(x_next - x_prev) > accuracy ) {
		count++;
		delta = x_next - x_prev;
		x_prev = x_next;
		x_next = this->search_tangent(x_next);
		//delta = x_next - x_prev;
	}
	double q = ( x_next - x_prev ) / delta;
	std::cout << "Metod of Newton's solution = " << x_next << " count = " << count << " S = " << round( 1 / ( 1 - q ))  << std::endl;
}


