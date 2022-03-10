#pragma once

#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

class Matrix {
        private:
            int row;
            int col;
        public:
            double** arr;
            Matrix(int _row, int _col);
            Matrix(int _row, int _col, double** _arr);
            Matrix(Matrix& A);
            Matrix(Matrix& A, int a);
	    Matrix(int row, int col, int a);
	    int get_row() const;
            int get_col() const;
	    void print() const;
	    bool check_diag(double* x);
	    bool check_sym() const;
	    void gauss();
	    void gauss_main_elem();
	    bool shuttle();
	    void yacobi(double accuracy);
	    void zeydel(double accuracy);
	    void gradient_down(double accuracy);
	    void discrepancy(double accuracy);
	    void swap(int dest, int source);
            void mul(int num, double a);
	    void mul(double a);
	    double scalar_multiply(Matrix& right) const;
	    void sort();
	    double count_incon() const;
	    int search_min(int num) const;
	    int search_max(int num) const;
	    void sub(int a, int b, double k);
	    void record(std::string name) const;
            bool mul_impr(const Matrix& right, Matrix& result);
	    bool add_matrix(const Matrix& A);
	    bool sub_matrix(const Matrix& left, const Matrix& right);
	    void print_reply() const;
	    double count_norm() const;
	    ~Matrix();
};

Matrix::Matrix(int _row, int _col) {
        row = _row;
        col = _col;
        arr = new double*[row];
        for(int i = 0; i < row; ++i) {
            arr[i] = new double[col];
        }
}

Matrix::Matrix(int _row, int _col, int a) {
        row = _row;
        col = _col;
        arr = new double*[row];
        for(int i = 0; i < row; ++i) {
            arr[i] = new double[col];
        }
	for(int i = 0; i < row; ++i) {
		for(int j = 0; j < col; ++j) {
			arr[i][j] = a;
		}
	}
}



Matrix::Matrix(int _row, int _col, double** _arr) {
        row = _row;
        col = _col;
        arr = new double*[row];
        for(int i = 0; i < row; ++i) {
            arr[i] = new double[col];
        }
        for(int i = 0; i < row; ++i) {
            for(int j = 0; j < col; ++j) {
                arr[i][j] = _arr[i][j];
            }
        }
}

Matrix::Matrix(Matrix& A) {
	row = A.get_row();
	col = A.get_col();
        for(int i = 0; i < row; ++i) {
		arr[i] = new double[col];

        }
        for(int i = 0; i < row; ++i) {
            for(int j = 0; j < col; ++j) {
		    std::cout << "i = " << i << " j = " << j << std::endl;
		    arr[i][j] = A.arr[i][j];
            }
        }
}

int Matrix::get_row() const {
        return row;
}

int Matrix::get_col() const {
        return col;
}


Matrix::~Matrix() {
        for(int i = 0; i < row; ++i) {
            delete[] arr[i];
        }
        delete[] arr;
}

void Matrix::print() const {
        for(int i = 0; i < row; ++i) {
            for(int j = 0; j < col; ++j) {
                std::cout << arr[i][j] << " ";
            }
            std::cout << std::endl;
        }
}

void Matrix::record(std::string name) const {
        std::ofstream fout(name);
        for(int i = 0; i < row; ++i) {
            for(int j = 0; j < col; ++j) {
                fout << arr[i][j] << " ";
            }
            fout << std::endl;
        }
        fout.close();
}

double Matrix::scalar_multiply(Matrix& right) const {
	if ( row != right.get_row() ) {
		std::cout << "Error in size\n";
		return -1;
	}
	double sum = 0;
	for(int i = 0; i < row; ++i) {
		sum += arr[i][0] * right.arr[i][0];
	}
	return sum;
}



bool Matrix::mul_impr(const Matrix& right, Matrix& result) {              // Matrix C = A * B;
        if ( col != right.row ) {
            std::cout << "WRONG! Incorrect matrix dimensions are set\n";
            return false;
        }
        int M = result.row;
        int N = result.col;
        int num = 0;
        double c = 0;
        for(int j = 0; j < N; ++j) {
            double Y[col];           //столбец
            for(int q = 0; q < col; ++q) {
                Y[q] = right.arr[q][num];
            }
            num++;
            for(int i = 0; i < M; ++i) {            //ускорение перестановкой циклов, вынос переменных
                c = 0;
                int k = 7;
                while ( k < col ) {                 //loop unrolling
                    c += arr[i][k-7] * Y[k-7];
                    c += arr[i][k-6] * Y[k-6];
                    c += arr[i][k-5] * Y[k-5];
                    c += arr[i][k-4] * Y[k-4];
                    c += arr[i][k-3] * Y[k-3];
                    c += arr[i][k-2] * Y[k-2];
                    c += arr[i][k-1] * Y[k-1];
                    c += arr[i][k] * Y[k];
                    k += 8;
                }
                k -= 8;
                if ( col < 8 ) {                                              //ch
                    for(int ik = 0; ik < col; ++ik) {
                        c += arr[i][ik] * Y[ik];
                    }
                }
                else if ( col % 8 != 0 ) {
                    for(int ik = k + 1; ik < col; ++ik) {
                        c += arr[i][ik] * Y[ik];
                    }
                }
                result.arr[i][j] = c;
            }
        }
        //result.print();
        return true;
        
    }

bool Matrix::add_matrix(const Matrix& A) {
        if ( ( row != A.get_row() ) & ( col != A.get_col() ) ) {
            return false;
        }
        for(int i = 0; i < row; ++i) {
            for(int j = 0; j < col; ++j) {
                arr[i][j] += A.arr[i][j];
            }
        }
        return true;
}

bool Matrix::sub_matrix(const Matrix& right, const Matrix& left) {
	if ( ( right.get_row() != left.get_row() ) | ( right.get_col() != left.get_col() )  ) {
			return false;
	}
	for(int i = 0; i < row; ++i) {
		for(int j = 0; j < col; ++j) {
			arr[i][j] = right.arr[i][j] - left.arr[i][j];
		}
	}
	return true;
}

void Matrix::mul(double a) {
	for(int i = 0; i < row; ++i) {
		for(int j = 0; j < col; ++j) {
			arr[i][j] *= a;
		}
	}
}

void Matrix::mul(int num, double a) { 		//col
	for(int i = 0; i < col; ++i) {
		arr[num][i] *= a;
	}
}

void Matrix::sub(int a, int b, double k) {
	if ( a == b ) {
		return;
	}
	for(int i = 0; i < col; ++i) {
		arr[a][i] = arr[a][i] - k * arr[b][i];
	}
}

void Matrix::swap(int dest, int source) {
	double temp;
	//std::cout << "SWAP\n";
	for(int i = 0; i < col; ++i) {
		temp = arr[dest][i];
		arr[dest][i] = arr[source][i];
		arr[source][i] = temp;
	}
}

int Matrix::search_min(int num) const {
	for(int i = num; i < row; ++i) {
		if ( arr[i][num] != 0 ) {
			return i;
		}
	}
	return -1;
}

void Matrix::sort() {
	for(int i = 0; i < row; ++i) {
		if ( arr[i][i] == 0 ) {
			//std::cout << "number sort = " << i << std::endl;
			int next = this->search_min(i);
			if ( next == -1 ) {
				return;
			}
			swap(i, next);
		}
	}
}

void Matrix::gauss() {
	int pointer = 0;
	//this->print();
	for(int j = 0; j < row; ++j) {
		this->sort();
		for(int i = 0; i < row; ++i) {
			if ( arr[j][pointer] != 0 ) {
				this->sub(i, j, arr[i][pointer] / arr[j][pointer]);
			}
		}
		if ( arr[j][pointer] != 0 ) {
                	this->mul(pointer, 1 / arr[j][pointer]);
		}
		pointer++;
		//std::cout << "check:\n";
		//this->print();
	}
}

bool Matrix::check_sym() const {
	for(int i = 0; i < row; ++i) {
		for(int j = 0; j < col; ++j) {
			if ( arr[i][j] != arr[j][i] ) {
				return false;
			}
		}
	}
	return true;
}

double diff_norm(Matrix& a, Matrix& b) {
	double sum = 0;
	for(int i = 0; i < a.get_row(); ++i) {
		sum += ( ( a.arr[i][0] - b.arr[i][0] ) * ( a.arr[i][0] - b.arr[i][0] ) );
	}
	return sqrt(sum);
}

void yacobi_recursive(Matrix& B, Matrix& X_prev, double accuracy, double q, Matrix& g) {
	Matrix X(B.get_row(), 1);
	B.mul_impr(X_prev, X);
	X.add_matrix(g);
	//if ( diff_norm(X, X_prev) < ( ( 1 - q ) / q * accuracy ) ) {
	if ( diff_norm(X, X_prev) < accuracy ) {
		X.print();
		return;
	}
	else {
		yacobi_recursive(B, X, accuracy, q, g);
	}
}

void Matrix::yacobi(double accuracy) {			//add check_norm
	Matrix B(row, col - 1, 0.0);
	Matrix x0(row, 1);
	for(int i = 0; i < row; ++i) {
		this->sort();
		double k = -1 / arr[i][i];
		for(int j = 0; j < col - 1; ++j) {
			B.arr[i][j] = arr[i][j] * k;
		}
		B.arr[i][i] = 0;
		x0.arr[i][0] = arr[i][col - 1] / arr[i][i];

	}
	double norm = B.count_norm();
	if ( norm < 1 ) {
		yacobi_recursive(B, x0, accuracy, norm, x0);
	}
	else {
		std::cout << "The display is not compressive" << std::endl;
	}
}


void zeydel_recursive(Matrix& B, Matrix& X_prev, double accuracy, double q, Matrix& g) {
        Matrix X(B.get_row(), 1);
	for(int i = 0; i < B.get_row(); ++i) {
		for(int j = 0; j < B.get_row(); ++j) {
			if ( j >= i ) {
				X.arr[i][0] = X_prev.arr[j][0] * B.arr[i][j] + X.arr[i][0];
			}
			else {
				X.arr[i][0] = X.arr[j][0] * B.arr[i][j] + X.arr[i][0];
			}
		}
		X.arr[i][0] += g.arr[i][0];
	}
        if ( diff_norm(X, X_prev) < accuracy ) {
                X.print();
                return;
        }
        else {
                zeydel_recursive(B, X, accuracy, q, g);
        }
}


void Matrix::zeydel(double accuracy) {                  //add check_norm
        Matrix B(row, col - 1, 0.0);
        Matrix x0(row, 1);
        for(int i = 0; i < row; ++i) {
                this->sort();
                double k = -1 / arr[i][i];
                for(int j = 0; j < col - 1; ++j) {
                        B.arr[i][j] = arr[i][j] * k;
                }
                B.arr[i][i] = 0;
                x0.arr[i][0] = arr[i][col - 1] / arr[i][i];

        }
        double norm = B.count_norm();
        if ( norm < 1 ) {
                zeydel_recursive(B, x0, accuracy, norm, x0);
        }
        else {
                std::cout << "The display is not compressive" << std::endl;
        }
}

//int nar = 0;

void gradient_down_recursive(Matrix& xk, Matrix& A, Matrix& b, double accuracy, int row) {
	Matrix pk(row, 1);
	A.mul_impr(xk, pk);
	Matrix rk(row, 1);
	rk.sub_matrix(b, pk);
	Matrix t(row, 1, 0);
	Matrix tk(row, 1);
	if ( diff_norm(rk, t) < accuracy ) {
		xk.print();
		return;
	}
	else {
		A.mul_impr(rk, tk);
		double ak = ( rk.scalar_multiply(rk) ) / ( rk.scalar_multiply(tk) );
		rk.mul(ak);
		xk.add_matrix(rk);
		gradient_down_recursive(xk, A, b, accuracy, row);
	}
}

void Matrix::gradient_down(double accuracy) {
	Matrix x0(row, 1, 1);
	Matrix A(row, row);
	Matrix b(row, 1);
	for(int i = 0; i < row; ++i) {
		for(int j = 0; j < col; ++j) {
			if ( j == row ) {
				b.arr[i][0] = arr[i][j];
			}
			else {
				A.arr[i][j] = arr[i][j];
			}
		}
	}
	bool fl = A.check_sym();
	if ( fl == false ) {
		std::cout << "Matrix isn't symmetry" << std::endl;
		return;
	}
	gradient_down_recursive(x0, A, b, accuracy, row);
}

void discrepancy_recursive(Matrix& x_k, Matrix& A, Matrix& b, int row, double accuracy) {
	Matrix Ax_k(row, 1);
	Matrix Y_k(row, 1);
	A.mul_impr(x_k, Ax_k);
	Y_k.sub_matrix(Ax_k, b);
	Matrix Ay_k(row, 1);
	A.mul_impr(Y_k, Ay_k);
	if ( Ay_k.scalar_multiply(Ay_k) == 0 ) {
		x_k.print();
		return;
	}
	double t_k = (Y_k.scalar_multiply(Ay_k)) / (Ay_k.scalar_multiply(Ay_k));
	Matrix x_next(row, 1);
	Y_k.mul(t_k);
	x_next.sub_matrix(x_k, Y_k);
	if ( diff_norm(x_next, x_k) < accuracy ) {
		x_next.print();
	}
	else {
		discrepancy_recursive(x_next, A, b, row, accuracy);
	}

}

void Matrix::discrepancy(double accuracy) {
	Matrix x0(row, 1, 0);
        Matrix A(row, row);
        Matrix b(row, 1);
        for(int i = 0; i < row; ++i) {
                for(int j = 0; j < col; ++j) {
                        if ( j == row ) {
                                b.arr[i][0] = arr[i][j];
                        }
                        else {
                                A.arr[i][j] = arr[i][j];
                        }
                }
        }
	bool fl = A.check_sym();
	if ( fl == false ) {
		std::cout << "Matrix isn't symmetry" << std::endl;
		return;
	}
	discrepancy_recursive(x0, A, b, row, accuracy);
}

double modul(double a) {
	if ( a >= 0 ) {
		return a;
	}
	return -1 * a;
}

double Matrix::count_incon() const {
	Matrix A(row, row);
	Matrix B(row, 2 * row);
	for(int i = 0; i < row; ++i) {
		for(int j = 0; j < row; ++j) {
			A.arr[i][j] = arr[i][j];
		}
	}
	for(int i = 0; i < row; ++i) {
		for(int j = 0; j < 2 * row; ++j) {
			if ( j >= row ) {
				if ( i == j - row ) {
					B.arr[i][j] = 1;
				}
				else {
					B.arr[i][j] = 0;
				}
			}
			else {
				B.arr[i][j] = arr[i][j];
			}
		}
	}
	B.gauss();
	Matrix C(row, row);
	Matrix D(row, row);
	for(int i = 0; i < row; ++i) {
		for(int j = 0; j < row; ++j) {
			if ( ( i != j ) & ( B.arr[i][j] > 0.00001 ) ) {
				std::cout << "Matrix inconsistent\n";
				return 0;
			}
			C.arr[i][j] = B.arr[i][row+j];
		}
	}
	return A.count_norm() * C.count_norm();
}


double Matrix::count_norm() const {
	double max = 0;
	for(int i = 0; i < row; ++i) {
		double sum = 0;
		for(int j = 0; j < row; ++j) {
			sum += modul(arr[i][j]);
		}
		if ( sum > max ) {
			max = sum;
		}
	}
	return max;
}


int Matrix::search_max(int num) const {
	double max = modul(arr[num][num]);
	int reply = num;
	for(int i = num; i < row; ++i) {
                if ( modul(arr[i][num]) > max ) {
                        max = modul(arr[i][num]);
			reply = i;
                }
        }
        if ( max == 0 ) {
		return -1;
	}
	return reply;
}


void Matrix::gauss_main_elem() {
	int pointer = 0;
        for(int j = 0; j < row; ++j) {
                int number = this->search_max(j);
		if ( number == -1 ) {
			return;
		}
		else {
			swap(j, number);
		}	
		for(int i = 0; i < row; ++i) {
                        if ( arr[j][pointer] != 0 ) {
                                this->sub(i, j, arr[i][pointer] / arr[j][pointer]);
                        }
                }
                if ( arr[j][pointer] != 0 ) {
                        this->mul(pointer, 1 / arr[j][pointer]);
                }
                pointer++;
        }
}

bool Matrix::check_diag(double* x) {
	for(int i = 0; i < row; ++i) {
		double diag = arr[i][i];
		double sum = 0;
		for(int j = 0; j < col - 1; ++j) {
			sum += modul(arr[i][j]);
		}
		if ( sum > 2 * diag ) {
			return false;
		}
	}
	double m;
	for(int i = 1; i < row; ++i) {
		m = arr[i][i-1] / arr[i-1][i-1];
	
		arr[i][i] = arr[i][i] - m * arr[i-1][i];
		arr[i][row] = arr[i][row] - m * arr[i-1][row];
		//std::cout << "arr[i][i] = " << arr[i][i] << " arr[i][row] = " << arr[i][row] << " m = " << m << std::endl;
	}
	x[row-1] = arr[row-1][row] / arr[row-1][row-1];
	for(int i = row - 2; i >= 0; i--) {
		x[i] = (arr[i][row] - arr[i][i+1] * x[i+1]) / arr[i][i];
	}
	return true;
}


bool Matrix::shuttle() {
	std::cout << "SHUTTLE\n";
	double x[row];
	bool flag = this->check_diag(x);
	if ( flag == false ) {
		return false;
	}
	for(int i = 0; i < row; ++i) {
		std::cout << x[i] << std::endl;
	}
	return true;
}


void Matrix::print_reply() const {
	int count = 0;
	double EPSILON = 0.000001;
	for(int i = 0; i < row; ++i) {
		for(int j = 0; j < col; ++j) {
			if ( modul(arr[i][j]) > EPSILON ) {
				count++;
			}
		}
		if ( count == 2 ) {
			std::cout << arr[i][col - 1] << std::endl;
			count = 0;
		}
		else if ( count == 0 ) { 
			std::cout << "Anything\n";
			count = 0;
		}
		else {
                        std::cout << "Matrix inconsistent\n";
                        return;
                }
	}
}



double** init(double** arr, int m, int n) {
    arr = new double*[m];
    for(int i = 0; i < m; ++i) {
        arr[i] = new double[n];
    }
    return arr;
}
