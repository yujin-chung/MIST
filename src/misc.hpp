/*MIST copyright 2016 Yujin Chung and Jody Hey */


#ifndef _misc_
#define _misc_

#include <vector>
#include <valarray>
#include "mt.h"
#include <unistd.h>

extern MersenneTwister mt;

using namespace std;

class Matrix
{
	valarray<double> mat;
	int nrow, ncol;

public:
	//~Matrix();
	void initialize(int nr, int nc);
	void initialize(double val, int nr, int nc);
	void initialize_scalar(double val);
	double val(int row, int col);
	void replace(valarray<double> newMat,int nr, int nc);
	void replace(Matrix newMat);
	void print();
	vector<int> get_nodeIDs_minMat ();
	valarray<double> getcol(int col);
	valarray<double> getrow(int row);
	void replace(double val, int row, int col);
	void del_row(int row);
	void del_col(int col);
	void add_row(valarray<double> arr);
	void add_col(valarray<double> arr);
	void update_distMat(vector<int> nodeIDs, vector<int> sizes);
	void update_distMat_IS(vector<int> nodeIDs, vector<int> sizes);
	int get_nrow();
	int get_ncol();
	void replace_row(Matrix newRow, int nr);
	void replace_row(valarray<double> arr, int nr);
	valarray<double> get_valarray();
};

Matrix matrixMultiplication(Matrix mat1, Matrix mat2);
Matrix matrix_pwProduct(Matrix mat1, Matrix mat2);

/*
class seq_double
{
	valarray<double> val;
	int size;

public:
	void initial_zero(int val_len);


};
*/
int choose(int n, int k);

double runiform();
int runiform_discrete(int upperlimit);
double rexpdist (double rate);
double rNormalDistribution (double mean, double stdev);

double log_incompleteUpperGamma (int a, double x);

void print_vectorInt(vector<int> vec);

unsigned int factorial(unsigned int x);
double factorial(double x);
double logfactorial(double x);
int choose(int n, int k);


/**  YC 10/19/2016
 *   To check if a file exists   **/
inline bool file_exist(const std::string& name) {
  return ( access( name.c_str(), F_OK ) != -1 );
};

#endif


