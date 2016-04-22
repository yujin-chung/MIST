/*MIST copyright 2016 Yujin Chung and Jody Hey */




#ifndef _MCMC_H_
#define _MCMC_H_

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <ctime>
// FIXME
// Arun will need this header files for random number generator
// #include <gsl/gsl_rng.h>
// #include <gsl/gsl_randist.h>
#include <math.h>
#include "mt.h"
#include <stdio.h>
#include <cstring>
#include <fstream>
#include <vector>

int swap_attempt (int* swapper, int* swappee, int iter, double u2, double R, int p);
void print_coldchain(std::ofstream& fp, int chains, std::vector<double> beta, std::vector<double> posterior, int swapped, int iters);
void print_chain(std::ofstream& fp, int chains, std::vector<double> beta, std::vector<double> posterior, int iters);
void fill_swap_arrays(int* swapper, int* swappee, int processes, int max_iter, int chains);
// FIXME
// Arun will need this function for random number generator
// double calc_posterior(gsl_rng *r, std::vector<double> x, std::vector<double> beta, double alpha, int which);
double calc_R(double aposterior, double abeta, double bposterior, double bbeta);

#endif



