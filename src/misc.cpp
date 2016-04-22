/* MIST copyright 2016 by Yujin Chung and Jody Hey */


#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "misc.hpp"
#include <time.h>

using namespace std;


void print_vectorInt(vector<int> vec)
{
	for(unsigned int i=0; i<vec.size(); i++)
		std::cout << vec.at(i) <<" ";
	std::cout << "\n";
}

int choose(int n, int k)
{
	int val = 0;
	if(k == n || k == 0)
		val = 1;
	else if(k <= n-k)
		val = n * choose(n-1,k-1)/k;
	else
		val = choose(n, n-k);
	return val;
}

// Copy from IMa2 code
/* (C) Copr. 1986-92 Numerical Recipes Software '$&'3$. */
double expint (int n, double x, int *islog)
{
	// YC 2014/2/20
	// Add the part of macro expansion here.
	int MAXIT = 100;
	double EULER  =0.5772156649;
	double FPMIN = 1.0e-30;
	double EPS= 3.0e-7;
	// End of the addition.

  int i, ii, nm1;
  double a, b, c, d, del, fact, h, psi, ans;
  *islog = 0;
  nm1 = n - 1;
  if (n < 0 || x < 0.0 || (x == 0.0 && (n == 0 || n == 1)))
  {
    std::cout << " expint() called from uppergamma() with bad value(s). n="
    		<< n << " x="<< x <<"\n";
  }
  else
  {
    if (n == 0)
    {
      ans = exp (-x) / x;
    }
    else
    {
      if (x == 0.0)
      {
        ans = 1.0 / nm1;
      }
      else
      {
        if (x > 1.0)
        {
          b = x + n;
          c = 1.0 /  FPMIN;
          d = 1.0 / b;
          h = d;
          for (i = 1; i <=MAXIT ; i++)
          {
            a = -i * (nm1 + i);
            b += 2.0;
            d = 1.0 / (a * d + b);
            c = b + a / c;
            del = c * d;
            h *= del;
            if (abs (del - 1.0) < EPS)
            {
              *islog = 1;

              //ans = h*exp(-x);
              ans = log (h) - x;
              return ans;
            }
          }

          // std::cout << " too many iterations in expint() called from uppergamma()";
        }
        else
        {
          ans = (nm1 != 0 ? 1.0 / nm1 : -log (x) - EULER);
          fact = 1.0;
          for (i = 1; i <= MAXIT; i++)
          {
            fact *= -x / i;
            if (i != nm1)
            {
              del = -fact / (i - nm1);
            }
            else
            {
              psi = -EULER;
              for (ii = 1; ii <= nm1; ii++)
                psi += 1.0 / ii;
              del = fact * (-log (x) + psi);
            }
            ans += del;
            if (abs (del) < abs (ans) * EPS)
              return ans;
          }
          // std::cout << " too many iterations in expint() called from uppergamma()";
        }
      }
    }
  }
  return ans;
}                               //expint

// Copy from IMa2 code
/* (C) Copr. 1986-92 Numerical Recipes Software '$&'3$. */
void gser (double *gamser, int a, double x, double *gln)
{
	unsigned int ITMAX=100000;
	double EPS=3.0e-7;

  unsigned int n;
  double sum, del, ap;

  // YC 2/20/2014
  // logfact[a-1] is log((a-1)!) = log(1)+log(2)+...+log(a-1).
  // I replaced the next statement by the following codes.
  // *gln = logfact[a - 1];
  double gln_tmp = 0.0;
  for(int i=1;i<= a-1; i++)
	  gln_tmp += log((double) i);
  *gln = gln_tmp;
  // End of edition.

  //FIXME
  // cout << "a = "<< a <<"; gln_tmp = " << gln_tmp <<"\n";

  if (x <= 0.0)
  {
    if (x < 0.0)
    	std::cout << " gser() called with negative x value, x: "<< x <<"\n";
     // IM_err (IMERR_UPPERGAMMA,
     //         " gser() called with negative x value, x: %lf ", x);
    *gamser = 0.0;
    return;
  }
  else
  {
    ap = a;
    del = sum = 1.0 / a;
    for (n = 1; n <= ITMAX; n++)
    {
      ++ap;
      del *= x / ap;
      sum += del;


      if (abs(del) < abs(sum) * EPS)
      {
    	//  *gamser = sum * exp (-x + a * log (x));
        //*gamser = sum * exp (-x + a * log (x) - (*gln));
        //return;
    	  break;
      }
      else if(n == ITMAX)
      	std::cout << " too many iterations within gser() \n";
    }


    // IM_err (IMERR_UPPERGAMMA, " too many iterations within gser() ");
    *gamser = sum * exp (-x + a * log (x));
    //*gamser = sum * exp (-x + a * log (x) - (*gln));
    return;
  }
}

// Copy of IMa2
// gcflog is the same as gcf but returns the logarithm of gammcf
// saves a little time and stops some underflows
void gcflog (double *gammcflog, double a, double x, double *gln)
{
	// Copy of the macro expansion in IMa2
	double FPMIN= 1.0e-30;
	int ITMAX= 10000;
	double EPS =3.0e-7;
	// End of addition

	int i;
	double an, b, c, d, del, h;


	// YC 2/20/2014
    // logfact[a-1] is log((a-1)!) = log(1)+log(2)+...+log(a-1).
    // I replaced the next statement by the following codes.
	// --- original code -------------//
    // *gln = logfact[(int) a - 1];
	// --- new code -------------------//
    double gln_tmp = 0.0;
    for(int i=1;i<= a-1; i++)
    	gln_tmp += log((double) i);
    *gln = gln_tmp;
    // --- end of the new code ---------//
    // More note about 'gln'.
    // 'gln' is log of the gamma function, Gamma(a).
    // If "a" is a non-negative integer, then
    // Gamma(a) = (a-1)!. Thus, 'gln' is
    // log(1) + log(2) + ... + log(a-1).


    // REMOVE
    //cout << "a = "<< a << "; gln_tmp = " << gln_tmp <<"\n";

    b = x + 1.0 - a;
    c = 1.0 / FPMIN;
    d = 1.0 / b;
    h = d;
    for (i = 1; i <= ITMAX; i++)
    {
    	an = -i * (i - a);
    	b += 2.0;
    	d = an * d + b;
    	if (abs (d) < FPMIN)
    		d = FPMIN;
    	c = b + an / c;
    	if (abs (c) < FPMIN)
    		c = FPMIN;
    	d = 1.0 / d;
    	del = d * c;
    	h *= del;
    	if (abs (del - 1.0) < EPS)
    		break;
    	else if (i == ITMAX)
    	  std::cout << " too many iterations within gcflog()\n";
    }
    // IM_err (IMERR_UPPERGAMMA, " too many iterations within gcflog()");

    //*gammcf=exp(-x+a*log(x)-(*gln))*h;
    //*gammcflog = (-x + a * log (x) - (*gln)) + log (h);
    *gammcflog = (-x + a * log (x) ) + log (h);

}                               //gcflog

/**
 *	Incomplete upper gamma function,
 *	    \int_x^{\infty} t^{a-1}e^{-t} dt,
 *	is returned in log scale.
 *	@param a scale of the gamma distribution (rate = 1)
 *	@param x a numeric value
 */
double log_incompleteUpperGamma (int a, double x)
/*  Returns the log of what Mathematica calls the incomplete gamma function.   */
{

	// Remove
	// std::cout << "in log_incompleteUpperGamma()\n";
	// std::cout << "\t a= "<<a << "; x= "<< x<< "\n";

  int logindicator;
  double gamser, gammcf, gln, p;
  double temp;
  if (x < 0.0 || a < 0)
	  std::cout << " uppergamma arguments: a="<< a<< " x="<<x <<"\n";
    // IM_err (IMERR_UPPERGAMMA, "  step %d, uppergamma arguments: a  %d, x %lf",step, a, x);
  else if (a == 0) // the case of one coalescent event (2 lineages)
  {
	  // Exponential integration Ei in log scale
	  // Ei(x) = \int_x^{\infty} exp(-t)/t dt
    temp = expint (1, x, &logindicator);
    if (!logindicator)
      p = log (temp);

    else
      p = temp;
  }
  else
  {
    if (x < (a + 1.0))
    {
      gser (&gamser, a, x, &gln);
      //p = gln + log (1.0 - gamser);
      // p = log (1.0 - gamser);
      p = log (exp(gln) - gamser);
      //REMOVE
      // cout << "gamser = " << gamser << "\n";
    }
    else
    {
      //gcf(&gammcf,a,x,&gln);
      //p = gln + log(gammcf);
      gcflog (&gammcf, (double) a, x, &gln);
      // p = gln + gammcf;
      p = gammcf;
      // YC 2/21/2014
      // 'gln' is the Gamma function, Gamma(a), in log scale
      // and 'gammcf' is the Gamma distribution function, gamma(x,a), in log scale.

      //REMOVE
      // cout << "gammacf = " << gammcf << "; exp(gammacf) = "<< exp(gammcf)<<"\n";
    }
  }
  if (p < -1e200)
    p = -1e200;

  // REMOVE
  // cout << "p = " << p <<"\n";

  return p;
}                               //uppergamma





void Matrix::initialize(int nr, int nc)
{
	nrow = nr;
	ncol = nc;
	mat.resize(nrow * ncol);
}


void Matrix::initialize_scalar(double val)
{
	nrow = 1;
	ncol = 1;
	valarray<double> tmp(1);

}


void Matrix::initialize(double val, int nr, int nc)
{
	nrow = nr;
	ncol = nc;
	mat.resize((unsigned long) nrow * ncol);
	for(unsigned int i=0; i< (unsigned) nrow*ncol; i++)
		mat[i] = val;
	//valarray<double> tmp (val ,nrow*ncol);
	//mat = tmp;
}


void Matrix::replace(Matrix newMat)
{
	nrow = newMat.nrow;
	ncol = newMat.ncol;
	mat.resize(nrow*ncol);
	mat = newMat.mat;
}

void Matrix::add_row(valarray<double> arr)
{
	valarray<double> new_mat((nrow+1)*ncol);
	int i;
	for(i=0; i<nrow; i++)
		new_mat[slice(i * ncol, ncol,1)] = Matrix::getrow(i);
	new_mat[slice(nrow * ncol, ncol,1)] = arr;
	nrow++;
	mat.resize(nrow*ncol);
	mat = new_mat;
}
void Matrix::add_col(valarray<double> arr)
{
	valarray<double> new_mat(nrow*(ncol+1));
	int i;
	for(i=0; i<ncol; i++)
		new_mat[slice(i, nrow,ncol+1)] = Matrix::getcol(i);
	new_mat[slice(ncol, nrow,ncol+1)] = arr;
	ncol++;
	mat.resize(nrow*ncol);
	mat = new_mat;
}

void Matrix::del_row(int row)
{
	int i=0;
	nrow--;
	valarray<double> new_mat(nrow*ncol);
	for(i=0; i<nrow; i++)
		if(i<row)
			new_mat[slice(i * ncol, ncol,1)] = Matrix::getrow(i);
		else
			new_mat[slice(i * ncol, ncol,1)] = Matrix::getrow(i+1);
	mat.resize(nrow*ncol);
	mat = new_mat;
}

void Matrix::del_col(int col)
{
	int i=0;
	valarray<double> new_mat(nrow*(ncol-1));
	for(i=0; i<ncol-1; i++)
		if(i<col)
			new_mat[slice(i, nrow,ncol-1)] = Matrix::getcol(i);
		else
			new_mat[slice(i, nrow, ncol-1)] = Matrix::getcol(i+1);
	ncol--;
	mat.resize(nrow*ncol);
	mat = new_mat;
}

double Matrix::val(int row, int col)
{
	std::size_t loc = row * ncol +col;
	return mat[loc];
}

void Matrix::replace(double val, int row, int col)
{
	int loc = row * ncol +col;
	mat[loc] = val;
}

valarray<double> Matrix::getrow(int row)
{
	return mat[slice(row * ncol, ncol,1)];
}

valarray<double> Matrix::getcol(int col)
{
	valarray<double> newCol;
	newCol.resize(nrow);
	try{
		newCol = mat[slice(col, nrow, ncol)];
	}catch (std::exception &e) {
		std::cout << "In Matrix::getcol(). Can't access elements of Matrix"
				" - index out of bounds\n";
	}

	// REMOVE
	// std::cout << "The current matrix is ";
	// print();
	// std::cout << col <<"th column is ";
	// for(unsigned int i=0; i<newCol.size(); i++)
	// std::cout << newCol[i] <<" ";
	// std::cout << "\n";

	return newCol;
}

void Matrix::replace_row(Matrix newRow, int nr)
{
	if(newRow.get_ncol()*newRow.get_nrow() == ncol)
		mat[slice(nr * ncol, ncol,1)] = newRow.get_valarray(); //[slice(0,ncol,1)];
	else
	{
		std::cout << "In Matrix::replace_row(). The size of new row does not match.";
		// break;
	}
}

void Matrix::replace_row(valarray<double> arr, int nr)
{
	mat[slice(nr * ncol, ncol,1)] = arr;
}

void Matrix::replace(valarray<double> newMat, int nr, int nc)
{
	mat.resize(nr*nc);
	mat = newMat;
	nrow =nr;
	ncol = nc;
}

void Matrix::print()
{
	for(int i=0; i<nrow; i++)
	{
		for(int j=0; j<ncol; j++)
			cout << Matrix::val(i,j) <<" ";
		cout <<"\n";
	}
	cout <<"\n";
}
int Matrix::get_nrow()
{
	return nrow;
}

int Matrix::get_ncol()
{
	return ncol;
}

valarray<double> Matrix::get_valarray()
{
	return mat;
}

/*
 * Computing mat1 * mat2.
 */
Matrix matrixMultiplication(Matrix mat1, Matrix mat2)
{
	Matrix new_mat;
	new_mat.initialize(mat1.get_nrow(), mat2.get_ncol());
	double new_val =0.0;

	int i,j;
	for(i=0; i< mat1.get_nrow(); i++)
		for(j=0; j< mat2.get_ncol();j++)
		{
			valarray<double> row_tmp(0.0,mat1.get_ncol());
			valarray<double> col_tmp(0.0,mat2.get_nrow());
			row_tmp = mat1.getrow(i);
			col_tmp = mat2.getcol(j);
			new_val =(row_tmp * col_tmp).sum();
			new_mat.replace(new_val,i,j);
			// new_mat.replace((mat1.getrow(i) * mat1.getcol(j)).sum(),i,j);

			//FIXME
			// cout << "i=" << i <<", j="<<j<<", new_val=" <<new_val<<"\n";
			//  for(int a =0; a<row_tmp.size();a++)
			//	cout << row_tmp[a] <<" ";
			//cout <<"\n";
			//for(int a =0; a<col_tmp.size();a++)
			//	cout << col_tmp[a] <<" ";
			//cout <<"\n";
		}

	// REMOVE
	//cout << "in matrixMultiplication():\n";
	//cout << "the multiplication of the followings:\n";
	//mat1.print();
	//mat2.print();
	//cout << "is \n";
	//new_mat.print();
	//cout << "\n";

	return new_mat;
}



/*
 * Computing pointwise matrix multiplication mat1 * mat2.
 */
Matrix matrix_pwProduct(Matrix mat1, Matrix mat2)
{
	valarray<double> val = mat1.get_valarray() * mat2.get_valarray();
	Matrix newMat;
	newMat.replace(val, mat1.get_nrow(),mat1.get_ncol());
	return newMat;
}

/**
 * Generate random number following Uniform(0,1)
 */
double runiform ()
{

  	return mt.genrand_real3();
  //   return ((double) rand() / (RAND_MAX));
}

/**
 * Generate random number following Uniform {0,1,...,upperlimit-1}
 */
int runiform_discrete(int upperlimit)
{
	return static_cast<int>(runiform() * upperlimit);
}

/**
 * Gaussian random number generator.
 * This is a copy of function 'normdev' from IMa2 code.
 * It uses the polar form of the Box-Muller transformation.
 */
// FIXME Just note. the Ziggurat algorithm is even more efficient
// than the Box–Muller transform. So, we may change this code later.
double rNormalDistribution (double mean, double stdev)
{
  static int iset = 0;
  static double gset;
  double fac, rsq, v1, v2;
  double rescale;
  if (iset == 0)
  {
    do
    {
      v1 = 2.0 * runiform() - 1.0;
      v2 = 2.0 * runiform() - 1.0;
      rsq = v1 * v1 + v2 * v2;
    }
    while (rsq >= 1.0 || rsq == 0.0);
    fac = sqrt (-2.0 * log (rsq) / rsq);
    gset = v1 * fac;
    iset = 1;
    rescale = ((v2 * fac * stdev) + mean);
    return rescale;
  }
  else
  {
    iset = 0;
    rescale = ((gset * stdev) + mean);
    return rescale;
  }
}     /// End of function 'rNormalDistribution'

/**
 * Generate random number following
 * Exponential distribution with mean 1/rate
 */
double rexpdist (double rate)
{
  return -log(runiform()) / rate;
}


unsigned int factorial(unsigned int x)
{
	if(x <= 2)
		return x;
	else
		return x*factorial(x-1);
}

double factorial(double x)
{
  if(x == 0.0)
    return 1;
  else if(x <= 2)
    return x;
  else
    return x*factorial(x-1);
}

// log(factorial(x)) = sum(log(1)+..log(x))
double logfactorial(double x)
{
  if(x == 0.0)
    return 0;
  else if(x <= 2)
    return log(x);
  else
    return log(x)+logfactorial(x-1);
}
