/* MIST copyright 2016 by Yujin Chung and Jody Hey */


/*
 * histograms.hpp
 *
 *  Created on: May 14, 2014
 *      Author: chung
 */

#ifndef HISTOGRAMS_HPP_
#define HISTOGRAMS_HPP_

#include "IM.hpp"
#include "popTree.hpp"
#include "Chain.hpp"

class plotpoint                /* (x axis are the possible values,  y axis is the counts) */
{
  double x;
  double y;
};


class hist
{	  
  char str[12];
  Eigen::MatrixXd x;
  Eigen::MatrixXd y;
  Eigen::MatrixXd freq; // used in DE to count the number of parameter values falling in the corresponding bin.
  // Eigen::MatrixXd p;
	  double yscaleadjust;
	  double xscaleadjust;
	  int before;
	  int after;
	  double mean;
	  double lower95, upper95;

public:
	  void initiateX(unsigned int gridsize, double priorMax, double priorMin);
	  void initiateY(unsigned int size){y.resize(1,size); y.setZero();}
	  void initiate_freq(unsigned int size){freq.resize(1,size); freq.setZero();}

  double getX(unsigned int i){return x(0,i);}
  double getY(unsigned int i){return y(0,i);}
  double get_freq(unsigned int i){return freq(0,i);}
  Eigen::MatrixXd getX(){return x;}
  Eigen::MatrixXd getY(){return y;}
	  double getMean(){return mean;}
	  double getLower95(){return lower95;}
	  double getUpper95(){return upper95;}

  void printHistorgram(){std::cout << "x = "<< x <<"\n" <<"y=" << y <<"\n";return;}
  void printHistorgram_prob(){std::cout << "x = "<< x <<"\n" <<"p=" << y/y.sum() <<"\n";return;}

  void addY(unsigned int i, double val){y(0,i) += val; return;}
  void replaceY(unsigned int i, double val){y(0,i) = val; return;}
  void add_freq(unsigned int i){freq(0,i) += 1; return;}
  unsigned int find_binID(double val, double numBins, double priorMax, double priorMin);
  void compute_bin_MeanPosterior_forDE();
  void computeMeanIntervals();

};

class marginal
{
  unsigned int gridsize;
  unsigned int nJointPoints;
  unsigned int nPara;
  unsigned int nZeroPara; /// the number of parameters whose upperbounds are zero.
  std::vector<double> priorsMax;
  std::vector<hist> hp;
  Eigen::MatrixXd jointDensity;
  std::vector<std::vector<double> > points2Integrate; // for each parameter, a vector of values to integrate
  std::vector<double> binWidth_small;
  std::vector<double> binWidth_large; // used to compute the marginal distribution
  
public:
  void initiate(IM im);
  void initiate_forDE(IM im, unsigned int numPara, std::vector<double> paraMax);
  void initiate_hist_forDE();
  void initiate_histXvalues_binWidth_small();
  void initiate_Points2Integrate_largeBinWidth();
  
  void set_gridsize_nJointPoints(unsigned int gs, unsigned int njp){gridsize = gs; nJointPoints = njp; return;}

  unsigned int get_nJointPoints(){return nJointPoints;}
  unsigned int get_nPara(){return nPara;}
  
  // Eigen::MatrixXd computeJointPoints_para(std::vector<unsigned int> ids_para);
  Eigen::MatrixXd compute_listJointPoints(std::vector<unsigned int> ids_para2integrate);
  void computeJointPoints();
  void computeJointDensity(IM im, popTree* poptree, Chain coldCh);
  void computeJointDensity_MPI(IM im, popTree* poptree, Chain coldCh, unsigned int nProc, unsigned int crr_procID); // parallelized over parameter values
  void computeJointDensity_MPI_overSample(IM im, popTree* poptree, Chain coldCh, unsigned int nProc, unsigned int crr_procID); // parallelized over loci
  void computeMarginalDensity(IM im, popTree* poptree, Chain coldCh);
  void computeMarginalDensity_MPI(IM im, popTree* poptree, Chain coldCh,unsigned int nProcs,unsigned int crrProcID);
  void computeMarginalDensity_MPI_overSample(IM im, popTree* poptree, Chain coldCh,unsigned int nProcs,unsigned int crrProcID);

  void add(Eigen::MatrixXd para, double posterior);

  void prepareComputingHistrograms_allProc(IM im);
  void prepareComputingHistrograms(IM im);
  void computeHistograms(IM im, popTree* poptree, Chain coldCh);
  void computeMeanIntervals();
  void computeMeanIntervals_forDE();
  
  void print_jointDensity(){ std::cout <<jointDensity <<"\n";}
  void printMeanIntervals();
  void printHistogram(unsigned int para){ hp.at(para).printHistorgram(); return;}
  void printHistogram_prob(unsigned int para){ hp.at(para).printHistorgram_prob(); return;}
  void printHistogram();
  
  void saveHistogram();
  void saveHistogram_forDE();
  void saveMeanIntervals();

  void shareJointPoints();
};


struct Param
{
  double priorMax;
  double priorMin;
  double priorMean;

  std::vector<plotpoint> xy;         /* if used points to an array of gridize
                                 * elements, each a plotpoint */
  /*
  strn str;
  int b;                        // period when this parameter first appears
  int e;                        // period when this ends
  struct weightposition wp;     // info on where the weights are for
    */                            // calculating the integration

public:
  void setup_param(IM im);
};










#endif /* HISTOGRAMS_HPP_ */
