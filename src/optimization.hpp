/*MIST copyright 2016 Yujin Chung and Jody Hey */

#ifndef OPTIMIZATION_HPP_
#define OPTIMIZATION_HPP_

#include <chrono>
#include "IM.hpp"
#include "popTree.hpp"
#include "Chain.hpp"
#include "histograms.hpp"

// Class for a maximum a posterior probability (MAP) estimate
class MaxPosterior
{
  unsigned int nPara;
  Eigen::MatrixXd para_atCrr;
  Eigen::MatrixXd para_atPrev;
  
  double durationOfSplitting_atPrev;
  double durationOfSplitting_atCrr;
  
  std::vector<double> priorsMax;
  std::vector<long double> posterior_atCrr;
  std::vector<long double> posterior_atPrev;
  double diff_posteriors;
  long double logPosteriorMax;
  long double logPosteriorMin;
  double maxDist; // the maximum distance between the current optimum and others
  unsigned int ID_max;

  //-- for differential evolution --//
  unsigned int nParaVectors;
  unsigned int nIters;
  // YC 1/12/2015
  // The parameters, F and CR, for differential evolution
  // are the same as those of IMa2
  double F; // the differential weigth [0,2]
  double CR; // the cross over probability [0,1]
  double a,b,c; // three agents

  unsigned int howOften_checkpoint;
  // how often to write a checkpoint.
  unsigned int checkpoint;
  // 0 for no checkpoint (default)
  // 1 for writing checkpoint
  // 2 for reading the checkpoint and starting from it but not writing checkpoint anymore
  // 3 for reading the checkpoint and writing checkpoint.


  marginal marginals;


  
  std::chrono::microseconds totalComputingTime_eigen;
  std::chrono::microseconds totalComputingTime_eigen_subMatOfEigenVectors;
  std::chrono::microseconds totalComputingTime_condiProb;
  unsigned int totalNum_eigenFunctionCalls;
  unsigned int totalNum_condiProbFunctionCalls;

 
public:
  void initiate(IM im, popTree* poptree, Chain coldCh, unsigned int nProcs, unsigned int crr_procID);
  void read_checkpoint();
  void write_checkpoint();

  std::chrono::microseconds get_totalComputingTime_eigen(){return totalComputingTime_eigen;} 
  std::chrono::microseconds get_totalComputingTime_condiProb(){return totalComputingTime_condiProb;}
  unsigned int get_totalNum_eigenFunctionCalls(){return totalNum_eigenFunctionCalls;}
  unsigned int get_totalNum_condiProbFunctionCalls(){return totalNum_condiProbFunctionCalls;}

  double computeJointDensity_MPI_overSubSample(Eigen::MatrixXd demographicPara, IM im, popTree* poptree, Chain coldCh, unsigned int nProc, unsigned int crr_procID);
  long double computeLogJointDensity_MPI_overSubLoci(Eigen::MatrixXd demographicPara, IM im, popTree* poptree, Chain coldCh, unsigned int nProcs, unsigned int crr_procID);
  long double computeLogJointDensity_MPI_overSubLoci_ESS(Eigen::MatrixXd demographicPara, IM im, popTree* poptree, Chain coldCh, unsigned int nProcs, unsigned int crr_procID);
  long double computeLogJointDensity_mutationScalars_MPI_overSubLoci(Eigen::MatrixXd demographicPara, IM im, popTree* poptree, Chain coldCh, unsigned int nProcs, unsigned int crr_procID);
  double computeJointDensity_MPI_overSubSample_selfNormalized(Eigen::MatrixXd demographicPara, IM im, popTree* poptree, Chain coldCh, unsigned int nProcs, unsigned int crr_procID);


  void DE(IM im, popTree* poptree, Chain coldCh, unsigned int nProcs, unsigned int crr_procID);
  void DE_eachIter(IM im, popTree* poptree, Chain coldCh, unsigned int nProcs, unsigned int crr_procID);


};


#endif /* OPTIMIZATION_HPP_ */
