/*MIST copyright 2016 Yujin Chung and Jody Hey */

#ifndef OPTIMIZATION_HPP_
#define OPTIMIZATION_HPP_

#include <chrono>
#include "IM.hpp"
#include "popTree.hpp"
#include "Chain.hpp"
#include "histograms.hpp"
#include "cppoptlib/solver/bfgssolver.h"
#include "cppoptlib/solver/lbfgsbsolver.h"
#include <cmath>

// Class for a maximum a posterior probability (MAP) estimate
class MaxPosterior
{
  unsigned int nPara;
  Eigen::MatrixXd para_atCrr;
  Eigen::MatrixXd para_atPrev;
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

  unsigned int optimizer;
  // 0 if differential evolution
  // 1 if L-BFGS-B solver

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
  unsigned int get_nPara(){return nPara;}
  std::vector<double> get_priorsMax(){return priorsMax;}
  Eigen::MatrixXd get_para_atCrr(){return para_atCrr;}
  Eigen::MatrixXd get_para_atPrev(){return para_atPrev;}

  double computeJointDensity_MPI_overSubSample(Eigen::MatrixXd demographicPara, IM im, popTree* poptree, Chain coldCh, unsigned int nProc, unsigned int crr_procID);
  long double computeLogJointDensity_MPI_overSubLoci(Eigen::MatrixXd demographicPara, IM im, popTree* poptree, Chain coldCh, unsigned int nProcs, unsigned int crr_procID);
  long double computeLogJointDensity_MPI_overSubLoci_ESS(Eigen::MatrixXd demographicPara, IM im, popTree* poptree, Chain coldCh, unsigned int nProcs, unsigned int crr_procID);
  long double computeLogJointDensity_mutationScalars_MPI_overSubLoci(Eigen::MatrixXd demographicPara, IM im, popTree* poptree, Chain coldCh, unsigned int nProcs, unsigned int crr_procID);
  double computeJointDensity_MPI_overSubSample_selfNormalized(Eigen::MatrixXd demographicPara, IM im, popTree* poptree, Chain coldCh, unsigned int nProcs, unsigned int crr_procID);


  void DE(IM im, popTree* poptree, Chain coldCh, unsigned int nProcs, unsigned int crr_procID);
  void DE_eachIter(IM im, popTree* poptree, Chain coldCh, unsigned int nProcs, unsigned int crr_procID);

  // void LBFGSB(IM im, popTree* poptree, Chain coldCh, unsigned int nProcs, unsigned int crr_procID);

};




template<typename T, int D>
class computeLogJointDensity : public cppoptlib::BoundedProblem<T, D>
{
public:
  using typename cppoptlib::Problem<T, D>::TVector;
  // using MatrixType = Eigen::MatrixXd<T, Eigen::Dynamic, Eigen::Dynamic>;

protected:
  const unsigned int currentProcessID;
  Chain ch;
  const unsigned int numProcesses;
  popTree* poptree;
  const IM im;
  MaxPosterior MAPestimate;
  
public:
    computeLogJointDensity(MaxPosterior &MAPestimate_, const IM &im_, popTree* &poptree_, Chain &ch_, const unsigned int &numProcesses_, const unsigned int &currentProcessID_): MAPestimate(MAPestimate_), im(im_), poptree(poptree_), ch(ch_), numProcesses(numProcesses_), currentProcessID(currentProcessID_) {}
  
  // objective
  T value(const TVector &paraVec) 
  {
    double logJointDensity = 0;

    Eigen::MatrixXd para;
    para.resize(1, MAPestimate.get_nPara());
    for(unsigned int i=0; i< paraVec.size(); i++)
      para(0,i) = paraVec(i);

    if(ch.get_multiLocusSpecific_mutationRate() ==1) // variable mutation rate scalars
      {
	    logJointDensity =  -1 * (double) MAPestimate.computeLogJointDensity_mutationScalars_MPI_overSubLoci(para, im, poptree, ch, numProcesses, currentProcessID); 
      }
    else // constant mutation rate scalars
      {
	    logJointDensity =  -1 * (double) MAPestimate.computeLogJointDensity_MPI_overSubLoci(para, im, poptree, ch, numProcesses, currentProcessID); 
	    }	

//////////////////// xing, 11/29/2017, comment these cout
////////////////// The following cout would be called whenever the computeLogJointDensity object is instantiated, and it would be created many many times (whenever the value needs to be computed). 
////////////////// Do not need to print the value here, for comparison purpose, it can be printed in lbfgsbsolver
    // std::cout <<" paraVect = " << paraVec <<"\n";
//    std::cout <<"para = "<< para <<"\n";
//    std::cout <<" logJointDensity = "<< logJointDensity <<"\n\n";
    
    if(std::isnan(logJointDensity)==1)
      logJointDensity = numeric_limits<double>::infinity();

    return logJointDensity;
  }
  /*
  T value(const TVector &paraVec) 
  {
    double logJointDensity = computeLogJointDensity_mutationScalars_MPI_overSubLoci(paraVec, MAPestimate, im, poptree, ch, numProcesses, currentProcessID);
    return logJointDensity;
  }
  */

  // double computeLogJointDensity_mutationScalars_MPI_overSubLoci(TVector paraVec, MaxPosterior MAPestimate, IM im, popTree* poptree, Chain coldCh, unsigned int nProcs, unsigned int crr_procID);
};

// typedef computeLogJointDensity<double, > TcomputeLogJointDensity;
// typedef typename computeLogJointDensity<double, int>::TVector TVector;


#endif /* OPTIMIZATION_HPP_ */
