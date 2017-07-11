/* MIST copyright 2016 by Yujin Chung and Jody Hey */

#ifndef _Chain_
#define _Chain_

#include <vector>
#include <string>
//FIXME
	// For MPI, we need to uncommand out
// #include <gsl/gsl_rng.h>
// #include <gsl/gsl_randist.h>
#include <math.h>
#include <fstream>
#include <time.h>
#include <chrono>
#include "IM.hpp"
#include "coaltree.hpp"

class Chain
{
private :
  //--- For M and L mode ---//
  unsigned int n_loci; // Total number of loci
  unsigned int n_MCMCgen; 
  //  number of sampled (saved or read) iterations/generations to run the chain
  unsigned int lociInParallel; 
  // 1 if loci are in parallel
  // 0 if chains are in parallel in M mode; 0 if samples are in parallel in L mode
  unsigned int locusID_start; 
  // If 'lociInParallel' =0, then 'locusID_start' = 0.
  unsigned int locusID_end; 
  // If 'lociInParallel' =0, then locusID_end = n_loci-1.
  unsigned int numSubLoci; 
  // If 'lociInParallel' =0, then numSubLoci = n_loci
  unsigned int nPairsMut;
  unsigned int nProcesses; // The total number of processes.
  unsigned int cpuID; 

  unsigned int nChains; // The number of chains per process. M mode only.


  
  unsigned int priorType; // The kind of importance sampling distribution
  // 1 if one peak
  // 2 if two peaks
  // 3 if 3 peaks
  // 4 if 4 peaks
  double changeOfRatePoint;

  std::vector<double> rNormalDist; // for debugging only

 
  unsigned int totalNum_coalEvents; // Total number of coalescent events from all the loci.
  std::vector<nodeSimple*> list_trees;
  std::vector<std::vector<unsigned int> > treeIDs; // shared across processes
  std::vector<node*> trees_atPrev;
  std::vector<std::vector<std::list<double> > > coalTimes;  // shared across processes
  // coalescent times (the age of node) from top to bottom
  std::vector<std::vector<std::vector<unsigned int> > > tipIDs;
  std::vector<std::vector<double> > kappa; /// Rate of transitions and transversions
  std::vector<double> kappa_atPrev;
  std::vector<std::vector<double> > mutationScaler; /// Mutation rate scalers
  std::vector<double> mutationScaler_atPrev;
  std::vector<std::vector<double> > logLikelihood;
  std::vector<std::vector<double> > logPrior_each;
  std::vector<double> logLikelihood_atPrev;
  std::vector<double> logPrior_trees; // shared across processes in L mode  
  std::vector<double> logPriorTrees_atPrev_lociParallel;
  std::vector<double> logJointLikelihood;
  double logJointLikelihood_atPrev_lociParallel;
  double logPriorTrees_atPrev;
  double popSizeMax;
  // FIXME remove mutrate
  double mutrate;
  double temperature; // The heating parameter in Metropolis-coupled MCMC. (0 < temperature < 1)
  /// Note that the cold chain (chainID = 0) has heat 1.
  /// For heated chain i,
  ///    "temperature" = 1/(1 + deltaT * (i-1)),
  /// where deltaT (called "temperature" in some literatures) controls
  /// the amount of the smoothness of the posterior distribution.
  /// That is, larger delta T, smaller temperature, larger posterior probability
  /// - YC 3/5/2014

  // std::vector< std::vector< unsigned int>> pairIDsMut;
  // the IDs of pairs of mutation scalars to update.
  // For an odd j, pairIDsMut.at(0).at(j) and pairsIDsMut.at(0).at(j+1) try to update at the current iteration.
  std::vector< unsigned int> pairIDsMut;

  double splittingTimeMax; /// the upper bound of the splitting times

  unsigned int rankTemp; // The rank of temperature. Cold chain has rank 0. Added by YC 9/12/2014

  
  double partialPosterior_log_numerator;
  double partialPosterior_log_denominator;
  
  double prior;
  double posterior;
  unsigned int chainID;
  unsigned int current_iteration;

  // Added by YC 8/7/2014
  Eigen::MatrixXi numTries_trees;
  Eigen::MatrixXd numAccpt_trees;

  //clock_t clock_start, clock_end;
  //double totalComputingTime_eigen;
  //unsigned long long totalComputingTime_condiProb;
  std::chrono::microseconds eachComputingTime_eigen;
  std::chrono::microseconds eachComputingTime_eigen_subMatOfEigenVectors;
  std::chrono::microseconds eachComputingTime_condiProb;
  std::chrono::microseconds eachComputingTime_condiProb_case1;
  std::chrono::microseconds eachComputingTime_condiProb_case2_1;
  std::chrono::microseconds eachComputingTime_condiProb_case2_2;
  std::chrono::microseconds eachComputingTime_condiProb_case3_1;
  std::chrono::microseconds eachComputingTime_condiProb_case3_2;
  unsigned int count_condiProbFunctionCalls;
  unsigned int count_case1;
  unsigned int count_case2_1;
  unsigned int count_case2_2;
  unsigned int count_case3_1;
  unsigned int count_case3_2;




  // ---- members used in the L mode only ----- //
  unsigned int numTotalSeq; // the number of total seuqneces 
  std::vector<unsigned int> SeqPop;
  unsigned int nTrees; // the number of trees per locus to use in L mode
  // this number should be smaller or equal to n_MCMCgen
  unsigned int nSubSample; // MCMC sample
  unsigned int sampleID_start;
  unsigned int sampleID_end;  
  unsigned int TreesWithMaxP;
  // 1 if trees with highest probabilities are selected 
  // 0 if some first trees are selected.

  int topoID_start;
  int topoID_end;

  unsigned int multiLocusSpecific_mutationRate; // 0 if constant mutation rate; 1 if mutation rate variation over loci.

  // The following objects are used when "TreesWithMaxP" =1.
  std::vector<double> MaxLogLik; 
  std::vector<unsigned int> MaxSampleID; // contains which samples are selected for L mode 

  // the followings are computed only once
  std::vector<std::vector<std::vector<double> > > summaryOfCoalTree;
  //-- Definition of "summaryOfCoalTree" --//
  // "summaryOfCoalTree" is computed only when the upper bound of migration rates is 0
  // For each unique tree topology, for each coalescent event of the corresponding tree,
  //    
  std::vector<std::vector<std::vector<unsigned int> > > states_observed; 
  std::vector<std::vector<std::vector<unsigned int> > > states_observed_freq;
    // For each given topology, for each time period (from tips to root) determined by coalescent events, 
    // 'states_observed' and 'states_observed_freq' contains the kinds of ids of lineages (rank or population id) and their frequency.
    // For example, if the first topology in 'list_trees' is ((1,1):5,2):6, 
    // then states_observed.at(0).at(0) = (1,2);
    //      states_observed.at(0).at(1) = (5,2);
    //      states_observed.at(0).at(2) = (6);
    // and states_observed_freq.at(0).at(0) = (2,2);
    //     states_observed_freq.at(0).at(1) = (1,1);
    //     states_observed_freq.at(0).at(2) = (1).
  std::vector<std::vector<unsigned int> > nKinds_lineages; // shared across processes
    // For each given topology, for each time period (from tips to root) determined by coalescent events, 
    // 'nKinds_lineages' contains the number of kinds of lineages. 
    // For example, if the first topology in 'list_trees' is ((1,1):5,2):6 (same as the above example),
    // then nKinds_lineages.at(0) = (2,2,1).
  std::vector<std::vector<Eigen::MatrixXd> > stateSpaces; // shared across processes
  std::vector<std::vector<std::vector<std::vector<unsigned int> > > > possiblePaths; // shared across processes
  std::vector<std::vector<std::vector<std::vector<unsigned int> > > >nPossibleCoalEvents; // for each unique tree topology, for each coalescent event, for each population where coalescent events happen, 
  std::vector<unsigned int> numTrees; // the size of an equivalence class.
  std::vector<std::vector<unsigned int> > numSameCoalEvents_inAncPop; // shared across processes
  
  std::vector<std::vector<std::vector<Eigen::MatrixXd> > > coeff4TransitionRateMat;
  // coeff4TransitionRateMat.at(i).at(j).resize(0): coefficients for pop1 size
  // coeff4TransitionRateMat.at(i).at(j).resize(1): coefficients for pop2 size
  // coeff4TransitionRateMat.at(i).at(j).resize(2): coefficients for migration rate 1
  // coeff4TransitionRateMat.at(i).at(j).resize(3): coefficients for migration rate 2

  // the followings are computed once whenever a different population tree is considered
  std::vector<std::vector<std::vector<Eigen::MatrixXcd> > > Qeigen;
  std::vector<std::vector<std::vector<Eigen::MatrixXcd> > > subMatV;
  // for ith topology, for jth coalescent event (backward in time)
  // subMatV.at(i).at(j).at(k):
  //  if k=0, submatrix of eigen vector matrix for possible paths from the previous event
  //  if k=1, submatrix of inverse eigenvector matrix for a coalescent event (possible paths of a coalescent event)
  //  if k=2, submatrix of inverse eigenvector matrix for population merging (all transient paths)
  std::vector<std::vector<Eigen::ArrayXcd> > eigenvalues;
  // for ith topology, for jth coalescent event (backward in time)
  // eigenvalues.at(i).at(j): eigenvalues

  std::vector<long double> expectationOfCoalProb; // L mode: parallel over sample 
  std::vector<long double> logExpectationOfCoalProb; // L mode: parallel over over loci only if same mutation rates
  std::vector<long double> logExpectationOfCoalProbSquared; // L mode: for ESS of importance sampling, parallel over over loci only if same mutation rates;
  std::vector<long double> partialLogJointCoalProb; // L mode: parallel over loci and variable mutation scalars

  double percentile_4upperBoundOfTMRCA;
  double upperBoundOfTMRCA;
  unsigned int nSamples_belowTheUpperBoundOfTMRCA;
  unsigned int nSamples_twoLargeImporatanceRatio;

  // YC 5/4/2017
  unsigned int Forest;
  // 0 if a full coalescent tree and population tree are considered (default)
  // 1 if the forest of coalescent trees are condisered
  // 2 if the forest of population trees are considered
  // 3 if the forests of coalescent trees and population trees are considered.
  unsigned int sizeCoalsubtree; // the size of subtrees of coalescent trees (default=0)
  unsigned int sizePopsubtree;  // the size of subtrees of population trees (default=0)
  unsigned int sizeForest; // the number of subtrees per locus per sample.
  
  // subtrees in Step 2 - YC 3/8/2017
  std::vector<std::vector<std::vector< node*> > > Coalsubtrees;
  // per-sample, per-locus subtrees
  // for example, if the MCMC sample size is 1, then subTrees.at(0) contains all subtrees of loci.
  // std::vector<nodeSimple*> listSubtrees; // the list of topo for subtree-version is the same as list_trees
  std::vector<std::vector<std::vector<unsigned int> > > subtreeIDs; // shared across processes
  // per-sample, per-locus
  std::vector<std::vector<std::vector<std::list<double> > > > subcoalTimes;  // shared across processes
  std::vector<std::vector<std::vector<std::vector<unsigned int> > > > subtipIDs;
  
  
public:
  ~Chain() {};
  void InitializeChain(IM im, int nMCMCgen, int c_id, unsigned int processID, unsigned int n_chains, unsigned int numprocesses);
  void InitializeChain_usePriorOnly(IM im, int nMCMCgen, int c_id, unsigned int processID, unsigned int n_chains, unsigned int numprocesses);   
  void deleteAll();
  void deleteAll_Lmode();
  void deleteTrees();
  void deleteListTrees();
  void delete_temporaryPara();
  
  double GetTemperature() { return temperature;}
  std::vector<double> getTreePriors() {return logPrior_trees;}
  std::vector<std::vector<double> > getKappa(){return kappa;}
  std::vector<std::vector<double> > getMutationScaler(){return mutationScaler;}
  double GetPosterior() { return posterior;}
  int GetChainID() { return chainID;}
  int GetCurrentIteration() { return current_iteration;}
  int GetNumIterations() { return n_MCMCgen; }
  unsigned int GetNumLoci() { return n_loci; }
  unsigned int getNumSubLoci(){return numSubLoci;}
  unsigned int get_nSubSample(){return nSubSample;}
  std::vector<nodeSimple*> GetListTrees() {return list_trees;}
  std::list<double> GetCoalTimes(unsigned int id_crrIter, unsigned int i) { return coalTimes.at(id_crrIter).at(i);}
  node* GetTreesAtPrev(unsigned int loc) { return trees_atPrev.at(loc);}
  unsigned int GetTreeIDs(unsigned int id_crrIter, unsigned int i) { return treeIDs.at(id_crrIter).at(i);}
  double GetLogLikelihoodAtPrev(unsigned int loc) { return logLikelihood_atPrev.at(loc);}
  std::vector<double> getLogLikelihoodAtPrev(){return logLikelihood_atPrev;}
  double GetLogPrior(int locus_id) { return logPrior_trees.at(locus_id);} 
  double GetLogPriorAtPrev() { return logPriorTrees_atPrev;} 
  double getLogPriorTrees_atPrev_lociParallel(unsigned int locID){return logPriorTrees_atPrev_lociParallel.at(locID);}
  double getLogJointLikelihood_atPrev_lociParallel(){return logJointLikelihood_atPrev_lociParallel;}
  std::vector<double> GetLogPrior() { return logPrior_trees;}
  std::vector<double> GetKappaAtPrev() { return kappa_atPrev;}
  double getKappaAtPrev(unsigned int loc){return kappa_atPrev.at(loc);}
  std::vector<double> GetMutationScalerAtPrev() { return mutationScaler_atPrev;}
  double getMutationScalerAtPrev(unsigned int loc) { return mutationScaler_atPrev.at(loc);}
  unsigned int getRankTemp(){return rankTemp;}	
  std::chrono::microseconds get_eachComputingTime_eigen(){return eachComputingTime_eigen;}
  std::chrono::microseconds get_eachComputingTime_condiProb(){return eachComputingTime_condiProb;}
  unsigned int get_countCondiProbFunctionCalls(){return count_condiProbFunctionCalls;}    
  double get_partialPosterior_log_numerator(){return partialPosterior_log_numerator;}
  double get_partialPosterior_log_denominator(){return partialPosterior_log_denominator;}
  std::vector<long double> get_expectationOfCoalProb(){return expectationOfCoalProb;}
  std::vector<long double> get_logExpectationOfCoalProb(){return logExpectationOfCoalProb;}
  std::vector<long double> get_logExpectationOfCoalProbSquared(){return logExpectationOfCoalProbSquared;}
  unsigned int get_multiLocusSpecific_mutationRate(){return multiLocusSpecific_mutationRate;}
  std::vector<long double> get_partialLogJointCoalProb(){return partialLogJointCoalProb;}

  

  void resizeLogLikelihood(unsigned int iterID){logLikelihood.at(iterID).resize(n_loci); return;}
  void resizeTreeIDs(unsigned int iterID){treeIDs.at(iterID).resize(n_loci); return;}
  void resizeCoalTimes(unsigned int iterID){coalTimes.at(iterID).resize(n_loci); return;}
  void resizeTipIDs(unsigned int iterID){tipIDs.at(iterID).resize(n_loci); return;}
  void resizeKappa(unsigned int iterID){kappa.at(iterID).resize(n_loci); return;}
  void resizeMutationScaler(unsigned int iterID){mutationScaler.at(iterID).resize(n_loci); return;}

  void SetChainID(int c_id);
  void SetListTrees(nodeSimple* topo) { list_trees.push_back(topo); return;}
  void SetTreeIDs(unsigned int id_crrIter, unsigned int i, unsigned int value) { treeIDs.at(id_crrIter).at(i) = value; return;}
  void SetLogLikelihood(unsigned int id_crrIter, unsigned int i, int value) {logLikelihood.at(id_crrIter).at(i) = logLikelihood_atPrev.at(i); return;}
  void SetLogLikelihood(unsigned int savingID, std::vector<double> logLikeValues){logLikelihood.at(savingID) = logLikeValues; return;}
  void SetLogLikelihood_lociMPI(unsigned int savingID, unsigned int lcID, double logLikeValues){logLikelihood.at(savingID).at(lcID) = logLikeValues; return;}
  void setEachLogPrior_lociMPI(unsigned int savingID, unsigned int lcID, double logprior){logPrior_each.at(savingID).at(lcID) = logprior; return;}
  void SetLogPriors(unsigned int id_crrIter, double logpriorprev1) {logPrior_trees.at(id_crrIter) = logpriorprev1; return;}
  void SetLogJointLike(unsigned int id_crrIter, double loglik) {logJointLikelihood.at(id_crrIter) = loglik; return;}
  void SetKappa(unsigned int id_crrIter, std::vector<double> kappa1 ) {kappa.at(id_crrIter) = kappa1; return;}
  void setKappa(unsigned int id_crrIter, unsigned int loc, double newKappa ) {kappa.at(id_crrIter).at(loc) = newKappa; return;}
  void SetMutationScaler(unsigned int id_crrIter, std::vector<double> mutationscaler1) {mutationScaler.at(id_crrIter) = mutationscaler1; return;} 
  void setMutationScaler(unsigned int id_crrIter, unsigned int loc, double newScaler) {mutationScaler.at(id_crrIter).at(loc) = newScaler; return;}
  void setZero_numTries_trees(){numTries_trees.setZero(); return;}
  void setZero_numAccpt_trees(){numAccpt_trees.setZero(); return;}
  void SetNumLoci(int n_loci);
  void SetTemperature(double n_temperature);
  void SetMutationRate(double mut_rate);
  void SetKappa();
  void SetNumIterations(int n_gen);
  void setCoalTimes(unsigned int id_crrIter, unsigned int loc, std::list<double> newCoalTime){coalTimes.at(id_crrIter).at(loc) = newCoalTime; return;}
  void setTipIDs(unsigned int id_crrIter, unsigned int loc, std::vector<unsigned int> newTipIDs){tipIDs.at(id_crrIter).at(loc) = newTipIDs; return;}
  void set_splittingTimeMax(double t){splittingTimeMax = t; return;}
  void set_lociInParallel(unsigned int x){lociInParallel=x;return;}
  void set_locusID_start(unsigned int x){locusID_start=x;return;}
  void set_locusID_end(unsigned int x){locusID_end=x;return;}
  void set_numSubLoci(unsigned int x){numSubLoci=x;return;}
  void set_nPairsMut(unsigned int x){nPairsMut=x;return;}


	void PrintChain();
	void print_chain_initial();
	void print_trees();
	void print_trees_lik();
	void print_trees_atIter(unsigned int iter);
	void print_states_atIter(unsigned int iter);
	void fileprint_states_atIter(unsigned int iter, std::ofstream& fp);
	void print_listTrees();
	void fileprint_listTrees(std::ofstream& fp);
	void print_coalTimes(unsigned int iter, unsigned int id_locus);
	void print_tipIDs(unsigned int iter, unsigned int id_locus);
  void print_savedStates(unsigned int iter);
  void print_treeAccptRate();
  std::string print2string_states_atCrr();

  void saveTrees2File();
  void saveOriginalTreeTopo2File();
  void saveCoalTimes2File();
  void saveKappa2File();
  void saveMutationScalers2File();
  void saveTipIDs();
  void saveTreeIDs();
  void saveListTopo2File();
  void saveLogPriorTrees();
  void saveLogJointLikelihoods();
  void saveLogEachLikelihoods();
  void saveEachLogPrior();
  void saveRNormalDist();

  void read_logJointPriorTrees();
  void read_treeIDs();
  void read_coalTimes(IM im);
  void read_kappa();
  void read_mutationScaler();
  void read_tipIDs(IM im);
  void read_listTrees();
  void read_logEachLikelihood_partial();
  void read_logJointLikelihood_partial();

  void read_LmodeInputFile(IM im);
  void read_eachLogPrior_partial();
  void read_logJointPriorTrees_partial();
  void read_treeIDs_partial();
  void read_coalTimes_partial(IM im);
  void read_tipIDs_partial(IM im);
  void read_mutationScaler_partial();
  void read_newickTrees(IM im);
  void read_newickTrees_partial(IM im);

  void initialize_trees_atPrev(vector<locus> lc);

  void FillUpTreesIterations();
  
  void UpdateChain(unsigned int id_crrIter, vector<locus> lc,unsigned int save, int currentid, int numprocesses, int chainid, unsigned int n_chains);
  //void UpdateChain_usePriorsOnly(unsigned int id_crrIter, vector<locus> lc,unsigned int save, int currentid, int coldchain, int numprocesses, int chainid);
  void UpdateChain_usePriorsOnly(unsigned int id_crrIter, vector<locus> lc,unsigned int save, int currentid, int numprocesses, int chainid, unsigned int n_chains);
  void UpdateTrees(unsigned int id_crrIter, vector<locus> lc,unsigned int save);
  
  void UpdateTrees_newProposal(unsigned int id_crrIter, vector<locus> lc, unsigned int save);
  void UpdateTrees_usePriorsOnly(unsigned int id_crrIter, vector<locus> lc,unsigned int save);
  void UpdateTrees_usePriorsOnly_newProposal(unsigned int id_crrIter, vector<locus> lc, unsigned int save);
#ifdef MPI_ENABLED
  // Replaced "UpdateTrees_MPI" updated by AS.
  void UpdateTrees_MPI_newProposal(unsigned int id_crrIter, vector<locus> lc, unsigned int save, int currentid, int chainid);
  void UpdateTrees_LociInMPI_newProposal(unsigned int id_crrIter, vector<locus> lc, unsigned int save, int currentid, int chainid);
  int UpdateTrees_MPI_newProposal_old(unsigned int id_crrIter, vector<locus> lc, unsigned int save, int currentid, int chainid, int numprocesses, int coldchain);
  // int UpdateTrees_MPI_usePriorsOnly_newProposal(unsigned int id_crrIter, vector<locus> lc, unsigned int save, int currentid, int chainid);
  void UpdateTrees_MPI_usePriorsOnly_newProposal(unsigned int id_crrIter, vector<locus> lc, unsigned int save, int currentid, int chainid);
  int UpdateTrees_MPI(unsigned int id_crrIter, vector<locus> lc, unsigned int save, int currentid, int chainid, int numprocesses, int coldchain);
  int UpdateTrees_MPI_usePriorsOnly(unsigned int id_crrIter, vector<locus> lc, unsigned int save, int currentid, int chainid);

  // void samplingUpdatedTrees_MPI(unsigned int id_crrIter, vector<locus> lc, unsigned int save, int currentid, int chainid, int coldchain)

	#endif

  void updateKappa(unsigned int id_crrIter, locus lc, unsigned int save, int currentid); //, int doisendinfo, int chainid, int coldchain);
  void updateKappa_usePriorsOnly(unsigned int id_crrIter, locus lc, unsigned int save, int currentid); // (unsigned int id_crrIter, locus lc, unsigned int save, int currentid, int doisendinfo, int chainid);
  void update_mutationScaler_Kappa(unsigned int id_crrIter, std::vector<locus> lc, unsigned int crrProcID, unsigned int nProc);
  void update_mutationScalar_lociInParallel(unsigned int id_crrIter, std::vector<locus> lc, unsigned int crrProcID, unsigned int nProc);
  void update_randomPairsMutationScalars_lociInParallel(unsigned int id_crrIter, std::vector<locus> lc, unsigned int crrProcID, unsigned int nProc);
void update_mutationScaler_Kappa_usePriorsOnly(unsigned int id_crrIter, std::vector<locus> lc);//(unsigned int id_crrIter, std::vector<locus> lc,unsigned int save, int currentid, int doisendinfo, int chainid);

  void collectAllUpdates(unsigned int savingID, unsigned int Lmode); // 1 process and 1 chain
  void collectAllUpdates_Lmode(unsigned int savingID);
  void collectAllUpdates_Lmode_subtrees( unsigned int savingID);
  void collectAllUpdates_bwProcs_Lmode();
  void collectAllUpdates_bwProcs_Lmode_subtrees();


  void AddTrees(std::vector<node*> trees);
  double logLikelihoodHKY(locus lc, int which_iteration, int which_locus);

  void compute_jointPrior_indepLoci();
  void compute_logLikelihood_indepLoci();
  void compute_rankTemp();
  void compute_totalnum_coalEvents(std::vector<locus> loci);
  void compute_coalTimes_tipIDs(unsigned int iter, unsigned int id_locus, node* tree);
  void compute_coalTimes_tipIDs_subtrees(unsigned int iter, unsigned int id_locus, unsigned int id_subtr, node* tree);

  double compute_logPriorTrees(node* newTrees, unsigned int id_locus);
  double compute_logPrior_usingTrees_atPrev();
  double compute_totalLogLik();
  
  double MHratio_trees(unsigned int id_locus, node* newTree, double newlogLik, double logRatioPriors, double slideDist);
  double MHratio_trees_usePriorOnly(unsigned int id_locus, node* newTree, double newlogLik, double logRatioPriors, double slideDist);
  double MHratio_trees_newProposal(unsigned int id_locus, node* newTree, double newlogLik, double logRatioPriors, double slideDist);
  double MHratio_trees_usePriorOnly_newProposal(unsigned int id_locus, node* newTree, double newlogLik, double logRatioPriors, double slideDist);
  
  Eigen::MatrixXi get_numTries_trees(){return numTries_trees;}
  Eigen::MatrixXd get_numAccpt_trees(){return numAccpt_trees;}
  
  //AS: adding set functions for sharecoaltree
  

  double compute_sumOfInversePrior();

  // ----- L mode ------ //
  void initializeLmode(IM im, unsigned int crrProcID, unsigned int nProcs);

  // the followings are called once
  void compute_observedStates_fromSubtrees(unsigned int id_listTopo, unsigned int id_coalEvent, std::vector<nodeSimple*> subtrees);
  void compute_observedStates_fromTopo();
  Eigen::MatrixXd compute_stateSpaces_recursion(std::vector<unsigned int> freq,unsigned int nPops);
  void compute_stateSpaces(unsigned int nPops);
  unsigned int find_initialState(unsigned int nPops);
  Eigen::MatrixXd find_minimumFreq_forNextTimePeriod(unsigned int id_samples,unsigned int  id_period);
  void compute_possiblePaths(unsigned int nPops);
  void compute_possiblePaths_nPossibleCoalEvents(unsigned int nPops);
  void compute_numTrees(unsigned int nPops);
  void compute_numSameCoalEvents_inAncPop();
  void prepare_Lmode(popTree* poptree, IM im);
  
  
  // the followings are called whenever different population trees are considered
  void compute_eigenValuesVectors_subMatOfEigenVectors(popTree* poptree);
  void compute_eigenValuesVectors_subMatOfEigenVectors_MPI(popTree* poptree);
  void compute_summaryOfCoalTrees(popTree* poptree);
  void compute_eigenValuesVectors(popTree* poptree);
  void compute_eigenValuesVectors_rateMat(popTree* poptree); // old version
  void compute_eigenValuesVectors_rateMat_MPI(popTree* poptree, unsigned int crrProc_ID); // old version
  void compute_eigenValuesVectors_rateMat_MPI(popTree* poptree); // old version
  void compute_coefficientsOfElementsInTransitionRateMat(unsigned int nPops);
  void share_treeIDs_coalTimes_logPriorTrees(unsigned int crrProc_ID);
  void share_stateSpaces_nKindsLineages_possiblePaths_numTrees(unsigned int crrProc_ID);
  void share_stateSpaces_nKindsLineages_possiblePaths_numTrees_nPossibleCoalEvents(unsigned int crrProc_ID);
  void share_eigenValuesVectors(unsigned int crrProc_ID);
  double compute_conditionalProb(unsigned int id_sample, unsigned int id_locus, popTree* poptree); 
  double compute_conditionalProb_MPI(unsigned int id_sample, unsigned int id_locus, popTree* poptree, unsigned int crrProcID); 
  double compute_logConditionalProb(unsigned int id_sample, unsigned int id_locus, popTree* poptree, unsigned int crrProcID); // currently used
  double compute_logConditionalProb_subtrees(unsigned int id_sample, unsigned int id_locus, popTree* poptree, unsigned int crrProcID);
  long double compute_logConditionalProb_zeroMig(unsigned int id_sample, unsigned int id_locus, popTree* poptree, unsigned int crrProcID); // currently used
  long double compute_logConditionalProb_zeroMig_subtrees(unsigned int id_sample, unsigned int id_locus, popTree* poptree, unsigned int crrProcID);
  double compute_logConditionalProb_old(unsigned int id_sample, unsigned int id_locus, popTree* poptree, unsigned int crrProcID); 
  double compute_jointPosteriorDensity(popTree* poptree,IM im);
  double compute_jointPosteriorDensity_MPI(popTree* poptree, IM im, unsigned int crrProcID);
  double compute_partialJointPosteriorDensity_MPI_overSample(popTree* poptree, IM im, unsigned int crrProcID, unsigned int nProcs); // old version
  double compute_partialJointPosteriorDensity_overSubSample_old(popTree* poptree, IM im, unsigned int crrProcID, unsigned int nProcs);
  double compute_partialJointPosteriorDensity_overSubSample_expectationOfJoint(popTree* poptree, IM im, unsigned int crrProcID, unsigned int nProcs);
  void compute_partialJointPosteriorDensity_overSubSample(popTree* poptree, IM im, unsigned int crrProcID, unsigned int nProcs);
  void compute_partialJointPosteriorDensity_overSubLoci(popTree* poptree, IM im, unsigned int crrProcID, unsigned int nProcs);
  void compute_partialJointPosteriorDensity_overSubLoci_ESS(popTree* poptree, IM im, unsigned int crrProcID, unsigned int nProcs);
  void compute_partialJointPosteriorDensity_mutationScalars_overSubLoci(popTree* poptree, IM im, unsigned int crrProcID, unsigned int nProcs);
  void compute_partialJointPosteriorDensity_overSubSample_selfNormalized(popTree* poptree, IM im, unsigned int crrProcID, unsigned int nProcs);
  void find_print_theMaxHeightsTrees();

  // subtrees 
  void getSubtree(unsigned int subSize, node* tree, unsigned int sampleID, unsigned int subLocusID);
  
	// old version
	unsigned int compute_size_stateSpaces(unsigned int freq, unsigned int nPops);
	std::string find_nextCoalLineage(unsigned int id_samples,unsigned int  id_period);

	// --- End of the list of functions used in L mode ----//
};



class MCMC
{
private:
  unsigned int nProcesses;
  unsigned int n_chains; // the number of chains per process - YC 5/29/2014
  std::vector<Chain> chains;
  unsigned int process_id;
  unsigned int n_MCMCgen;
  unsigned int nBurning;
  unsigned int thinning;
  unsigned int samplingFromPriorOnly;
  double upperkappa;

  unsigned int lociInParallel; // for M mode only
  // 1 if loci are in parallel
  // 0 if chains are in parallel
  unsigned int locusID_start;
  unsigned int locusID_end;
  unsigned int numSubLoci;
  unsigned int nPairsMut;
  
  // unsigned int coldChainID;
  
  Eigen::MatrixXd AccepRate_swapChains;
  std::vector<unsigned int> order_globalChainIDs; // the order of temperatures of chains
  
  Eigen::MatrixXd AccepRate_swapChains_local;
  Eigen::MatrixXd AccepRate_swapChains_global;
  
  
  Eigen::MatrixXd numAccpt_trees_allChains_local;
  Eigen::MatrixXi numTries_trees_allChains_local;

  Eigen::MatrixXd numAccpt_trees_allChains_global;
  Eigen::MatrixXi numTries_trees_allChains_global;

  std::vector<std::string> local_trees;
  std::vector<std::string> global_trees;
  Eigen::MatrixXd local_mcmcStates;
  Eigen::MatrixXd global_mcmcStates;

  Eigen::MatrixXd local_totalLogLik;
  Eigen::MatrixXd global_totalLogLik;
  
public:
	
  void InitializeMCMC(unsigned int crrProcID, IM im, unsigned int nProcs);
  void InitializeMCMC_usePriorOnly(unsigned int n_process_id, IM im, unsigned int nProcs);
  
  void deleteHotChains();
  void deleteTrees();
  
  void updateTreeAcceptanceRate();
  void UpdateAllChains(unsigned int id_crrIter, vector<locus> lc, unsigned int save, unsigned int swapper, unsigned int swappe);
  void UpdateAllChains(int n_current_iteration, vector<locus> lc);
  void UpdateAllChains_Burning(vector<locus> lc);

  // saving step in MCMC
  void collectUpdates(unsigned int savingID);
  void collectUpdates_lociInParallel(unsigned int savingID);
  void collectUpdates_chainsInParallel(unsigned int savingID);
  void collectUpdates_onSameProc(unsigned int savingID);
  void collectUpdates_fromColdChain_onSameProc(unsigned int savingID, unsigned int coldChainID);
  void collectUpdates_fromColdChain_onDiffProcs(unsigned int savingID, unsigned int procID_wColdChain,unsigned int coldChainID);
  void collectUpdates_bwProcs_lociInParallel(unsigned int savingID);
  void collectUpdates_logPrior_lociInParallel(unsigned int savingID);
  void collectUpdates_logLikelihood_lociInParallel(unsigned int savingID);

  void SetNChains(unsigned int n) { n_chains = n; return;}
  void setZero_numTriesTrees_chain(unsigned int chainID){chains.at(chainID).setZero_numTries_trees(); return;}
  void setZero_numAccptTrees_chain(unsigned int chainID){chains.at(chainID).setZero_numAccpt_trees(); return;}

  void saveCoalTimes2File_chain0(){chains.at(0).saveCoalTimes2File();return;}
  void saveKappa2File_chain0(){chains.at(0).saveKappa2File();return;}
  void saveMutationScalers2File_chain0(){chains.at(0).saveMutationScalers2File();return;}
  void saveTipIDs_chain0(){chains.at(0).saveTipIDs();return;}
  void saveTreeIDs_chain0(){chains.at(0).saveTreeIDs();return;}
  void saveListTopo2File_chain0(){chains.at(0).saveListTopo2File();return;}
  void saveLogPriorTrees_chain0(){chains.at(0).saveLogPriorTrees();return;}
  void saveLogJointLikelihoods_chain0(){chains.at(0).saveLogJointLikelihoods();return;}
  void saveLogEachLikelihoods_chain0(){chains.at(0).saveLogEachLikelihoods();return;}
  void saveEachLogPrior_chain0(){chains.at(0).saveEachLogPrior();return;}
  void saveRNormalDist_chain0(){chains.at(0).saveRNormalDist(); return;} // for debugging only

  void print_allLogLik();
  void print_coldChainState();
  void print_globalMCMCstate();
  void printMCMCstate_eachChain(unsigned int chainID, unsigned int crrIter);
  void PrintMCMCState();
  void PrintMCMCState(int chain_ID); // Arun's old version
  void printAllChains(unsigned int id_iter);
  void print_mcmc_initial();
  void print_savedTrees(unsigned int upto_IterID);
  void print_AccepRate_swapChains(){std::cout <<AccepRate_swapChains <<"\n";return;}
  void print_treeAccptRate_chain(unsigned int chainID){chains.at(chainID).print_treeAccptRate();return;}
  void print_globalTreeAcceptanceRate();
  void print_AccepRate_swapChains_global();
  
  void fileprintAllChains(unsigned int id_iter, std::ofstream& fp);
  void FilePrintMCMCState(std::ofstream& fp1, std::ofstream& fp2, int swapping);
  
  //  unsigned int get_coldChainID(){return coldChainID;}
  int GetProcessID() { return process_id; }
  double GetChainTemperature(unsigned int chain_ID);
  double GetChainPosterior(int chain_ID);
  //AS: changed this - for details see function
  Chain getColdChain();
  int getColdChain(int process_id);
  Chain getChain(int chain_id) { return chains.at(chain_id); }
  unsigned int getRankTempOfChain(unsigned int chainID){return chains.at(chainID).getRankTemp();}
  Eigen::MatrixXi get_numTriesTrees_chain(unsigned int chainID){return chains.at(chainID).get_numTries_trees();}
  Eigen::MatrixXd get_numAccptTrees_chain(unsigned int chainID){return chains.at(chainID).get_numAccpt_trees();}
  
  void computeGlobalTreeAcceptanceRate();
  void compute_localAcceptanceRate_chainSwap(unsigned int rankTemp1, unsigned int rankTemp2, unsigned int swapSuccess);
  void compute_globalAcceptanceRate_chainSwap();

  void collect_localCurrentMCMCstates();


  void swapChains(unsigned int swapper_globalChainID, unsigned int swappee_globalChainID);
  void ChangeMCMCState(int chain_ID, double new_temperature);
  int SwapChainHeats(int n_chains);
  void SwapChainHeats(int chain1, int chain2);
  int SwapChainHeats_WithinProcess(int chain1, int chain2);
  int SwapChainHeats_BwProcesses(int swapper, int swappee, unsigned int n_chain, int currentid); // Yujin's new version - updated Arun's version --- YC 8/13/2014
  int SwapChainHeats_BwProcesses(int *swapper, int *swappee, unsigned int n_chains, int step, int currentid, int swapA, int swapB); // Arun's verion. I simplified this function. See another "SwapChainHeats_BwProcesses" --- YC 8/13/2014
  int AttemptSwap(int chain1, int chain2);
  double CalcSwapProbability(int chain1, int chain2);
  double GetMHSwapTerm(int chain_ID);
  double GetUpperKappa() {return upperkappa;}
  
  void SetUpperKappa(double upper_kappa) { upperkappa = upper_kappa;}
  //  void set_coldChainID(double newID){coldChainID = newID; return;}
  
  //AS: adding this for sharing coalescent tree info within the same node
  void sharecoaltree(unsigned int id_crrIter, unsigned int loc);
  
  void addSwapChainsTry(unsigned int chainID1, unsigned int chainID2, unsigned int successSwap);
  void update_orderGlobalChainIDs(unsigned int chainID1, unsigned int chainID2);
  unsigned int get_orderGlobalChainIDs(unsigned int chainID){return order_globalChainIDs.at(chainID);}

	~MCMC() {};
};

void saveKappa2File(std::vector<std::vector<double> > kappa);
void saveMutationScalers2File(std::vector<std::vector<double> > mutationScaler);
	
#endif
