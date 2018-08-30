/*MIST copyright 2016 Yujin Chung and Jody Hey */

/*
 * IM.hpp
 *
 */

#ifndef IM_HPP_
#define IM_HPP_

///AS: adding a global define to indicate MPI
#define MPI_ENABLED

#include <vector>
#include <valarray>
#include <string>
#include "Eigen/Dense"
// #include "popTree.hpp"
#include "misc.hpp"



#ifdef MPI_ENABLED
#include <mpi.h>
#endif

#include <fstream>

// using namespace std;

/***
 * for one locus. For multiple loci (or general case),
 * "class IM" has a member of 'vector<locus> loci'.
 */
class locus
{
  vector<vector<unsigned int> > seq_uniq; // (n_sites_uniq) by (n_genes) matrix
  vector<unsigned int> freq_uniqueSeq; // size: n_sites_uniq
  unsigned int n_sites_uniq; // the number of unique sites
  unsigned int nSNPs; // the number of SNPs (used only for IS)
  valarray<double> pi; // empirical stationary distribution of (A,C,G,T)
  unsigned int modelType; // mutation/substitution type (0: IS, 1: JC, 2: HKY)

  int n_sites;
  unsigned int n_geneCopies; /// The number of gene copies
  unsigned int done_computingSumOfMuationsOnBranches;
  unsigned int multiLocusSpecific_mutationRate; 
  // 0 if all the loci have the same mutation rates; 1 if mutation rate per site is constant; 2 if mutation rate per locus are various.

public:

  vector<string> popNames;
  vector<string> tipNames;

  void initializeLikelihoodModel(string model);
  void initialize_multiLocusSpecific_mutationRate(unsigned int lociRate){multiLocusSpecific_mutationRate = lociRate;return;}

  void print_data();
  void print_pi();
  void print_likelihoodModel();

  Matrix get_distMat();
  unsigned int get_uniqSeq(int idx_site, int idx_gene);
  vector<unsigned int> getSite_uniqSeq(unsigned int idx_site){return seq_uniq.at(idx_site);}
  double get_pi(int i);
  unsigned int get_n_sites_uniq();
  valarray<double> get_pi();
  vector<unsigned int> get_freq_uniqueSeq();
  double getstandfactor (double kappa);
  unsigned int getLikelihoodModel(){return modelType;}
  unsigned int get_nGeneCopies(){return n_geneCopies;} // the number of genes
  unsigned int get_nSites(){return n_sites;}
  unsigned int get_doneComputingSumOfMuationsOnBranches(){return done_computingSumOfMuationsOnBranches;}
  unsigned int get_multiLocusSpecific_mutationRate(){return multiLocusSpecific_mutationRate;}
  unsigned int get_modelType(){return modelType;}
  unsigned int get_sizeOfSeqUniq(){return seq_uniq.size();}
  unsigned int get_nSNPs(){return nSNPs;}
  vector<vector<unsigned int> > get_seq_uniq(){return seq_uniq;}

  void set_nGeneCopies(unsigned int n){n_geneCopies = n; return;}
  void set_nSites(unsigned int n){n_sites = n; return;}
  void set_n_sites_uniq(int n){n_sites_uniq = n; return;}
  void set_seq_uniq(vector<vector<unsigned int> > seq){seq_uniq = seq;return;}
  void set_freq_uniqueSeq(vector<unsigned int> freq){freq_uniqueSeq = freq;return;} 
  double transitionPr_HKY(double t, int idx_gene_fr, int idx_gene_to,double mutrate, double kappa);
  int IDsameSeq(vector<unsigned int> site);

  
  void compute_pi_uniqSeq(vector<vector<char> > data);
  void compute_pi_uniqSeq_HKY(vector<vector<char> > data);
  void compute_uniqSegSites_IS(vector<vector<char> > data);
};


/**
 * YC 2/24/2014
 *	This class currently contains the DNA sequences (called 'loci')
 *	and the number of loci (called 'n_loci') as private members,
 *	but we may use this class to contain other analysis settings
 *	such as the basic setting for MCMC (the number of generations,
 *	prior options etc) and for parallel computing (the number of
 *	processes)
 */
class IM
{
private:
  vector<locus> loci;
  unsigned int n_loci; /// the number of loci
  unsigned int n_chains; /// the number of MCMC chains
  unsigned int n_MCMCgen; /// the number of (total) MCMC generations
  unsigned int n_sampledTrees; // the number of sampled trees after burning and thinning step.
  unsigned int n_Burning;
  unsigned int thinning;
  unsigned int printFreq; // print the log-likelihood valuse of each chain 
                         // every 'printFreq' iterations and print
                         // acceptance rates every 10*'printFreq' iterations.
  double popSizeMax; // the upper bound of population size (4Nu)
  double splittingTimeMax; /// the upper bound of the splitting times
  double migRateMax; 

  unsigned int MLmodes; 
  // 1 if M mode is only called;
  // 2 if L mode is only called. Need the result files of M mode.
  // 3 if L mode is running for the given population tree only.
  // 4 for L mode for grid search. Returns posterior means.
  // 5 for MLmodes=2 but parallized over samples rather than loci

  unsigned int lociInParallel; // for M mode only
  // 1 if loci are in parallel
  // 0 if chains are in parallel in  M model; samples are in parallel for L mode.
  unsigned int nPairsMut; // the number of pairs of mutations to update in MCMC
  
  unsigned int priorType; // The kind of importance sampling distribution
  // 1 if one peak
  // 2 if two peaks
  // 3 if 3 peaks
  // 4 if 4 peaks
  // For L mode: 0 if improper prior was used in MCMC; 1 if non-uniform prior was used
  double changeOfRatePoint;
 
  unsigned int samplingFromPriorOnly; ///	1 if trees and other mutation parameters are sampled from their prior
  /// distribution, not from the posterior distribution;
  /// 0 if sampling from the posterior distribution (the standard MCMC).

  std::string poptree_string;

  //------------------------------------------//
  // Options for demographic models in L mode
  //------------------------------------------//
  unsigned int ancPop; 
  // 0 if island model (no ancestral population);
  // 1 if isolation model
  unsigned int samePopulationSizes;
  // 1 if the same population sizes; 0 if different population sizes
  unsigned int sameMigrationRates;
  // 1 if the same migration rates; 0 if different population sizes
  unsigned int nParaVectors;

  unsigned int newickTreeFormat;
  char* newickTreeFileName;
  char* SeqPopFileName;

  unsigned int multiLocusSpecific_mutationRate; // 0 if same mutation rates across loci; 1 if multilucs-specific mutation rate
  
  //----- Settings for L mode -----//

  unsigned int nTrees2compute; //the number of sampled trees for each locus to compute.
  unsigned int TreesWithMaxP;
  // 1 if trees with highest probabilities are selected 
  // 0 if some first trees are selected.


  unsigned int howOften_checkpoint;
  // how often to write a checkpoint.
  unsigned int checkpoint;
  // 0 for no checkpoint (default)
  // 1 for writing checkpoint
  // 2 for reading the checkpoint and starting from it but not writing checkpoint anymore
  // 3 for reading the checkpoint and writing checkpoint.


  unsigned int gridsize;
  unsigned int nJointPoints;
  Eigen::MatrixXd truePara;

public:
	//Loci(int nloci);
  ~IM() {};
  unsigned int initialization (int argc, char *argv[], unsigned int processID);
  // If the returning value is 0, then stop the software. If 1, execute MCMC.
  void print_data_pi();
  void initialize_nloci(int nloci); // this function can be combined into "void initialization()" later - YC 2/24/2014
  void initialize_MCMCsetting(int nChains, int nTreeSample, int samplingFromPriors); // this function can be combined into "void initialization()" later - YC 2/24/2014
  
  int get_nChains();
  int get_nMCMCgen();
  unsigned int get_nSampledTrees(){return n_sampledTrees;}
  unsigned int get_nBurning(){return n_Burning;}
  unsigned int get_thinning(){return thinning;}
  unsigned int get_printFreq(){return printFreq;}
  unsigned int get_samplingFromPriorOnly(){ return samplingFromPriorOnly;}
  unsigned int get_MLmodes(){return MLmodes;}
  unsigned int get_nPairsMut(){return nPairsMut;}
  unsigned int get_LikelihoodModelType(){return loci.at(0).get_modelType();}
  // returns the likelihood model of the first locus
  unsigned int get_multiLocusSpecific_mutationRate(){return multiLocusSpecific_mutationRate;}
  unsigned int get_newickTreeFormat(){return newickTreeFormat;}
  char* get_newickTreeFileName(){return newickTreeFileName;}
  char* get_SeqPopFileName(){return SeqPopFileName;}
  unsigned int get_nLoci() {return n_loci;}
  vector<locus> getLoci() {return loci;}
  unsigned int get_nGeneCopies() {return loci.at(0).get_nGeneCopies();}
  double get_popSizeMax() {return popSizeMax;}
  double get_splittingTimeMax(){ return splittingTimeMax;}
  double get_migRateMax(){return migRateMax;}
  Eigen::Vector3d get_paraMax();
  unsigned int get_gridsize(){return gridsize;}
  unsigned int get_nJointPoints(){return nJointPoints;}
  Eigen::MatrixXd get_truePara(){return truePara;}
  std::string get_poptree_string(){return poptree_string;}
  unsigned int get_ancPop(){return ancPop;} 
  unsigned int get_samePopulationSizes(){return samePopulationSizes;}
  unsigned int get_sameMigrationRates(){return sameMigrationRates;}
  unsigned int get_lociInParallel(){return lociInParallel;}
  unsigned int get_TreesWithMaxP(){return TreesWithMaxP;}
  unsigned int get_priorType(){return priorType;}
  double get_changeOfRatePoint(){return changeOfRatePoint;}
  unsigned int get_nTrees2compute(){return nTrees2compute;}
  unsigned int get_nParaVectors(){return nParaVectors;}
  unsigned int get_checkpoint(){return checkpoint;}
  unsigned int get_howOften_checkpoint(){return howOften_checkpoint;}
};

#endif /* IM_HPP_ */
