/* MIST copyright 2016 by Yujin Chung and Jody Hey */


#ifndef _coaltree_
 #define _coaltree_


#include <vector>
#include <valarray>
#include <list>
#include <memory>
#include "Eigen/Dense"
#include "IM.hpp"
#include "misc.hpp"
#include "popTree.hpp"

using namespace std;



/// Defining
class node
{
private:

  unsigned int isRoot; // 1 if the node is a root
  unsigned int label; // never used in MCMC. Used in L mode when a conditional distribution
  // of a coalescent tree given a population tree is computed.
  // labels start from 0.
  double totalCoalRate;
  unsigned int siblingOrder; // 0 if this node is the first child; 1 if the second child; 2 for initial value before a tree is constructed or for root node
  int n_mut;     /// The number of mutations on the edge
  
  // Computing the likelihood for the node (subtree) including its parent node with an edge to it.
  Eigen::MatrixXd lik_oneside_HKY(locus lc, double mutrate, double kappa);
  node* randomPick_node(); // Randomly pick a node to detach (used in a proposal of a new tree in MCMC)
  node* search_node2detach(int nodeID);
	node* search_node2detach_fromRoot(int nodeID);
	node* get_node(); // to copy the current node

  void print_node();
  void print_numMutations();

	void slider(double dist); // slide the current node (used when a new tree is proposed in MCMC)
	int find_sisterID();
	int isSameTree(node* node); //
	node* go2root();
	// vector<double> get_coalescentTimes(vector<double> coalTimes);

  std::vector<unsigned int> cumMutations;
  std::vector< std::vector<unsigned int>> cumulativeMutations;
  std::vector<double> waitingTime; // for likelihood calculation
  double numMutations; // the number of mutations (for a given site) of the ancestral edge of this node. used for likelihood computation under infinite-site model.
  double totalLengthsLineages; // 0 if it is undefined
  double numMutations_overLocus; // the number of mutations among the given locus occurred on the ancetral edge of this node.
  unsigned int mutationOnThisBranch; // 1 if a mutation happened on this branch; 0 otherwise.
  unsigned int rejectThisTree_IS; // 1 if this tree is incompatible; 0 otherwise

public:
	node(); // constructor
	void initialization();
	//node(node* tree); // deep copy constructor
	// node(const node& newNode); // copy constructor
	// virtual ~node();
	//Matrix	lik; // 4 by (# unique sites) matrix
	Eigen::MatrixXd lik; // 4 by (# unique sites) matrix
	unsigned int isTip; // 1 if the node is a tip
	unsigned int isLikelihoodNULL; // 1 if 'lik' has not been computed or initialized; 0 if 'lik' is computed.
  unsigned int tipID; // starting from 1.
	node *par; //Pointer to the parent node
	                  // We may not need this parent node
	node *desc[2]; /// Pointers to the daughter nodes
	double age; /// Time to the node from the tip
	            // Note that the branch length from the left daughter to the node
	            // is age - desc[0].age.
	unsigned int popID; /// The population ID including the node
	           /// connecting the edge and its daughters is.
			   /// This popID is used in L mode only, when
	           /// Pr(coalescent tree|data) is computed.
	// node* deepCopy_par();

  void initial_forComputingNumMutationOverBranches(locus lc, unsigned int site); // for IS
  void initial_forComputingNumMutationOverBranches_new(locus lc, unsigned int site); // for IS
  void initial_forComputingCumulativeNumMutationOverBranches(locus lc); // for IS
  void initial_forComputingCumulativeNumMutationOverBranches(locus lc, unsigned int site);


  void deleteCoalTree();
  void assign_siblingOrder(unsigned int i) {siblingOrder = i;}
  
  unsigned int get_siblingOrder(){return siblingOrder;}
  Eigen::MatrixXd get_likMatrix();
  double get_totalCoalescentRate() { return totalCoalRate; }
  double get_totalLengthsLineages(){return totalLengthsLineages;}
  std::vector<double> get_waitingTime(){return waitingTime;}
  unsigned int get_isTip(){return isTip;}
  unsigned int get_isRoot(){return isRoot;}
  void set_isRoot(unsigned int ir){isRoot=ir;}
  unsigned int get_tipID(){return tipID;}
  double get_numMutations(){return numMutations;}
  double get_numMutations_overLocus(){return numMutations_overLocus;}
  double get_age(){return age;}
  std::vector<unsigned int> get_cumMutations(){return cumMutations;}
  std::vector<unsigned int> get_cumulativeMutations(unsigned int s){return cumulativeMutations.at(s);}
  unsigned int get_rejectThisTree_IS(){return rejectThisTree_IS;}

  void set_numMutations(double x){numMutations = x; return;}
  void set_rejectThisTree_IS(unsigned int x){rejectThisTree_IS=x; return;}
  void setZero_totalLengthsLineages(){totalLengthsLineages=0;return;}

  void print_coaltree();
  void print_nodeMembers();
  void print_labelsPopIDs();

  void saveTree(ofstream& filename);
  void saveNode(ofstream& filename);


  std::string convert2string_newick();
  std::string convert2string_newick_root();
  void convertFromNewick(string tree, unsigned int root);


  unsigned int size_tree();
  void maxTipID(int &maxID);
  node* findNode_fromRoot(int tipid);
  node* findNode(int tipid, int &found);
  double get_minCoalTime(double coalT);
  list<double> get_coalescentTimes(list<double> coalTimes);
  // list<double> get_coalescentTimes_rec(list<double> coalTimes);
  Eigen::VectorXi get_crrState(int noPops);
  void get_tipState(Eigen::VectorXi &state, int nPops);
  
  node* deepCopy();
  node* deepCopy_root();
  void deepCopy(node* tr);
  void deepCopy_root(node* tr);
  void deepCopy_obj(node* tr);
  
  void truncate(int &numNewTips, double time);
  void truncateFirstCoalNode(int newtipID);
  void updatePopIDs(popTree* poptree, double eventTime);
  node* get_truncatedTree();
  void newLabels();
  
  int moreNextState(int nPops);
  void makeNULL_ancestorsLik();
  void assignPopulations2Tips(locus lc);
  void assignPopulations2Tips(std::vector<unsigned int> SeqPop);
  int assignPopulation(Eigen::VectorXi popAssign, int &idx, Eigen::VectorXi &state, int nPops);
  Eigen::VectorXi get_totalNumEachKind(int noPops);
  void compute_totalCoalescentRate();
  void compute_totalCoalescentRate(std::list<double> coalTimes);
  void compute_likMatrix_HKY(locus lc, double mutrate, double kappa);// Computing the likelihood
  node* propose_coaltree(double slidedist); // propose a new tree from the current tree
  void rescale_treeHeight(double treeHeightScaler);
  double conditionalPr_singlePop(double popSize);
  double logLikelihood_HKY (locus lc, double mutrate, double kappa); // returning log-likelihood under HKY model
  double logLikelihood_IS (locus lc, double mutrate, unsigned int updateMut);
  double compute_prior_indepLoci(locus lc, double popSizeMax, unsigned int priorType, double changeOfRatePoint);

  void computeTotalLengthsOfLineages(); // for IS model
  void computeNumMutationsOverBranches_site(); // for IS model
  void computeNumMutationsOverBranches_site_new(); // for IS model
  // double computeNumMutationsOverBranches(locus lc); // for IS model
  double computeWaitingTime_numMutationsOverLocus();
  double computeWaitingTime_new();
  double compute_partitionFunctionOfPossionDist_numMutations();
  unsigned int purity_cumMutations();
  unsigned int purity_cumulativeMutations(unsigned int s);
  void computingCumulativeNumMutations_waitingTime(locus lc);
  double computeWaitingTime_numMutationsOverLocus_new();
  // double computingWaitingTime_numMutationsOverLocus(locus lc);
  void computingWaitingTime_numMutationsOverLocus(unsigned int site);
  double computingWaitingTime_numMutationsOverLocus_recursive(unsigned int char2purify,unsigned int s);
  double computeWaitingTime_numMutationsOverLocus_recursive(unsigned int char2purify);
  void computeCumNumMutationsOverBranches();
  void computeCumulativeNumMutationsOverBranches(locus lc);
  void computeCumulativeNumMutationsOverBranches(unsigned int site);

  void increase_numMutationsOverLocus_byOne(){numMutations_overLocus++; return;};

  int validTree();
  //void get_UPGMA (locus lc);
#ifdef MPI_ENABLED	
  void MPIsend_coaltimes_tipIDs();
  void send_size();
  void MPIsend_coaltree();
  void send_tipids();
  unsigned int receive_size();
  void send_topoinfo();
#endif

  // subtree
  std::vector<node*> getSubtree(unsigned int subSize);
};


/// Defining
class nodeSimple
{
private:
  unsigned int isRoot; // 1 if the node is a root
  unsigned int isTip; // 1 if the node is a tip
  nodeSimple *par; //Pointer to the parent node
                   // We may not need this parent node
  nodeSimple *firstChild; /// Pointers to the daughter nodes
  nodeSimple *secondChild; /// Pointers to the daughter nodes
  unsigned int popID; /// The population ID including the node
		           /// connecting the edge and its daughters is.
				   /// This popID is used in L mode only, when
		           /// Pr(coalescent tree|data) is computed.
  unsigned int rank; // larger for older node.
  /* If tip, rank =0.
   */
  unsigned int size;

  // Computed the value and use in Chain::compute_summaryOfCoalTrees() 
  unsigned int AreTipsFromSamePop; 
  unsigned int nodePopID; // the population ID from which the (descent) tips are
                          // If tips are from different populations, nodePopID =0
  double age;
  double popSize;
  /* If tip, popSize = 0
   */
  int nodeLabel;
  /* If tip, nodeLabel = popID.
   */

  unsigned int totalNumSeq; // the total number of sequences or tips
  /* If a full tree is considered, this number is same as the size of tree.
     If a subtree is considered, this number is larger than the size of tree.
   */


public:
  void assignSize(unsigned int s){size = s;}
  
  unsigned int get_isTip(){return isTip;}
  unsigned int getSize(){return size;}
  unsigned int getPopID(){return popID;}
  unsigned int getRank(){return rank;}
  nodeSimple* getFirstChild(){return firstChild;}
  nodeSimple* getSecondChild(){return secondChild;}
  nodeSimple* getPar(){return par;}
  unsigned int get_nodePopID(){return nodePopID;}
  unsigned int get_AreTipsFromSamePop(){return AreTipsFromSamePop;}
  double get_age(){return age;}
  int get_nodeLabel(){return nodeLabel;}
  // int get_nodeLabel(){return nodeLabel;}
  
  void convert_oldversion(node* tree);
  void convert(node* tree, std::list<double> coaltimes, unsigned int nLineages); // updated by YC 5/8/2014
  unsigned int sameTopo(node* tree);
  unsigned int sameTopo(nodeSimple* topo);
	void computeSizes();
  std::string convert2Newick_recursion(std::list<double> coalT, std::vector<unsigned int> tipIDs, double parentAge);
  std::string convert2Newick(std::list<double> coalT, std::vector<unsigned int> tipIDs);
  std::string convert2Newick_originalTopo(std::vector<unsigned int> tipIDs);
  std::string convert2Newick_originalTopo_recursion(std::vector<unsigned int> tipIDs);
  std::string convert2Newick_topo();
  std::string convert2Newick_topo_root();
  void convertFromNewick(string topo, unsigned int root);

	void deleteTopo();


	void print_topo();
	void fileprint_topo(std::ofstream& fp);
	void print_nodeMembers();
	void print_labelsPopIDs();
	void saveTree(ofstream& filename);
	void saveNode(ofstream& filename);
	void tree2string();

	int size_tree();
	void maxTipID(int &maxID);
	node* findNode_fromRoot(int tipid);
	node* findNode(int tipid, int &found);
	double get_minCoalTime(double coalT);
	list<double> get_coalescentTimes(list<double> coalTimes);
	Eigen::VectorXi get_crrState(int noPops);
	void get_tipState(Eigen::VectorXi &state, int nPops);

	//auto_ptr<node> deepCopy_root();
	node* deepCopy();
	node* deepCopy_root();
	void deepCopy(node* copiedTree, node* tr);
	void deepCopy_root(node* copiedTree, node* tr);
	void deepCopy_obj(node copiedTree, node* tr);

	void truncate(int &numNewTips, double time);
	void truncateFirstCoalNode(int newtipID);
	void updatePopIDs(popTree* poptree, double eventTime);
	node* get_truncatedTree();
	void newLabels();

	int moreNextState(int nPops);
	void makeNULL_ancestorsLik();
	void assignPopulations2Tips(locus lc);
  //  void assignNodeLabel();
	int assignPopulation(Eigen::VectorXi popAssign, int &idx, Eigen::VectorXi &state, int nPops);
  void assign_age_nodeLabel_popSize(std::vector<double> coalT, std::vector<double> allPopSize, double splittingTime);
	Eigen::VectorXi get_totalNumEachKind(int noPops);
	void compute_totalCoalescentRate();
	void compute_likMatrix_HKY(locus lc, double mutrate, double kappa);// Computing the likelihood
	node* propose_coaltree(double slidedist); // propose a new tree from the current tree
	double conditionalPr_singlePop(double popSize);
	double logLikelihood_HKY (locus lc, double mutrate, double kappa);
	int validTree();


	#ifdef MPI_ENABLED
	void MPIreceive_coaltree(std::list<double> ct, int coldprocess);
	unsigned int receive_topoinfo(int coldprocess);
	unsigned int receive_size(int coldprocess);
	void send_size(int coldprocess);
	void MPIreceive_tipids(int coldprocess);
	#endif

  // used in L mode
  void compute_nSamples_nSymNode(Eigen::MatrixXi &nSamples, unsigned int &nSymNode);
  unsigned int compute_nTrees_sameTopo(unsigned int nPops);
  void compute_areTipsFromSamePop();
  std::vector<unsigned int> get_nodePopID_withRank(unsigned int R);

  long double compute_logProb_zeroMig(std::vector<double> coalT, std::vector<double> allPopSize, double splittingTime);
  unsigned int compute_nCoalEvents(int label);
  unsigned int compute_maxNumLin(int label);
  std::list<double> compute_intervalCoalT(int label, double splittingTime, std::list<double> intervalCoalT);

  void MPIreceive_coaltree(int senderID);
  void MPIsend_nodeSimple(unsigned int receiverID);
  

};



/// Defining
class node_old
{
private:

	unsigned int label; // never used in MCMC. Used in L mode when a conditional distribution
					// of a coalescent tree given a population tree is computed.
					// labels start from 0.
	double totalCoalRate;
	unsigned int siblingOrder; // 0 if this node is the first child; 1 if the second child; 2 for initial value before a tree is constructed or for root node


	// Computing the likelihood for the node (subtree) including its parent node with an edge to it.
	Matrix lik_oneside_HKY(locus lc, double mutrate, double kappa);
	node* randomPick_node(); // Randomly pick a node to detach (used in a proposal of a new tree in MCMC)
	node* search_node2detach(int nodeID);
	node* search_node2detach_fromRoot(int nodeID);
	node* get_node(); // to copy the current node
	void print_node();
	void slider(double dist); // slide the current node (used when a new tree is proposed in MCMC)
	int find_sisterID();
	int isSameTree(node* node); //
	node* go2root();
	// vector<double> get_coalescentTimes(vector<double> coalTimes);

public:
	node_old(); // constructor
	void initialization();
	//node(node* tree); // deep copy constructor
	// node(const node& newNode); // copy constructor
	// virtual ~node();
	Matrix lik; // 4 by (# unique sites) matrix
	int isRoot; // 1 if the node is a root
	int isTip; // 1 if the node is a tip
	unsigned int isLikelihoodNULL; // 1 if 'lik' has not been computed or initialized; 0 if 'lik' is computed.
	int tipID; // starting from 1.
	node_old *par; //Pointer to the parent node
	                  // We may not need this parent node
	node_old *desc[2]; /// Pointers to the daughter nodes
	double age; /// Time to the node from the tip
	            // Note that the branch length from the left daughter to the node
	            // is age - desc[0].age.
	int n_mut;     /// The number of mutations on the edge
	int popID; /// The population ID including the node
	           /// connecting the edge and its daughters is.
			   /// This popID is used in L mode only, when
	           /// Pr(coalescent tree|data) is computed.
	// node* deepCopy_par();
	void deleteCoalTree();
	void assign_siblingOrder(unsigned int i) {siblingOrder = i;}

	unsigned int get_siblingOrder(){return siblingOrder;}
	Matrix get_likMatrix();
	double get_totalCoalescentRate() { return totalCoalRate; }

	void print_coaltree();
	void print_nodeMembers();
	void print_labelsPopIDs();
	void saveTree(ofstream& filename);
	void saveNode(ofstream& filename);
	void tree2string();

	int size_tree();
	void maxTipID(int &maxID);
	node* findNode_fromRoot(int tipid);
	node* findNode(int tipid, int &found);
	double get_minCoalTime(double coalT);
	list<double> get_coalescentTimes(list<double> coalTimes);
	Eigen::VectorXi get_crrState(int noPops);
	void get_tipState(Eigen::VectorXi &state, int nPops);

	//auto_ptr<node> deepCopy_root();
	node* deepCopy();
	node_old* deepCopy_root();
	void deepCopy(node* copiedTree, node* tr);
	void deepCopy_root(node_old* copiedTree, node_old* tr);
	void deepCopy_obj(node_old copiedTree, node_old* tr);

	void truncate(int &numNewTips, double time);
	void truncateFirstCoalNode(int newtipID);
	void updatePopIDs(popTree* poptree, double eventTime);
	node* get_truncatedTree();
	void newLabels();

	int moreNextState(int nPops);
	void makeNULL_ancestorsLik();
	void assignPopulations2Tips(locus lc);
	int assignPopulation(Eigen::VectorXi popAssign, int &idx, Eigen::VectorXi &state, int nPops);
	Eigen::VectorXi get_totalNumEachKind(int noPops);
	void compute_totalCoalescentRate();
	void compute_likMatrix_HKY(locus lc, double mutrate, double kappa);// Computing the likelihood
	node* propose_coaltree(double slidedist); // propose a new tree from the current tree
	double conditionalPr_singlePop(double popSize);
	double logLikelihood_HKY (locus lc, double mutrate, double kappa);
	int validTree();
	void get_UPGMA (locus lc);
};



//typedef vector<node*> trees_loci; // coalescent trees of loci

void print_coaltree(node *root);
int size_node(node *nd);


// Fixme: Delete it later
//int * get_loc_minMat (double **mat, int size);
// void get_UPGMA (locus* lc);
// vector<node*> initial_coalTree_HKY(vector<locus> loci);
//void initial_coalTree_HKY(vector<node*> tr, vector<locus> loci);
node* initial_coalTree_IS(locus lc);
node* initial_coalTree_IS_fromIMa2(locus lc);
node* initial_coalTree_HKY(locus lc);


void saveTrees2File(std::vector<std::vector<node* > > trees);
void saveCoalTimes2File(std::vector<std::vector<node* > > trees);


#endif
