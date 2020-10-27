/* MIST copyright 2016 by Yujin Chung and Jody Hey */


#ifndef POPTREE_HPP_
#define POPTREE_HPP_


#include <vector>
#include <string>
#include "Eigen/Dense"
#include "IM.hpp"

// 2020-10-27 YC creating a node for population tree
class popNode
{
  private:
  
  unsigned int isRoot; // 1 if root node; 0 otherwise
  unsigned int isTip; // 1 if tip; 0 otherwise
  unsigned int popID; 
	popNode *par;
  popNode *desc[2];
	double age;
	double populationSize;

public:
  
  popNode();
  unsigned int get_isRoot(){return isRoot;}
  unsigned int get_isTip(){return isTip;}
  unsigned int get_popID(){return popID;}
  
};




class Migration
{
private:
  
  // 2020-10-27 YC
  // 'popTree' structure has been modified as new class "popNode" was introduced
  unsigned int fromPopID; // the source popID (backward in time)
  std::vector<unsigned int> toPopIDs; // the list of recipients (backward in time)
  std::vector<double> migRates; // the list of migration rates.
  
  //---------------------------------//
	unsigned int popID; /// migration to the population with this popID. Backward in time.
	double migRate;

public:
	void initialization(unsigned int ID, double rate);

	unsigned int get_popID(){return popID;}
	double get_migRate(){return migRate;}

	void print();
};

class popTree
{
private:

  // 2020-10-27 YC
  // 'popTree' structure has been modified as new class "popNode" was introduced
  unsigned int nPopsAtTips; // the number of populations at present
  
  popNode *rootNode; // the root node of a population tree.
  std::vector<popNode*> popNodeList; // the list of nodes of a population tree in the order of their population IDs (starting from 1).

  unsigned int nEpochs; // the number of epochs (defined by the splitting times)
  std::vector<double> splittingTimes; // the list of splitting times in ascending order (smallest to largest). The length of splittingTimes is the same as nEpochs.
  std::vector< std::vector<unsigned int>> popIDs_eachEpoch; // the list of population ids of each epoch
  std::vector< std::vector<Migraion*>> mig_eachEpoch; // the list of migrations of each epoch.
  
  
  

  

  //--------------------------------------------//
  
	unsigned int isRoot;
	unsigned int isTip;
	unsigned int popID;
	popTree *par;
	popTree *desc[2];
	double age;
	double populationSize;

	// Added by YC 5/9/2014
	std::vector<popTree*> pop2mig;
	std::vector<double> migRate;

	// old version
	std::vector<Migration> mig;

  unsigned int no_assignedChildren; // 0 if no children assigned (desc[2] is empty); 1 if desc[0] is assigned; 2 if both are assigned.

  unsigned int ancPop; 
  // 0 if island model (no ancestral population);
  // 1 if isolation model
  unsigned int samePopulationSizes;
  // 1 if the same population sizes; 0 if different population sizes
  unsigned int sameMigrationRates;
  // 1 if the same migration rates; 0 if different population sizes

  
  
	void print_node();

	std::vector<int> find_popIDs2migrate();
	std::vector<int> find_popIDs2migrate_moveDown(std::vector<int> popIDs2migrate,double lowerB, double upperB);
	std::vector<int> find_popIDs2migrate_moveUpDown(std::vector<int> popIDs2migrate,double lowerB, double upperB);
	void find_popSizeTips(Eigen::VectorXd &popsizes);
	double find_migRate(unsigned int popID_fr, unsigned int popID_to);


public:

  
  popTree();
  //~popTree();
  void deletePopTree();

  void initialization(IM im);
  void initialize_popTree(IM im, unsigned int processID);
  void initialize_popTree_recursion(IM im, std::string newickTree, Eigen::Vector3d paraMax);
  void initialize_migrations_recursion(IM im,std::vector<popTree*> pops, double rateMax);
  void initialize_migrations(IM im,double rateMax, unsigned int processID);


  void assign_isRoot_isTip(unsigned int root, unsigned int tip) {isRoot = root; isTip = tip;}
  void assign_popID (unsigned int id){popID = id;}
  void assign_age(double t){age = t;}
  void assign_populationSize(double size) {populationSize = size;}
  void assign_desc(unsigned int i, popTree* tr){desc[i] = tr;}
  void assign_par(popTree* tr){par = tr;}
  void set_isTip(unsigned int tip){isTip = tip; return;}
  
  // Added by YC 5/15/2014
  void replace_migRate(std::vector<double> newMigRate){migRate = newMigRate;}
  void replace_pop2mig(std::vector<popTree*> pops){pop2mig = pops; return;}
  
  // Added by YC 5/9/2014
  void add_pop2mig(popTree* pop){pop2mig.push_back(pop);}
  void add_migRate(double rate){migRate.push_back(rate);}

  unsigned int get_popID(){return popID;}
  unsigned int get_isRoot(){return isRoot;}
  double get_popSize(){return populationSize;}
  double get_age(){return age;}
  popTree* get_parentNode(){return par;}
  popTree* get_pop2mig(unsigned int i){return pop2mig.at(i);}
  double get_migRate(unsigned int i){return migRate.at(i);}
  unsigned int get_no_assignedChildren(){return no_assignedChildren;}
  unsigned int get_ancPop(){return ancPop;} 
  popTree* get_firstChild(){return desc[0];}
  popTree* get_secondChild(){return desc[1];}

  void print_poptree();
  void print_popSize();
  void print_migRate();
  void print_allPopTree();

 
  popTree* deepCopy_root();
  
  popTree* go2root();
  unsigned int size();
  double findPopSize(int popid, int &found);
  double get_minSplittingTime(double time);
  int getID_beforeSplit(int id,double SplittingT, int &found);
  
  Eigen::MatrixXd getTransitionRateMat(Eigen::MatrixXi space, Eigen::VectorXi totalN);
  
  void truncate_updateIDs(double time, int eventID);
  
  double computeJointPrior(Eigen::Vector3d paraMax);

 
  double find_migrationRate(unsigned int popID_from, unsigned int popID_to);
  double find_popSize(unsigned int pop_ID);
  
  // void replacePara(unsigned int ancPop,Eigen::MatrixXd listPara);
  void replacePara(Eigen::MatrixXd listPara);
};







#endif /* POPTREE_HPP_ */
