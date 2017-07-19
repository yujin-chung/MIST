/* MIST copyright 2016 by Yujin Chung and Jody Hey */


#ifndef POPTREE_HPP_
#define POPTREE_HPP_


#include <vector>
#include <string>
#include "Eigen/Dense"
#include "IM.hpp"

class Migration
{
private:
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
  unsigned int isRoot;
  unsigned int isTip;
  unsigned int popID;
  popTree *par;
  std::vector<popTree*> desc;
  // popTree *desc[2];  // old version
  double age; // backward in time
  double populationSize;
  
	// Added by YC 5/9/2014
  std::vector<popTree*> pop2mig;
  std::vector<double> migRate;
  std::vector<std::vector<double> > epoch2mig;

  // Added by YC 7/19/2017
  /** splittingTimes.at(i) and pops.at(i) contains the age of popID i+1 and the pointer to popID i+1, respectively.
   * only root node has this information.
   */
  std::vector<double> splittingTimes;
  std::vector<popTree* > pops;
  
  
  // old version
  std::vector<Migration> mig;
  
  unsigned int no_assignedChildren; // 0 if no children assigned (desc[2] is empty); 1 if desc[0] is assigned; 2 if both are assigned. It can be larger than 2.

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
  void initialize_popTree_recursion(IM im, std::string newickTree, Eigen::Vector3d paraMax, popTree* root);
  // void initialize_migrations_recursion(IM im,std::vector<popTree*> list_pops, double rateMax);
  //  void initialize_migrations_recursion(IM im, popTree* pop, double rateMax);
  void initialize_migrations(IM im,double rateMax, unsigned int processID);

  
  void assign_isRoot_isTip(unsigned int root, unsigned int tip) {isRoot = root; isTip = tip;}
  void assign_popID (unsigned int id){popID = id;}
  void assign_age(double t){age = t;}
  void assign_populationSize(double size) {populationSize = size;}
  void assign_desc(unsigned int i, popTree* tr){desc.at(i) = tr;}
  void assign_par(popTree* tr){par = tr;}
  void set_isTip(unsigned int tip){isTip = tip; return;}
  
  // Added by YC 5/15/2014
  void replace_migRate(std::vector<double> newMigRate){migRate = newMigRate;}
  void replace_pop2mig(std::vector<popTree*> pops){pop2mig = pops; return;}
  
  // Added by YC 5/9/2014
  void add_pop2mig(popTree* pop){pop2mig.push_back(pop);}
  void add_migRate(double rate){migRate.push_back(rate);}
  void add_epoch2mig(std::vector<double> epoch){epoch2mig.push_back(epoch);}
  void add_splittingTimes(unsigned int loc, double age){splittingTimes.at(loc) = age; return;}
  void add_pops(unsigned int loc, popTree* pop){pops.at(loc) = pop; return;}
  
  unsigned int get_popID(){return popID;}
  unsigned int get_isRoot(){return isRoot;}
  double get_popSize(){return populationSize;}
  double get_age(){return age;}
  popTree* get_parentNode(){return par;}
  popTree* get_pop2mig(unsigned int i){return pop2mig.at(i);}
  double get_migRate(unsigned int i){return migRate.at(i);}
  unsigned int get_no_assignedChildren(){return no_assignedChildren;}
  unsigned int get_ancPop(){return ancPop;}
  popTree* get_child(unsigned int i){return desc.at(i);}
  popTree* get_firstChild(){return desc.at(0);}
  popTree* get_secondChild(){return desc.at(1);}
  int get_size_of_splittingTimes(){return splittingTimes.size();}
  int get_size_of_pops(){return pops.size();}
  int get_size_of_migRate(){return migRate.size();}
  double get_epoch2mig(unsigned int i, unsigned int j){return epoch2mig.at(i).at(j);}

  void resize_splittingTimes(unsigned int newSize){splittingTimes.resize(newSize); return;}
  void resize_pops(unsigned int newSize){pops.resize(newSize); return;}

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
