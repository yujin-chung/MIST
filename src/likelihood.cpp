/* MIST copyright 2016 by Yujin Chung and Jody Hey */


#include <iostream>
#include "Chain.hpp"
#include "coaltree.hpp"
#include<valarray>
using namespace std;


double locus::get_pi(int i)
{
	return pi[i];
}


/**
 * Computing the HKY transition probabilities.
 * See eq (1.20) of "Computational molecular evolution" by Z. Yang for the formula.
 */
double locus::transitionPr_HKY(double t, int base_fr, int base_to, double mutrate, double kappa)
{
	double prob = 0.0; // transition probability from base_fr to base_to
	double pi_y = pi[1]+pi[3]; // frequency of pyrimidines (C,T)
	double pi_r = pi[0]+pi[2]; // frequency of purines (A,G)
	double scaled_t = t*mutrate;

	prob =  pi[base_to];
	if(base_fr == base_to)
	{
		if(base_fr == 0 || base_fr == 2) // from base A or G
		{
			prob +=  exp(-scaled_t)*pi[base_to]*pi_y/pi_r
					+ exp(-(pi_r*kappa+pi_y)*scaled_t)*(pi_r - pi[base_to])/pi_r;
		}
		else // from base C or T
		{
			prob +=  exp(-scaled_t)*pi[base_to]*pi_r/pi_y
					+ exp(-(pi_y*kappa+pi_r)*scaled_t)*(pi_y - pi[base_to])/pi_y;
		}
	}
	else if(base_fr+base_to == 2)// transitions between A and G
		prob += exp(-scaled_t)*pi[base_to]*pi_y/pi_r - exp(-(pi_r*kappa+pi_y)*scaled_t)*pi[base_to]/pi_r;
	else if(base_fr+base_to == 4) // transitions between C and T
		prob += exp(-scaled_t)*pi[base_to]*pi_y/pi_r - exp(-(pi_r*kappa+pi_y)*scaled_t)*pi[base_to]/pi_r;
	else // transversions
		prob *=(1-exp(-scaled_t));

	return prob;
}

void node::makeNULL_ancestorsLik()
{
	isLikelihoodNULL = 1;
	if(!isRoot)
		par->makeNULL_ancestorsLik();
	return;
}


/*
 * Computing the likelihood for the subtree including the current node, its parent node with an edge to it
 * and its all descendants.
 */
Eigen::MatrixXd node::lik_oneside_HKY(locus lc, double mutrate, double kappa)
{

	Eigen::MatrixXd lik_child(4,lc.get_n_sites_uniq()), lik_br(1,4);
	lik_child.setZero();
	lik_br.setZero();
	//lik_child.initialize(4,lc.get_n_sites_uniq());
	//lik_br.initialize(1,4);
	//std::cout << "here11\n";
	// If the left child node is a tip
	if(isTip == 1)
	{
		for(unsigned int i=0; i<lc.get_n_sites_uniq(); i++)
		{
			for(unsigned int j=0; j<4; j++)
			{
				// Transition probability from the parent node (base j) to the tip for site i
				lik_child(j,i) = lc.transitionPr_HKY(par->age, j,lc.get_uniqSeq(i,tipID-1), mutrate, kappa);
				//lik_child.replace(lc.transitionPr_HKY(par->age, j,lc.get_uniqSeq(i,tipID-1), mutrate, kappa),j,i);
			}
		}
	}
	else
	{

		if(isLikelihoodNULL)
			compute_likMatrix_HKY(lc,mutrate,kappa);
		for(int j=0; j<4; j++)
		{
			for(int k=0; k<4; k++)
			{
				// Transition probability from the parent node (base j) to the current child (base k) for site i
				lik_br(0,k) = lc.transitionPr_HKY(par->age - age, j, k, mutrate, kappa);
				//lik_br.replace(lc.transitionPr_HKY(par->age - age, j, k, mutrate, kappa), 0,k);
			}
			lik_child.row(j) =lik_br * lik;
			//lik_child.replace_row(matrixMultiplication(lik_br, lik),j);
		}
	}

	return lik_child;
}


/*
 * Computing the likelihood for the subtree including the current node, its parent node with an edge to it
 * and its all descendants.
 */
Matrix node_old::lik_oneside_HKY(locus lc, double mutrate, double kappa)
{

	Matrix lik_child, lik_br;
	lik_child.initialize(4,lc.get_n_sites_uniq());
	lik_br.initialize(1,4);

	// If the left child node is a tip
	if(isTip == 1)
	{
		for(unsigned int i=0; i<lc.get_n_sites_uniq(); i++)
		{
			for(unsigned int j=0; j<4; j++)
			{
				// Transition probability from the parent node (base j) to the tip for site i
				lik_child.replace(lc.transitionPr_HKY(par->age, j,lc.get_uniqSeq(i,tipID-1), mutrate, kappa),j,i);
			}
		}
	}
	else
	{
		// REMOVE
		//  std::cout << "in lik_oneside_HKY()\n";
		//  std::cout << "\tThe current subtree is ";
		//  print_coaltree();
		//  std::cout << "\t isLikelihoodNULL: " << isLikelihoodNULL <<"\n";
		//  if(!isLikelihoodNULL)
		// 	 lik.print();

		if(isLikelihoodNULL)
			compute_likMatrix_HKY(lc,mutrate,kappa);
		for(int j=0; j<4; j++)
		{
			for(int k=0; k<4; k++)
			{
				// Transition probability from the parent node (base j) to the current child (base k) for site i
				lik_br.replace(lc.transitionPr_HKY(par->age - age, j, k, mutrate, kappa), 0,k);
			}
			lik_child.replace_row(matrixMultiplication(lik_br, lik),j);
		}
	}

	// REMOVE
	// cout << "in lik_oneside_HKY();\n";
	//		lik_child.print();

	return lik_child;
}


// Computing the likelihood for a subtree (the node and its descendants)
void node::compute_likMatrix_HKY(locus lc, double mutrate, double kappa)
{

	// REMOVE
	//std::cout << "in compute_likMatrix_HKY()\n";
	// std::cout << "tree is ";
	//print_coaltree();
	//std::cout << "node members are ";
	//print_nodeMembers();
	// std::cout << "isLikelihoodNULL:" <<  isLikelihoodNULL << "\n";
	// if(!isLikelihoodNULL)
	//  	lik.print();
	// std::cout <<"\n";

	if(isLikelihoodNULL)
	{
		Eigen::MatrixXd lik_left(4,lc.get_n_sites_uniq()), lik_right(4,lc.get_n_sites_uniq());
		lik_left.setZero();
		lik_right.setZero();
		// lik_left.initialize(4,lc.get_n_sites_uniq());
		// lik_right.initialize(4,lc.get_n_sites_uniq());
		if(desc[0]->isTip || desc[0]->isLikelihoodNULL)
		{
			lik_left =desc[0]->lik_oneside_HKY(lc,mutrate,kappa);
			//lik_left.replace(desc[0]->lik_oneside_HKY(lc,mutrate,kappa));
		}
		else
		{
			lik_left = desc[0]->lik;
			// lik_left.replace(desc[0]->lik);
		}
		if(desc[0]->isTip || desc[1]->isLikelihoodNULL)
		{
			lik_right = desc[1]->lik_oneside_HKY(lc,mutrate,kappa);
			//lik_right.replace(desc[1]->lik_oneside_HKY(lc,mutrate,kappa));
		}
		else
		{
			lik_right = desc[1]->lik;
			// lik_right.replace(desc[1]->lik);
		}
		lik.resize(lik_left.rows(),lik_left.cols());
		for(int i=0; i<lik.rows();i++)
			for(int j=0; j<lik.cols(); j++)
				lik(i,j) = lik_left(i,j) * lik_right(i,j);
		//lik.replace(matrix_pwProduct(lik_left,lik_right));
		isLikelihoodNULL =  0;

		// REMOVE
		// std::cout << "left child tree is";
		// desc[0]->print_coaltree();
		// std::cout <<"lik_left is\n";
		// lik_left.print();
		// std::cout << "Right child tree is";
		// desc[1]->print_coaltree();
		// std::cout << "lik_right is\n";
		// lik_right.print();
	}

	return;
}


// Computing the likelihood for a subtree (the node and its descendants)
void node_old::compute_likMatrix_HKY(locus lc, double mutrate, double kappa)
{


	if(isLikelihoodNULL)
	{
		Matrix lik_left, lik_right;
		lik_left.initialize(4,lc.get_n_sites_uniq());
		lik_right.initialize(4,lc.get_n_sites_uniq());
		if(desc[0]->isTip || desc[0]->isLikelihoodNULL)
			lik_left.replace(desc[0]->lik_oneside_HKY(lc,mutrate,kappa));
		else
			lik_left.replace(desc[0]->lik);
		if(desc[0]->isTip || desc[1]->isLikelihoodNULL)
			lik_right.replace(desc[1]->lik_oneside_HKY(lc,mutrate,kappa));
		else
			lik_right.replace(desc[1]->lik);
		lik.replace(matrix_pwProduct(lik_left,lik_right));
		isLikelihoodNULL =  0;

		// REMOVE
		// std::cout << "left child tree is";
		// desc[0]->print_coaltree();
		// std::cout <<"lik_left is\n";
		// lik_left.print();
		// std::cout << "Right child tree is";
		// desc[1]->print_coaltree();
		// std::cout << "lik_right is\n";
		// lik_right.print();
	}
	//else
	//{
		// REMOVE
		// std::cout <<"For the tree ";
		// print_coaltree();
		//std::cout <<"The likelihood matrix is already computed:\n";
		// lik.print();
	//}

	return;


}



double locus::getstandfactor (double kappa)
{
  int i, j;
  double p = 0;
  for (i = 0; i < 4; i++)
  {
    for (j = 0; j < 4; j++)
    {
      if (i != j)
      {
        if (i + j == 2 || i + j == 4) // transitions (A<->G, C<->T)
          p += get_pi(i) * get_pi(j) * kappa;
         /*HKY*/
        else // transversions
          p += get_pi(i) * get_pi(j);
      }
    }
  }
  return p;
}



void node::computeTotalLengthsOfLineages()
{
  //std::cout << "In node::computeTotalLengthsOfLineages\n";

  numMutations_overLocus =0.0;
  totalLengthsLineages = 2*age;
  if(isTip == 0)
    {
      desc[0]->computeTotalLengthsOfLineages();
      desc[1]->computeTotalLengthsOfLineages(); 
      if(desc[0]->get_isTip() ==0)
	{      
	  totalLengthsLineages -= desc[0]->get_age();
	  totalLengthsLineages += desc[0]->get_totalLengthsLineages();
	}
      if(desc[1]->get_isTip() ==0)
	{
	  totalLengthsLineages -= desc[1]->get_age();
	  totalLengthsLineages += desc[1]->get_totalLengthsLineages();
	}     
    }
}

void node::computingCumulativeNumMutations_waitingTime(locus lc)
{  
  // std::cout << "computingCumulativeNumMutations_waitingTime\n";

  unsigned int ns=lc.get_n_sites_uniq();
  waitingTime.resize(ns);
  unsigned int stop=0;
  for(unsigned int site =0; site<lc.get_n_sites_uniq() && stop==0; site++)
    {
      initial_forComputingCumulativeNumMutationOverBranches(lc,site);
      computeCumulativeNumMutationsOverBranches(site);

      computingWaitingTime_numMutationsOverLocus(site);
      if(waitingTime.at(site)==-1)
	stop=1;
    }
  return;
}

void node::initial_forComputingCumulativeNumMutationOverBranches(locus lc, unsigned int site)
{
  // std::cout << "In node::initial_forComputingCumulativeNumMutationOverBranches\n";
  if(site==0)
    {      
      unsigned int ns=lc.get_n_sites_uniq();
      cumulativeMutations.resize(ns);
    }
  cumulativeMutations.at(site).resize(2);
  cumulativeMutations.at(site).at(0) = 0; cumulativeMutations.at(site).at(1) = 0;
  if(isTip ==1)
    { 
      cumulativeMutations.at(site).at(lc.get_uniqSeq(site,tipID-1)) = 1;
    }
  else
    {
      desc[0]->initial_forComputingCumulativeNumMutationOverBranches(lc,site);
      desc[1]->initial_forComputingCumulativeNumMutationOverBranches(lc,site);
    }
}


// YC's version
// cumMutations = (x,y): number of tips whose labels 0 or 1, respectively below this node.
void node::computeCumulativeNumMutationsOverBranches(unsigned int site)
{
  if(isTip == 0)
    {
      desc[0]->computeCumulativeNumMutationsOverBranches(site);
      desc[1]->computeCumulativeNumMutationsOverBranches(site);
      
      std::vector<unsigned int> cm0 = desc[0]->get_cumulativeMutations(site);
      std::vector<unsigned int> cm1 = desc[1]->get_cumulativeMutations(site);
      cumulativeMutations.at(site).at(0) = cm0.at(0)+cm1.at(0);
      cumulativeMutations.at(site).at(1) = cm0.at(1)+cm1.at(1);
    }
  return;
}

void node::computingWaitingTime_numMutationsOverLocus(unsigned int site)
{
  if(site>0 && waitingTime.at(site-1)==-1)
    {
      waitingTime.at(site)=-1;
    }
  else
    {
      if(isRoot ==1 )
	{
	  std::vector<unsigned int> cm0 = desc[0]->get_cumulativeMutations(site);
	  std::vector<unsigned int> cm1 = desc[1]->get_cumulativeMutations(site);
	  unsigned int purity = purity_cumulativeMutations(site);
	  // REMOVE
	  /*
	  std::cout << "site=" << site <<"\n";
	  std::cout << "cm0 = " << cm0.at(0) << cm0.at(1) << " and cm1= " << cm1.at(0) << cm1.at(1) <<"\n";
	  std::cout << "purity = " << purity <<"\n";
	  */
	  if(purity ==1)
	    {
	      waitingTime.at(site) = 2*age - desc[0]->get_age() - desc[1]->get_age();
	      /*
	      std::cout << "age = " << age <<" desc[0]->get_age() = "
			<< desc[0]->get_age() << " desc[1]->get_age() = " 
			<< desc[1]->get_age() <<"\n";		
	      std::cout << "waitingTime.at(site) = "<<waitingTime.at(site) <<"\n";
	      */
	      numMutations_overLocus++;
	    }
	  else if(purity ==2)
	    {
	      std::vector<unsigned int> cm;
	      unsigned int char2purify;
	      unsigned int childID2go;
	      if(cm0.at(0)==0|| cm0.at(1)==0)
		{
		  cm = cm0; childID2go = 1;
		}
	      else if(cm1.at(0)==0 || cm1.at(1)==0)
		{
		  cm = cm1; childID2go = 0;		  
		}
	      if(cm.at(0)==0)
		char2purify = 0; // char2purify refers derived alleles
	      else 
		char2purify = 1; 
	      // REMOVE
	      // std::cout << "childID2go = " << childID2go
	      //	<< " char2purify = " << char2purify <<"\n";
	      
	      waitingTime.at(site) = desc[childID2go]->computingWaitingTime_numMutationsOverLocus_recursive(char2purify, site);	
	    }
	  else
	    {
	      waitingTime.at(site) = -1;
	    }
	}
    }
  
  return;
}


void node::initial_forComputingCumulativeNumMutationOverBranches(locus lc)
{
  if(lc.get_nSNPs() != 0)
    {
      unsigned int ns=lc.get_n_sites_uniq();
      cumulativeMutations.resize(ns);
      for(unsigned int site =0; site<lc.get_n_sites_uniq(); site++)
	{
	  cumulativeMutations.at(site).resize(2);
	  cumulativeMutations.at(site).at(0) = 0; cumulativeMutations.at(site).at(1) = 0;
	  if(isTip ==1)
	    { 
	      cumulativeMutations.at(site).at(lc.get_uniqSeq(site,tipID-1)) = 1;
	    }
	  else
	    {
	      desc[0]->initial_forComputingCumulativeNumMutationOverBranches(lc);
	      desc[1]->initial_forComputingCumulativeNumMutationOverBranches(lc);
	    }
	}
    }
  return;
}


// YC's version
// cumMutations = (x,y): number of tips whose labels 0 or 1, respectively below this node.
void node::computeCumulativeNumMutationsOverBranches(locus lc)
{
  unsigned int ns=lc.get_n_sites_uniq();
  for(unsigned int site =0; site<lc.get_n_sites_uniq(); site++)
    {
      if(isTip == 0)
	{
	  desc[0]->computeCumulativeNumMutationsOverBranches(lc);
	  desc[1]->computeCumulativeNumMutationsOverBranches(lc);
	  
	  std::vector<unsigned int> cm0 = desc[0]->get_cumulativeMutations(site);
	  std::vector<unsigned int> cm1 = desc[1]->get_cumulativeMutations(site);
	  cumulativeMutations.at(site).at(0) = cm0.at(0)+cm1.at(0);
	  cumulativeMutations.at(site).at(1) = cm0.at(1)+cm1.at(1);
	}
    }

  return;
}

void node::initial_forComputingNumMutationOverBranches_new(locus lc, unsigned int site)
{
  cumMutations.resize(2);
  cumMutations.at(0) = 0; cumMutations.at(1) = 0;
  if(isTip ==1)
    { 
      cumMutations.at(lc.get_uniqSeq(site,tipID-1)) = 1;
      // numMutations = lc.get_uniqSeq(site,tipID-1);
    }
  else
    {
      // numMutations = -1;
      desc[0]->initial_forComputingNumMutationOverBranches_new(lc,site);
      desc[1]->initial_forComputingNumMutationOverBranches_new(lc,site);
    }
  return;
}

void node::initial_forComputingNumMutationOverBranches(locus lc, unsigned int site)
{
  if(isTip ==1)
    {
      numMutations = lc.get_uniqSeq(site,tipID-1);
    }
  else
    {
      numMutations = -1;
      desc[0]->initial_forComputingNumMutationOverBranches(lc,site);
      desc[1]->initial_forComputingNumMutationOverBranches(lc,site);
    }
  return;
}


// YC's version
// cumMutations = (x,y): number of tips whose labels 0 or 1, respectively below this node.
void node::computeCumNumMutationsOverBranches()
{
  if(isTip == 0)
    {
      desc[0]->computeCumNumMutationsOverBranches();
      desc[1]->computeCumNumMutationsOverBranches();
      
      std::vector<unsigned int> cm0 = desc[0]->get_cumMutations();
      std::vector<unsigned int> cm1 = desc[1]->get_cumMutations();
      cumMutations.at(0) = cm0.at(0)+cm1.at(0);
      cumMutations.at(1) = cm0.at(1)+cm1.at(1);

      // if(desc[0]->get_numMutations() == -1)
      //	desc[0]->computeNumMutationsOverBranches_site_new();
      // if(desc[1]->get_numMutations() == -1)
      //	desc[1]->computeNumMutationsOverBranches_site_new();
    }
  return;
}

// YC's version
// If numMutations = 0 or 1, then it means the number of mutations on this (sub)tree including the parent edge of this node;
//                   2, then it means the mutation happened on the parent edge of this node;
//                   3, then it means that this tree is compatible and the mutation happened on a child branch of this node;
//                   4, then it means that this tree is incompatible.
void node::computeNumMutationsOverBranches_site_new()
{
  if(isTip == 0)
    {
      if(desc[0]->get_numMutations() == -1)
	desc[0]->computeNumMutationsOverBranches_site_new();
      if(desc[1]->get_numMutations() == -1)
	desc[1]->computeNumMutationsOverBranches_site_new();
      double desc0_nMut = desc[0]->get_numMutations();
      double desc1_nMut = desc[1]->get_numMutations();
      std::cout << "In node::computeNumMutationsOverBranches_site_new();";
      std::cout << "the subtree is ";
      print_coaltree();
      std::cout << "desc0_nMut =" << desc0_nMut << " desc1_nMut =" << desc1_nMut <<"\n";
      if(desc0_nMut == 0 && desc1_nMut ==0)
	{
	  numMutations =0;
	}
      else if(desc0_nMut + desc1_nMut == 1) // the mutation happened on one of them
	{
	  numMutations = 3;
	  if(isRoot ==1)
	    {
	      desc[0]->set_numMutations(2);
	      desc[1]->set_numMutations(2);
	      //numMutations_overLocus++; //desc[0]->increase_numMutationsOverLocus_byOne();
	    }
	  /*
	  else
	    {
	      if(desc0_nMut ==1 )
		{
		  desc[0]->set_numMutations(2);
		  //desc[0]->increase_numMutationsOverLocus_byOne();
		}
	      else
		{
		  desc[1]->set_numMutations(2);
		  //desc[1]->increase_numMutationsOverLocus_byOne();
		}
	    }
	  */
	}
      else if(desc0_nMut ==1 && desc1_nMut == 1) // it claims the mutation affects the tips below this node
	{
	  numMutations = 1;
	}
      else if((desc0_nMut<=1 && desc1_nMut <= 3) || (desc0_nMut <=3 && desc1_nMut <=1) ) // (2,0 or 1), (0 or 1,2), (3,0 or 1), (0 or 1,3)
	{
	  numMutations = 3; 
	}
      else
	{
	  numMutations = 4;
	}
    }
  // REMOVE
  std::cout << "sizeoftree = " << size_tree() << " numMutations = " << numMutations << "\n";
  return;
}

void node::computeNumMutationsOverBranches_site()
{
  if(isTip == 0)
    {
      if(desc[0]->get_numMutations() == -1)
	desc[0]->computeNumMutationsOverBranches_site();
      if(desc[1]->get_numMutations() == -1)
	desc[1]->computeNumMutationsOverBranches_site();
      double desc0_nMut = desc[0]->get_numMutations();
      double desc1_nMut = desc[1]->get_numMutations();
      // std::cout << "desc0_nMut =" << desc0_nMut << " desc1_nMut =" << desc1_nMut <<"\n";
      if(desc0_nMut != numMutations)
	{
	  if(desc0_nMut == 2)
	    {
	      if(desc1_nMut== -1) 
		numMutations = 2;
	      else
		numMutations = desc1_nMut;
	    }
	  else if(desc1_nMut != -1 && desc1_nMut != 2 && (desc0_nMut != desc1_nMut))
	    {  // if (desc0_nMut, desc1_nMut) = (1, 0) or (0,1)
	      numMutations = 2;
	    }
	  else
	    {
	      numMutations = desc0_nMut;
	    }
	}
      else if(numMutations == 2)
	{
	  if(desc1_nMut != -1 && desc1_nMut != 2) 
	    numMutations = desc1_nMut;
	}
    }
  // std::cout << "sizeoftree = " << size_tree() << " numMutations = " << numMutations << "\n";
  return;
}

/*
double node::computeNumMutationsOverBranches(locus lc)
{  
  for(unsigned int site =0; site < lc.get_n_sites_uniq(); site++)
    {
      initial_forComputingNumMutationOverBranches(lc, site);
      computeNumMutationsOverBranches_site();
    }
}
*/
 
unsigned int node::purity_cumulativeMutations(unsigned int s)
{
  unsigned int purity = 0;
  std::vector<unsigned int> cm0 = desc[0]->get_cumulativeMutations(s);
  std::vector<unsigned int> cm1 = desc[1]->get_cumulativeMutations(s);
  if((cm0.at(0) == 0 && cm1.at(1) ==0) || (cm0.at(1) == 0 && cm1.at(0) ==0))
    {
      purity = 1;
    }
  else if((cm0.at(0) == 0 || cm1.at(1) ==0) || (cm0.at(1) == 0 || cm1.at(0) ==0))
    {
      purity = 2; // possibly 
    }
  return purity;
}

unsigned int node::purity_cumMutations()
{
  unsigned int purity = 0;
  std::vector<unsigned int> cm0 = desc[0]->get_cumMutations();
  std::vector<unsigned int> cm1 = desc[1]->get_cumMutations();
  if((cm0.at(0) == 0 && cm1.at(1) ==0) || (cm0.at(1) == 0 && cm1.at(0) ==0))
    {
      purity = 1;
    }
  else if((cm0.at(0) == 0 || cm1.at(1) ==0) || (cm0.at(1) == 0 || cm1.at(0) ==0))
    {
      purity = 2; // possibly 
    }
  return purity;
}

double node::computeWaitingTime_numMutationsOverLocus_recursive(unsigned int char2purify)
{
  double waitingTime = 0.0;
  if(isTip ==0)
    {
      unsigned int purity = purity_cumMutations();
      std::vector<unsigned int> cm0 = desc[0]->get_cumMutations();
      std::vector<unsigned int> cm1 = desc[1]->get_cumMutations();
      // REMOVE
      //std::cout << "cm0 = " << cm0.at(0) << cm0.at(1) << " and cm1= " << cm1.at(0) << cm1.at(1) <<"\n";
      //std::cout << "purity = " << purity <<"\n";
      if(purity ==1)
	{
	  if(cm0.at(char2purify) == 0)
	    {	      
	      waitingTime = age - desc[0]->get_age();
	      desc[0]->increase_numMutationsOverLocus_byOne();	      
	    }
	  else
	    {	      
	      waitingTime = age - desc[1]->get_age();
	      desc[1]->increase_numMutationsOverLocus_byOne();	      
	    }
	}
      else if( purity ==2)
	{
	  std::vector<unsigned int> cm;
	  unsigned int char2purify_new;
	  unsigned int childID2go;
	  if(cm0.at(0)==0|| cm0.at(1)==0)
	    {
	      cm = cm0;
	      childID2go = 1;
	    }
	  else if(cm1.at(0)==0 || cm1.at(1)==0)
	    {
	      cm = cm1;
	      childID2go = 0;		  
	    }
	  if(cm.at(0)==1)
	    char2purify_new = 0;
	  else 
	    char2purify_new = 1;


	  if(char2purify_new == char2purify)
	    waitingTime = desc[childID2go]->computeWaitingTime_numMutationsOverLocus_recursive(char2purify);	
	  else
	    {
	      waitingTime = -1;
	    }
	}
      else
	{
	  waitingTime = -1;
	}
    }
  return waitingTime;
}

double node::computingWaitingTime_numMutationsOverLocus_recursive(unsigned int char2purify, unsigned int s)
{
  // std::cout << "node::computingWaitingTime_numMutationsOverLocus_recursive\n";

  double wTime = 0.0;
  if(isTip ==0)
    {
      unsigned int purity = purity_cumulativeMutations(s);
      std::vector<unsigned int> cm0 = desc[0]->get_cumulativeMutations(s);
      std::vector<unsigned int> cm1 = desc[1]->get_cumulativeMutations(s);
      // REMOVE
      // std::cout << "cm0 = " << cm0.at(0) << cm0.at(1) << " and cm1= " << cm1.at(0) << cm1.at(1) <<"\n";
      // std::cout << "purity = " << purity <<"\n";
      if(purity ==1)
	{
	  if(cm0.at(char2purify) > 0) // this clades contains derived alleles.
	    {	      
	      wTime = age - desc[0]->get_age();
	      desc[0]->increase_numMutationsOverLocus_byOne();	      
	    }
	  else
	    {	      
	      wTime = age - desc[1]->get_age();
	      desc[1]->increase_numMutationsOverLocus_byOne();	      
	    }
	}
      else if( purity ==2)
	{
	  std::vector<unsigned int> cm;
	  unsigned int char2purify_new;
	  unsigned int childID2go;
	  if(cm0.at(0)==0|| cm0.at(1)==0)
	    {
	      cm = cm0;
	      childID2go = 1;
	    }
	  else if(cm1.at(0)==0 || cm1.at(1)==0)
	    {
	      cm = cm1;
	      childID2go = 0;		  
	    }
	  if(cm.at(0)==0)
	    char2purify_new = 0;
	  else 
	    char2purify_new = 1;


	  if(char2purify_new == char2purify)
	    wTime = desc[childID2go]->computingWaitingTime_numMutationsOverLocus_recursive(char2purify, s);	
	  else
	    {
	      wTime = -1;
	    }
	}
      else
	{
	  wTime = -1;
	}
    }
  return wTime;
}

/*
double node::computingWaitingTime_numMutationsOverLocus(locus lc)
{
  unsigned int ns=lc.get_n_sites_uniq();
  waitingTime.resize(ns);
  for(unsigned int site =0; site<lc.get_n_sites_uniq(); site++)
    {
      if(site>0 && waitingTime.at(site-1)==-1){
	waitingTime.at(site)=-1;
      }else{
	if(isRoot ==1 )
	  {
	    std::vector<unsigned int> cm0 = desc[0]->get_cumulativeMutations(site);
	    std::vector<unsigned int> cm1 = desc[1]->get_cumulativeMutations(site);
	    unsigned int purity = purity_cumulativeMutations(site);
	    // REMOVE
	    //std::cout << "cm0 = " << cm0.at(0) << cm0.at(1) << " and cm1= " << cm1.at(0) << cm1.at(1) <<"\n";
	    //std::cout << "purity = " << purity <<"\n";
	    if(purity ==1)
	      {
		waitingTime.at(site) = 2*age - desc[0]->get_age() - desc[1]->get_age();
		numMutations_overLocus++;
	      }
	    else if( purity ==2)
	      {
		std::vector<unsigned int> cm;
		unsigned int char2purify;
		unsigned int childID2go;
		if(cm0.at(0)==0|| cm0.at(1)==0)
		  {
		    cm = cm0;
		    childID2go = 1;
		  }
		else if(cm1.at(0)==0 || cm1.at(1)==0)
		  {
		    cm = cm1;
		    childID2go = 0;		  
		  }
		if(cm.at(0)==1)
		  char2purify = 0;
		else 
		  char2purify = 1;
		// REMOVE
		//std::cout << "childID2go = " << childID2go
		//	    << " char2purify = " << char2purify <<"\n";
		
		waitingTime.at(site) = desc[childID2go]->computingWaitingTime_numMutationsOverLocus_recursive(char2purify, site);	
	      }
	    else
	      {
		waitingTime.at(site) = -1;
	      }
	  }
      }
    }
  return;
}
*/

double node::computeWaitingTime_numMutationsOverLocus_new()
{
  double waitingTime = 0.0;
  if(isRoot ==1 )
    {
      std::vector<unsigned int> cm0 = desc[0]->get_cumMutations();
      std::vector<unsigned int> cm1 = desc[1]->get_cumMutations();
      unsigned int purity = purity_cumMutations();
      // REMOVE
      //std::cout << "cm0 = " << cm0.at(0) << cm0.at(1) << " and cm1= " << cm1.at(0) << cm1.at(1) <<"\n";
      //std::cout << "purity = " << purity <<"\n";
      if(purity ==1)
	{
	  waitingTime = 2*age - desc[0]->get_age() - desc[1]->get_age();
	  numMutations_overLocus++;
	}
      else if( purity ==2)
	{
	  std::vector<unsigned int> cm;
	  unsigned int char2purify;
	  unsigned int childID2go;
	  if(cm0.at(0)==0|| cm0.at(1)==0)
	    {
	      cm = cm0;
	      childID2go = 1;
	    }
	  else if(cm1.at(0)==0 || cm1.at(1)==0)
	    {
	      cm = cm1;
	      childID2go = 0;		  
	    }
	  if(cm.at(0)==1)
	    char2purify = 0;
	  else 
	    char2purify = 1;
	  // REMOVE
	  //std::cout << "childID2go = " << childID2go
	  //	    << " char2purify = " << char2purify <<"\n";

	  waitingTime = desc[childID2go]->computeWaitingTime_numMutationsOverLocus_recursive(char2purify);	
	}
      else
	{
	  waitingTime = -1;
	}
    }
  return waitingTime;
}

double node::computeWaitingTime_new()
{
  double waitingTime = 0.0;
  if(numMutations == 4)
    {
      waitingTime = -1;
    }
  else
    {
      if(numMutations == 3)
	{
	  double desc0_numMut = desc[0]->get_numMutations();
	  double desc1_numMut = desc[1]->get_numMutations();
	  if(isRoot ==1 && (desc0_numMut==2|| desc1_numMut==2) ) // root node
	    {
	      waitingTime = 2*age - desc[0]->get_age() + desc[1]->get_age();
	      numMutations_overLocus++;
	    }
	  else
	    {
	      if(desc0_numMut == 3 || desc0_numMut ==2)
		{
		  waitingTime = desc[0]-> computeWaitingTime_new();
		}
	      else if(desc1_numMut == 3 || desc1_numMut ==2)
		{
		  waitingTime = desc[1]-> computeWaitingTime_new();		  
		}
	      else if(desc0_numMut + desc1_numMut == 1) // (1,0), (0,1)
		{
		  unsigned int sisNode = 0;
		  if(par->desc[sisNode]->get_age() == age)
		    sisNode = 1;
		  double sisNode_numMut = par->desc[sisNode]->get_numMutations();
		  if(desc0_numMut != sisNode_numMut)
		    {
		      waitingTime = age - desc[0]->get_age();
		      desc[0]->increase_numMutationsOverLocus_byOne();
		    }
		  else if(desc1_numMut != sisNode_numMut)
		    {
		      waitingTime = age - desc[1]->get_age();
		      desc[1]->increase_numMutationsOverLocus_byOne();
		    }
		  else
		    {
		      std::cout << "Error in node::computeWaitingTime_numMutationsOverLocus_new()";
		      std::cout << "desc0_numMut = " << desc0_numMut  << " and desc1_numMut = "<<desc1_numMut <<" and sisNode_numMut = " << sisNode_numMut <<"\n";
		    }
		}
	      
	      if((desc0_numMut == 3 || desc0_numMut ==2) && (desc1_numMut == 3 || desc1_numMut ==2))
		{
		  std::cout << "Error in node::computeWaitingTime_numMutationsOverLocus_new()";
		  std::cout << "desc0_numMut = " << desc0_numMut  << " and desc1_numMut = "<<desc1_numMut <<"\n";
		}
	    }
	}
      if(numMutations ==2)
	{
	  if(par->get_isRoot() ==1 )
	    {
	      std::cout << "Error in node::computeWaitingTime_numMutationsOverLocus_new()";
	    }
	  else
	    waitingTime = par->get_age()-age;	  
	}
    }

  return waitingTime;
}

double node::computeWaitingTime_numMutationsOverLocus()
{
  double waitingTime = 0.0;
  if(isTip != 1)
    {
      double desc0_numMut = desc[0]->get_numMutations();
      double desc1_numMut = desc[1]->get_numMutations();
      if((desc0_numMut == 0 && desc1_numMut == 1) || (desc0_numMut == 1 && desc1_numMut == 0))
	{
	  if(isRoot == 1) // root node
	    {	    
	      if(desc[0]->get_isTip() == 1)
		waitingTime = age;
	      else
		waitingTime = age - desc[0]->get_age();
	      if(desc[1]->get_isTip() ==1)
		waitingTime += age;
	      else
		waitingTime += age-desc[1]->get_age();
	    }
	  else // not root node
	    {
	      if(par->get_numMutations() == desc[0]->get_numMutations())
		{
		  if(desc[1]->get_isTip() ==1)
		    waitingTime = age;
		  else
		    waitingTime = age - desc[1]->get_age();
		}
	      else
		{
		  if(desc[0]->get_isTip() ==1)
		    waitingTime = age;
		  else
		    waitingTime = age - desc[0]->get_age();
		}
	    }
	  // FIXME - YC 
	  // Here the way to compute the numMutations_overLocus is a bit different from IMa2
	  if(isRoot==1)
	    numMutations_overLocus++;
	  else
	    {
	      if(desc0_numMut == 1)
		desc[0]->increase_numMutationsOverLocus_byOne();
	      else
		desc[1]->increase_numMutationsOverLocus_byOne();
	    }
	    
	  //	  std::cout << " numMutations_overLocus = " <<  numMutations_overLocus << " In double node::computeWaitingTime_numMutationsOverLocus()\n";
	}
      else
	{
	  double waitingTime1 = desc[0]->computeWaitingTime_numMutationsOverLocus();
	  double waitingTime2 = desc[1]->computeWaitingTime_numMutationsOverLocus();
	  // REMOVE
	  // std::cout << "waitingTime1 = "<< waitingTime1 << " and waitingTime2 = " << waitingTime2 <<"\n";
	  if(waitingTime1 >0 && waitingTime2 >0 || waitingTime1 <0 || waitingTime2 < 0) // incompatible tree
	    {
	      waitingTime = -1;
	    }
	  else
	    {
	      waitingTime = waitingTime1+waitingTime2;
	    }
	}      
    }

  return waitingTime;
}


double node::compute_partitionFunctionOfPossionDist_numMutations()
{
  // REMOVE
  // std::cout << "numMutations_overLocus = " << numMutations_overLocus << "\n";
  double logNumMutations_overLocus = logfactorial(numMutations_overLocus);
  // double logNumMutations_overLocus = log(factorial(numMutations_overLocus));
  if(isTip != 1)
    {
      logNumMutations_overLocus += desc[0]->compute_partitionFunctionOfPossionDist_numMutations();      
      logNumMutations_overLocus += desc[1]->compute_partitionFunctionOfPossionDist_numMutations();
    }
  return logNumMutations_overLocus;
}

void node::print_numMutations() // for debugging
{
  
  if(isRoot ==1)
    {
      std::cout << "(";
      desc[0]->print_numMutations();
      std::cout << ",";
      desc[1]->print_numMutations();
      std::cout << "):" << numMutations;
      std::cout <<"\n";
    }
  else if(isTip ==1 )
    {
      std::cout << numMutations;
    }
  else
    {
      std::cout << "(";
      desc[0]->print_numMutations();
      std::cout << ",";
      desc[1]->print_numMutations();
      std::cout << "):" << numMutations;
    }
}

//double node::logLikelihood_IS (int ci, int li, double mutrate)

/* written by Rasmus to return the following value
 - (total length of gtree in time) * theta + SUM(log(t_i * theta)*k_i, over all branches)
  where k_i is the number of mutations on branch i. 
  This is equal to the logarithm of the product of a bunch of poisson variables, one for each branch, except that 
  one term is dropped out. [ -SUM(log(k_i !)  ]  is not included.  However this will be the same for two infinite 
  sites gtrees that fit the same data set.
 */
 /* under new parameterization
    - (total length of gtree in time) * mutrate + SUM(log(t_i * mutrate)*k_i, over all branches) - SUM(log(k_i !) 
    where k_i is the number of mutations on branch i. 
    as before, do not need to include  the term : -SUM(log(k_i !) 

    Can also do it as 
    - (total length of gtree in time) * mutrate + log(mutrate) * SUM(k_i) + Sum(log(t_i)*k_i, over all branches) - SUM(log(k_i !) 
    check out the method of getting gtree length - seems ok.
  */
double node::logLikelihood_IS (locus lc, double mutrate, unsigned int updateMut)
{  
  //if(updateMut==1)
  //  std::cout << "\nIn node::logLikelihood_IS()\n";
  // print_coaltree();

  // computeTotalLengthsOfLineages();
  //if(updateMut==1)
  //  std::cout << "computingTotalLengthsOfLineages..Done\n";

  double loglik = -1 * totalLengthsLineages * mutrate;
  
  // double waitingTime = 0.0;
  double sumlogNumMutationsFactorial = 0.0;
  if(lc.get_nSNPs() != 0)
    {
      if(updateMut==0)
	{
	  //initial_forComputingCumulativeNumMutationOverBranches(lc);
	  //computeCumulativeNumMutationsOverBranches(lc);
	  //computingWaitingTime_numMutationsOverLocus(lc);
	  
	  computingCumulativeNumMutations_waitingTime(lc);
	}

      /*
      if(updateMut==1)
	{
	  std::cout << "lc.get_n_sites_uniq() = " << lc.get_n_sites_uniq()<<"\n";
	  std::cout << "waitingTime.size() = " << waitingTime.size() <<"\n";
	}
      */

      unsigned int stop=0;
      for(unsigned int site =0; site<lc.get_n_sites_uniq() && stop==0; site++)
	{
	  // initial_forComputingNumMutationOverBranches_new(lc, site);
	  // computeCumNumMutationsOverBranches();
	  // REMOVE
	  // print_numMutations();

	  // waitingTime = computeWaitingTime_numMutationsOverLocus_new();

	  //if(updateMut==1)
	  //  std::cout <<"site="<< site << " waitingTime = " << waitingTime.at(site) << "\n";

	  if(waitingTime.at(site) == -1)
	    {
	      stop=1;
	      // std::cout << "waitingTime = " << waitingTime << "\n";
	    }
	  else if(waitingTime.at(site) > 0)
	    {
	      loglik += log(waitingTime.at(site) * mutrate);
	    }
	  // REMOVE   
	  /*
	  std::cout << "site = " << site <<": ";	 
	  for(unsigned int i=0; i<size_tree();i++)
	    std::cout << lc.get_uniqSeq(site,i) <<" ";	  
	  std::cout << " waitingTime = " << waitingTime.at(site) << "\n";
	  std::cout <<"loglik = " << loglik <<"\n";
	  */
	}
      if(stop == 1)
	{
	  loglik = -1*numeric_limits<double>::infinity();
	  rejectThisTree_IS = 1;
	}
      else
	{
	  sumlogNumMutationsFactorial = compute_partitionFunctionOfPossionDist_numMutations();
	  loglik -= sumlogNumMutationsFactorial;
	  
	  // REMOVE
	  /*
	  if(loglik == -1*numeric_limits<double>::infinity())
	    std::cout << "sumlogNumMutationsFactorial = " << sumlogNumMutationsFactorial <<"\n";
	  */
	}
    }

  // REMOVE
  /*
  if(loglik == -1*numeric_limits<double>::infinity())
    {
      std::cout << "lc.get_n_sites_uniq() = " << lc.get_n_sites_uniq()<<"\n";
      std::cout << "totalLengthsLineages= " <<totalLengthsLineages
		<< " and mutrate = " << mutrate <<"\n";
      std::cout << "loglik = " << loglik << " of tree ";
      print_coaltree();
      if(rejectThisTree_IS ==0)
	{
	  std::cout << "totalLengthsLineages= " << totalLengthsLineages 
		    << " mutrate = " << mutrate 
		    << " sumlogNumMutationsFactorial  = " << sumlogNumMutationsFactorial
		    << "\n";
	}
      
      for(unsigned int site =0; site<lc.get_n_sites_uniq(); site++)
	{
	  std::cout << "site = " << site <<": ";	 
	  for(unsigned int i=0; i<size_tree();i++)
	    std::cout << lc.get_uniqSeq(site,i) <<" ";	  
	  std::cout << " waitingTime = " << waitingTime.at(site) << "\n";
	}
    }
  */

  return loglik;
}      // likelihood IS


/**
 * Computing the likelihood value under HKY model
 * This function is called when the likelihood value of an initial genealogy is computed
 */
double node::logLikelihood_HKY (locus lc, double mutrate, double kappa)
{
  unsigned int nSites = lc.get_nSites();
	Eigen::MatrixXd pi_mat(1,4); // Empirical distribution
	pi_mat.setZero();
	//pi_mat.initialize(1,4);
	for(int i=0; i<4; i++)
	{
		pi_mat(0,i) = lc.get_pi()[i];
	}
	//pi_mat.replace(lc.get_pi(), 1, 4); // (1 by 4) matrix


	// FIXME YC 3/14/2014
	// "mutrate" is the number of mutations per locus.
	// "new_mutRate" is rescaled, so it is the number of mutations per site.
	// I don't understand why it is divided by getstandfactor().
	// See the original code "likelihoodHKY()" in "calc_prob_data.cpp" of Ima2.
	double new_mutRate = mutrate / ((double) nSites * lc.getstandfactor (kappa));


	compute_likMatrix_HKY(lc,new_mutRate,kappa);
	// (unique) site-likelihood conditional on root's nucleotide
	// (4 by n_sites_uniq) matrix
	Eigen::MatrixXd likMat(4,lc.get_n_sites_uniq());//  = get_likMatrix();
	//likMat.initialize(4,lc.get_n_sites_uniq());
	likMat = get_likMatrix();

	// (unique) site-likelihood
	// (1 by n_sites_uniq) matrix
	Eigen::MatrixXd lik(1,lc.get_n_sites_uniq());// =matrixMultiplication(pi_mat,likMat);
	//lik.initialize(1,lc.get_n_sites_uniq());
	lik =pi_mat*likMat;
	//lik =matrixMultiplication(pi_mat,likMat);


	double loglik = 0.0;
	for(unsigned int i=0; i<lc.get_n_sites_uniq(); i++)
	{
		// std::cout <<"at site "<<i<<", lik = "<<lik.val(0,i)<<"\n";
		try{
			double lik_val =lik(0,i);
			loglik += log(lik_val) * (double) (lc.get_freq_uniqueSeq()).at(i);
			//loglik += log(lik.val(0,i)) * (double) (lc.get_freq_uniqueSeq()).at(i);
		} catch (std::exception &e) {
			std::cout << "Can't access matrix lik or the frequency of sequences to calculate likelihood - vector index out of bounds\n";
		}

		// REMOVE later
		// std::cout << "i: " << i
		//		<< " log(lik.val(0,i))" << log(lik.val(0,i))
		//		<< " (lc.get_freq_uniqueSeq()).at(i)" << (lc.get_freq_uniqueSeq()).at(i)
		//		<< " partial loglik = " << loglik << "\n";

	}

	// REMOVE later
	// std::cout << "\t Final loglik = " << loglik << "\n";

	return loglik;

}

/**
 * Computing the likelihood value under HKY model
 */
double node_old::logLikelihood_HKY (locus lc, double mutrate, double kappa)
{
  unsigned int nSites = lc.get_nSites();
	Matrix pi_mat; // Empirical distribution
	pi_mat.initialize(1,4);
	pi_mat.replace(lc.get_pi(), 1, 4); // (1 by 4) matrix

	// FIXME YC 3/14/2014
	// "mutrate" is the number of mutations per locus.
	// "new_mutRate" is rescaled, so it is the number of mutations per site.
	// I don't understand why it is divided by getstandfactor().
	// See the original code "likelihoodHKY()" in "calc_prob_data.cpp" of Ima2.
	double new_mutRate = mutrate / ((double) nSites * lc.getstandfactor (kappa));

	compute_likMatrix_HKY(lc,new_mutRate,kappa);
	// (unique) site-likelihood conditional on root's nucleotide
	// (4 by n_sites_uniq) matrix
	Matrix likMat;//  = get_likMatrix();
	likMat.initialize(4,lc.get_n_sites_uniq());
	likMat = get_likMatrix();

	// (unique) site-likelihood
	// (1 by n_sites_uniq) matrix
	Matrix lik;// =matrixMultiplication(pi_mat,likMat);
	lik.initialize(1,lc.get_n_sites_uniq());
	lik =matrixMultiplication(pi_mat,likMat);

	// REMOVE
	//std::cout << "in logLikelihood_HKY()\n";
	//pi_mat.print();
	//likMat.print();
	//lik.print();

	double loglik = 0.0;
	for(unsigned int i=0; i<lc.get_n_sites_uniq(); i++)
	{
		// std::cout <<"at site "<<i<<", lik = "<<lik.val(0,i)<<"\n";
		try{
			loglik += log(lik.val(0,i)) * (double) (lc.get_freq_uniqueSeq()).at(i);
		} catch (std::exception &e) {
			std::cout << "Can't access matrix lik or the frequency of sequences to calculate likelihood - vector index out of bounds\n";
		}

		// REMOVE later
		// std::cout << "i: " << i
		//		<< " log(lik.val(0,i))" << log(lik.val(0,i))
		//		<< " (lc.get_freq_uniqueSeq()).at(i)" << (lc.get_freq_uniqueSeq()).at(i)
		//		<< " partial loglik = " << loglik << "\n";

	}

	// REMOVE later
	// std::cout << "\t Final loglik = " << loglik << "\n";

	return loglik;

}
