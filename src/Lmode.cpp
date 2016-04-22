/* MIST copyright 2016 by Yujin Chung and Jody Hey */

/*
 * Lmode.cpp
 *
 */

#include "coaltree.hpp"
#include "popTree.hpp"
#include "Chain.hpp"
#include "misc.hpp"
#include <math.h>
#include <iostream>
#include <complex>
#include <limits>
#include <chrono>
#include "Eigen/Dense"
#include <complex>

void node::assignPopulations2Tips(locus lc)
{
	if(isTip)
	{
		// REMOVE
		// std::cout << "tipID = " << tipID <<"\n";
		// std::cout << "popNames = " << lc.popNames.at(tipID-1) <<"\n";

		stringstream convert(lc.popNames.at(tipID-1));
		convert >> popID;
		label = popID-1;
	}
	else
	{
		desc[0]->assignPopulations2Tips(lc);
		desc[1]->assignPopulations2Tips(lc);
	}
}

 void node::get_tipState(Eigen::VectorXi &state, int nPops)
{
	 // REMOVE
	//	std::cout << "in node::get_tipState()."
	//			" node tree is ";
	//	print_coaltree();
	//	std::cout << "isTip == " << isTip <<"\n";
	//	std::cout << "label is " << label << "; popID is " << popID <<"\n";
	 //	std::cout << "size of tree is " << size_tree() <<"\n";
	 //	std::cout << "age is " << age <<"\n";

	if(isTip==1 && popID < nPops)
	{
		// REMOVE
		// std::cout << "in node::get_tipState()."
		// 		" node tree is ";
		// print_coaltree();
		// std::cout << "label is " << label << "; popID is " << popID <<"\n";

		state(label*(popID-1)+label)++;
	}
	else if(age ==0 && size_tree()==2 && popID < nPops)
	{
		state(desc[0]->label*(popID-1)+desc[0]->label) ++;
		state(desc[1]->label*(popID-1)+desc[1]->label) ++;
	}
	else if(isTip == 0 && age > 0.0)
	{
		desc[0]->get_tipState(state,nPops);
		desc[1]->get_tipState(state,nPops);
	}

	return;
}



Eigen::VectorXi node::get_totalNumEachKind(int noPops)
{
	// REMOVE
	// print_coaltree();

	Eigen::VectorXi state(noPops);
	state.setZero();

	if(isTip)
	{
		// REMOVE
		// std::cout << "popID is "<< popID <<"\n";

		state(label)++;

		//REMOVE
		// for(unsigned int i=0; i< (unsigned) noPops; i++)
		// 	std::cout << "state at i = "<<  state.at(i)<<"\n";
	}
	else
	{
		Eigen::VectorXi state1 = desc[0]->get_totalNumEachKind(noPops);
		Eigen::VectorXi state2 = desc[1]->get_totalNumEachKind(noPops);
		state += state1+state2;

		//for(unsigned int i=0; i< (unsigned) noPops; i++)
		//{
			// REMOVE
			// std::cout << state.at(i) << " "<< state1.at(i) << " " << state2.at(i)<<"\n";
		//   state.at(i) += state1.at(i)+state2.at(i);
		//}
	}
	return state;
}


/***
 *  The size of the state space.
 *  Note that a state is the number of each kind in each population
 *  @param noLineages the total number of each kind
 *  @param noPops the number of populations at a given time
 */
int compute_spaceSize(Eigen::VectorXi noLineages, int noPops)
{
	int size = 0;
	if(noPops == 2)
	{
		size = (noLineages(0)+1) * (noLineages(1)+1);
	}
	else if(noPops == 3) // FIXME YC 3/21/2014
	{					 // what if the number of populations >= 4
		size = 0;
		for(int n1 =0; n1<=noLineages(0); n1++)
			for(int n2 =0; n2<=noLineages(1); n2++)
				for(int n3 =0; n3<=noLineages(2); n3++)
				{
					Eigen::VectorXi new_noLineages =  noLineages;
					new_noLineages(0) -= n1;
					new_noLineages(1) -= n2;
					new_noLineages(2) -= n3;
					size += compute_spaceSize(new_noLineages, noPops-1);

					// REMOVE
					// std::cout <<compute_spaceSize(new_noLineages, noPops-1) <<"\n";
				}
	}
	return size;
}

/***
 * A state denote the number of lineages in each population.
 * This function generates the state space.
 * @param noLineages the total number of each kind
 * @param noPops the number of populations at a given time
 */
Eigen::MatrixXi getSpace(Eigen::VectorXi noLineages, int noPops)
{

	int spaceSize = compute_spaceSize(noLineages,noPops);
	Eigen::MatrixXi space(spaceSize,noLineages.size()*(noPops-1));

	// REMOVE
	//std::cout << "In Eigen::MatrixXi getSpace(). spaceSize = " << spaceSize <<"\n";
	//std::cout << "noLineages is..\n" << noLineages
	//		<< "\n noPops = " << noPops <<"\n";

	if(noPops == 2)
	{
		int idx_row = 0;
		for(int n1=0; n1<= noLineages(0); n1++)
			for(int n2=0; n2<= noLineages(1); n2++)
			{
				if(noLineages.size()==2)
				{
					space.row(idx_row) = Eigen::Vector2i(noLineages(0)-n1,noLineages(1)-n2);
					idx_row++;
				}
				else
				{
					for(int n3=0; n3<= noLineages(2); n3++)
					{
						space.row(idx_row) = Eigen::Vector3i(noLineages(0)-n1,noLineages(1)-n2,noLineages(2)-n3);
						idx_row++;
					}
				}

				// REMOVE
				// std::cout << "n1 = " << n1 << "; n2= " << n2 << "; (n1+1)*(n2+1)-1 = "<<(n1+1)*(n2+1)-1  <<"\n\t "
				//		<< Eigen::Vector2i(noLineages.at(0)-n1,noLineages.at(1)-n2)<<"\n";
			}
	}
	else if(noPops >= 3)
	{
		Eigen::MatrixXi subSpace1 = getSpace(noLineages, noPops-1);
		// REMOVE
		//std::cout << "The first part of space is \n" << subSpace1 <<"\n";

		int idx_row = 0;
		for(int i=0; i<subSpace1.rows(); i++)
		{
			Eigen::VectorXi row = subSpace1.row(i);
			Eigen::MatrixXi subSpace2 = getSpace(noLineages-row, noPops-1);
			// REMOVE
			//std::cout << "row is \n" << row
			//		<< "\n space for pop2 is \n" << subSpace2 <<"\n";

			int subSpaceSize2 = subSpace2.rows();
			space.block(idx_row,noPops,subSpaceSize2,noLineages.size()*(noPops-2)) = subSpace2;
			for(int j=0; j<subSpaceSize2; j++)
				space.block(idx_row+j,0,1,noLineages.size()) = row.transpose();


			idx_row += subSpaceSize2;
		}
	}
	return space;
}

/***
 * Return (-1,-1, 0)                                   if no migration can happen
 *        (populationID a migrant is from,
 *        	   to population ID a migrant is to,
 *        	   		the number of candidate migrants)   if a migration can happen
 * Note that it is considered backward in time.
 * Note that the population Ids here start from 0.
 */
Eigen::Vector3i canMigrate(Eigen::VectorXi state1, Eigen::VectorXi state2, int nPops, Eigen::VectorXi totalN)
{
	Eigen::VectorXi diff = state1-state2;

	// REMOVE
	//std::cout << "diff= " << diff <<"\n";

	int nomig = 0;
	int id_fromPop = -1;
	int id_toPop = -1;
	int noLin_mig = 0;
	int idx_migrants = -1; // This index will use to compute the number of candidate migrants
                           // when id_fromPop = nPop-1

	int idx_pop = 0;
	while(nomig==0 && idx_pop < nPops-1)
	{
		Eigen::VectorXi seg = diff.segment(idx_pop*nPops,nPops);
		if(seg.squaredNorm() > 1)
			nomig = 1;
		else
		{
			int i=0;
			while(i<nPops && nomig==0)
			{
				// REMOVE
				//std::cout << "idx_pop = " << idx_pop
				//		<< " i = "<< i << "\n";

			if(diff(idx_pop*i+i)==1 && id_fromPop==-1) // difference value of 1 should be found only once.
			{
				if(idx_migrants == -1 || idx_migrants == i)
				{
					id_fromPop = idx_pop;
					noLin_mig =	state1(idx_pop*i+i);
					idx_migrants = i;
				}
				else
					nomig = 1;
			}
			else if(diff(idx_pop*i+i)==-1 && id_toPop==-1) // difference value of -1 should be found only once.
			{
				if(idx_migrants == -1|| idx_migrants == i)
				{
					idx_migrants = i;
					id_toPop = idx_pop;
				}
				else
					nomig = 1;
			}
			else if(diff(idx_pop*i+i) != 0)
					nomig = 1;
			i++;
			}
		}
		idx_pop++;
	}

	// REMOVE
	//std::cout << "nomig = " << nomig
	//		<< " id_fromPop = " << id_fromPop
	//		<< " id_toPop = " << id_toPop << "\n";

	Eigen::Vector3i popIDs;
	if(nomig==1)
	{
		popIDs << -1, -1, 0;
	}
	else
	{
		if(id_fromPop == -1 && id_toPop >= 0)
		{
			noLin_mig = totalN(idx_migrants);
			for(int i=0; i< nPops-1; i++)
				noLin_mig -= state1(i*nPops+idx_migrants);
			popIDs << nPops-1, id_toPop, noLin_mig;
		}
		else if(id_fromPop >=0 && id_toPop == -1)
			popIDs << id_fromPop, nPops-1, noLin_mig;
		else if(id_fromPop >=0 && id_toPop >=0)
			popIDs << id_fromPop, id_toPop, noLin_mig;
		else
		{
			popIDs << -1, -1, 0;
			std::cout << "Error in canMigrate().\n";
		}

	}
	return popIDs;
}



double popTree::find_migRate(unsigned int popID_fr, unsigned int popID_to)
{
	double rate = -1.0;

	unsigned int crrID = popID;
	if(crrID == popID_fr)
	{
		unsigned int i=0;
		unsigned int found = 0;
		while(found == 0 && i<mig.size())
		{
			if(mig.at(i).get_popID() == popID_to)
			{
				found = 1;
				rate = mig.at(i).get_migRate();
			}
			else
				i++;
		}
	}
	else if(!isTip)
	{
		rate = desc[0]->find_migRate(popID_fr,popID_to);
		if(rate <0)
			rate = desc[1]->find_migRate(popID_fr,popID_to);
	}

	return rate;
}


void popTree::find_popSizeTips(Eigen::VectorXd &popSizes)
{
	if(isTip)
		popSizes(popID-1) = populationSize;
	else
	{
		desc[0]->find_popSizeTips(popSizes);
		desc[1]->find_popSizeTips(popSizes);
	}
}


Eigen::MatrixXd popTree::getTransitionRateMat(Eigen::MatrixXi space, Eigen::VectorXi totalN)
{
	int nrow = space.rows();
	Eigen::MatrixXd rateMatrix(nrow+1,nrow+1);
	rateMatrix.setZero();

	int npops = size();

	// transient states
	for(int r=0; r<nrow; r++)
		for(int c=0; c<nrow; c++)
		{
			if(r != c)
			{
				Eigen::VectorXi s1 = space.row(r);
				Eigen::VectorXi s2 = space.row(c);
				Eigen::Vector3i mig_popIDs = canMigrate(s1,s2,npops,totalN);
				// REMOVE
				//std::cout << "s1 = " << s1 <<"\n"
				//		<<  "s2 = " << s2 <<"\n"
				//		<<"mig_popIDs = " << mig_popIDs <<"\n";

				if(mig_popIDs(0) < 0)
					rateMatrix(r,c) = 0;
				else
				{
					double migrate = find_migRate((unsigned) mig_popIDs(0)+1, (unsigned) mig_popIDs(1)+1);
					// REMOVE
					//std::cout << "mig_popIDs " << mig_popIDs <<"\n";
					//std::cout << "migrate = " << migrate <<"\n";

					rateMatrix(r,c) = mig_popIDs(2)*migrate;

				//int noL_mig = space(mig_popIDs(0))
				}
			}
		}

	// absorbing state
	Eigen::VectorXd popSizes(npops);
	popSizes.setZero();
	find_popSizeTips(popSizes);
	// REMOVE
	// std::cout << "popSizes = " << popSizes <<"\n";


	int nLineages = 0;
	for(int r=0; r<nrow; r++)
	{
		Eigen::VectorXi stateAtLastPop = totalN;
		for(unsigned int i=0; i< (unsigned) (npops-1); i++)
		{
			stateAtLastPop -= space.row(r);

			nLineages = space.block(r,i*npops, 1, npops).sum();
			if(nLineages >=2 )
			  {
			    // modified by YC 6/10/2014
			    rateMatrix(r,nrow) += nLineages*(nLineages-1)/popSizes(i);
			    
			    // old version
			    // rateMatrix(r,nrow) += nLineages*(nLineages-1)*2/popSizes(i);
			  }

			// REMOVE - old version
			//for(unsigned int j=0; j<npops; j++)
			//{
			//	if(space(r,i*npops+j) >= 2)
			//		rateMatrix(r,nrow) += space(r,i*npops+j)*(space(r,i*npops+j)-1)*2/popSizes(i);
			//}
		}
		nLineages = stateAtLastPop.sum();
		if(nLineages >=2 )
		  {
		    // modified by YC 6/10/2014
		    rateMatrix(r,nrow) += nLineages*(nLineages-1)/popSizes(npops-1);
		    // old version
		    // rateMatrix(r,nrow) += nLineages*(nLineages-1)*2/popSizes(npops-1);
		  }

		// REMOVE - old version
		//for(unsigned int j=0; j<npops; j++)
		//	if(stateAtLastPop(j) >= 2)
		//		rateMatrix(r,nrow) += stateAtLastPop(j)*(stateAtLastPop(j)-1)*2/popSizes(npops-1);

		rateMatrix(r,r) = -rateMatrix.row(r).sum();
	}



	// REMOVE
	//std::cout << "Rate Matrix is... \n" << rateMatrix <<"\n";

	return rateMatrix;
}

/***
 * Remove the first coalescent node.
 */
void node::truncateFirstCoalNode(int newtipID)
{
	if(age == 0.0 && size_tree() == 2)
	{
		desc[0]->deleteCoalTree();
		desc[1]->deleteCoalTree();
		delete desc[0];
		delete desc[1];
		desc[0] = 0; desc[1] = 0;
		isTip = 1;
		tipID = newtipID;
	}
	else if(isTip == 0)
	{
		desc[0]->truncateFirstCoalNode(newtipID);
		desc[1]->truncateFirstCoalNode(newtipID);
	}
	return;
}

/***
 * Update popIDs of each node, since the populations lineages belong to may change before the speciation.
 * Moreover, the part of the population tree after the speciation will be truncated out,
 * so the popIDs should be relabelled (start from 1).
 */
void node::updatePopIDs(popTree* poptree, double eventTime)
{
	if(age == 0.0)
	{
		int found = 0;
		popID = poptree->getID_beforeSplit(popID, eventTime, found) - 2;
	}
	else if(isTip == 0)
	{
		popID -= 2;
		desc[0]->updatePopIDs(poptree, eventTime);
		desc[1]->updatePopIDs(poptree, eventTime);
	}

	return;
}

/***
 * @param time the smallest coalescent time on the given tree
 */
void node::truncate(int &numNewTips, double time)
{
	if(age > time && !isTip)
	{
		age -= time;
		desc[0]->truncate(numNewTips, time);
		desc[1]->truncate(numNewTips, time);
	}
	else if(age < time)
	{
		popID = 1; // initialize the population ID
		numNewTips++;

		if(!isTip)
			std::cout << "Error in node::truncate(double time).\n"
					"The given time may not be the smallest coalesnt time.\n";
	}
	else
	{
		// Note that the first coalesced node is not pruned out yet.

		popID = 1; // initialize the population ID
		age = 0;
		numNewTips++;

		if(size_tree()!=2)
			std::cout << "Error in node::truncate(double time).\n"
					"The given time may not be the smallest coalescent time or "
					"the given tree may be a binary tree.\n";
	}
}

/***
 * If eventID == 0, then it means speciation. PopID should be updated.
 */
void popTree::truncate_updateIDs(double time, int eventID)
{
	if(age > time && isTip==0)
	{
		age -= time;
		if(eventID == 0)
			popID -= 2;
		if(desc[0]->isTip == 0)
			desc[0]->truncate_updateIDs(time,eventID);
		if(desc[1]->isTip == 0)
			desc[1]->truncate_updateIDs(time,eventID);
	}
	else if(age == time)
	{
		// the first splitting node is truncated here.
		age = 0.0;
		isTip = 1;
		desc[0] = 0;
		desc[1] = 0;
		if(eventID == 0)
			popID -= 2;
	}
	else if(age == 0.0 && eventID ==0)
		popID -= 2;

}

int getID(Eigen::MatrixXi space, Eigen::VectorXi state)
{
	int id = -1;
	int idx_r = 0;
	while(id<0 && idx_r <space.rows())
	{
		if(state == space.row(idx_r).transpose())
		{
			id = idx_r;
			// REMOVe
			// std::cout << space.row(idx_r) <<"\n";
		}
		idx_r++;
	}
	return id;
}

int node::moreNextState(int nPops)
{
	// REMOVE
	//print_coaltree();

	int foundNextState = 0;

	if(isTip || (age == 0.0 && size_tree()==2))
	{
		// REMOVE
		//std::cout << "The original popID = " << popID <<"\n";

		if(popID < nPops)
		{
			foundNextState = 1;
			popID++;

			// REMOVE
			//std::cout << "The new popID = " << popID <<"\n";
		}
	}
	else if(!isTip && age >0.0)
	{
		foundNextState = desc[0]->moreNextState(nPops);
		if(!foundNextState)
			foundNextState = desc[1]->moreNextState(nPops);
	}

	return foundNextState;
}


Eigen::MatrixXi getNextPopAssignments(int numNewTips, int nPops)
{
	int nrow = std::pow(nPops,numNewTips);
	Eigen::MatrixXi popAssignment(nrow,numNewTips);

	//std::cout << "nrow = " << nrow <<"\n";

	if(numNewTips == 1)
	{
		for(int t1=1; t1 <= nPops; t1++)
			popAssignment(t1-1,0) = t1;
	}
	else if(numNewTips == 2)
	{
		for(int t1=1; t1 <= nPops; t1++)
			for(int t2=1; t2 <= nPops; t2++)
			{
				popAssignment.row((t1-1)*nPops+t2-1) << t1, t2;
				//std::cout << "t1 = " <<t1 <<"; t2=" <<t2 <<"; (t1-1)*nPops+t2-1=" <<(t1-1)*nPops+t2-1 <<"\n";
			}
	}
	else if(numNewTips > 2)
	{
		Eigen::MatrixXi subPopAssgn = getNextPopAssignments(numNewTips-1,nPops);
		int sub_nr = subPopAssgn.rows();
		for(int t1=1; t1 <= nPops; t1++)
		{
			popAssignment.block((t1-1)*sub_nr,1,sub_nr,numNewTips-1) = subPopAssgn;
			for(int i=0; i<sub_nr; i++)
				popAssignment((t1-1)*sub_nr+i,0) = t1;
		}
	}

	// REMOVE
	// std::cout << "In getNextPopAssignments():\n"
	//		<< "numNewTips = " << numNewTips <<", nPops = " << nPops
	//		<<"\npopAssignment:\n" << popAssignment <<"\n";


	return popAssignment;
}

/***
 *  Assign population IDs to tips and return the ID of the population where the first coalescent event happens
 *
 */
int node::assignPopulation(Eigen::VectorXi popAssign, int &idx, Eigen::VectorXi &state, int nPops)
{
	//REMOVE
	//std::cout << "In node::assignPopulation().\n Tree node is ";
	//print_coaltree();

	int popID_firstCoal = 0;

	if(isTip == 1)
	{
		popID = popAssign(idx);
		idx++;
		if(popID < nPops)
			state(label*(popID-1)+label)++;

		//print_coaltree();
		//std::cout << "popID = " <<popID<< "; state is " << state.transpose() << "\n";
	}
	else if(age == 0.0 && size_tree() ==2)
	{
		popID = popAssign(idx);
		idx++;
		if(popID < nPops)
		{
			state(desc[0]->label*(popID-1)+desc[0]->label) ++;
			state(desc[1]->label*(popID-1)+desc[1]->label) ++;
		}

		popID_firstCoal = popID;

		//print_coaltree();
		//std::cout << "popID = " <<popID<< "; state is " << state.transpose() << "\n";

	}
	else if(!isTip && age>0.0)
	{
		int tempID = 0;
		tempID = desc[0]->assignPopulation(popAssign, idx, state, nPops);
		if(tempID > 0)
			popID_firstCoal = tempID;
		tempID = desc[1]->assignPopulation(popAssign, idx, state, nPops);
		if(tempID > 0)
			popID_firstCoal = tempID;
	}

	return popID_firstCoal;
}

void node::newLabels()
{
	if(isTip == 1)
		label = popID -1;
	else if(isTip == 0)
	{
		desc[0]->newLabels();
		desc[1]->newLabels();
	}
}

double popTree::findPopSize(int popid, int &found)
{
	double popSize = 0.0;
	if(popID == (unsigned) popid)
	{
		 popSize = populationSize;
		 found = 1;
	}
	else if(isTip == 0)
	{
		if(found == 0 )
			popSize = desc[0]->findPopSize(popid,found);
		if(found == 0)
			popSize = desc[1]->findPopSize(popid,found);
	}
	return popSize;
}

/***
 *  Compute the conditional probability	of a coalescent tree in a single population
 */
double node::conditionalPr_singlePop(double popSize)
{
	//std::cout << "In node::conditionalPr_singlePop\n";
	//std::cout << "popSize = " << popSize <<"\n";
	//std::cout << "Tree is ";
	//print_coaltree();

	if(totalCoalRate == 0.0)
		compute_totalCoalescentRate();
	int nLineages = size_tree();

	double prob = std::pow(2/popSize,nLineages-1) * std::exp(-2*totalCoalRate/popSize);

	return prob;
}

int getNumNextOrignalStates(Eigen::VectorXi nextLumpedState, Eigen::VectorXi totalN)
{
	int noStates = 1;
	int npops = totalN.size();
	for(int id_pop=0; id_pop< npops-1; id_pop++)
	{
		for(int id_kind =0; id_kind< npops; id_kind++)
		{
			int n = totalN(id_kind);
			int k = nextLumpedState(id_pop*npops+id_kind);
			if(n>0 && k>0 && n!=k)
				noStates *= choose(n,k);
		}
		totalN -= nextLumpedState.segment(id_pop*npops,npops);
	}
	return noStates;
}

double compute_conditionalProb_recursion(node* coaltree, popTree* poptree, Eigen::VectorXi crrState)
{
	// REMOVE
	//std::cout << "\nIn compute_conditionalProb_recursion()\n";
	//std::cout << "The input coalescent tree is "; coaltree->print_coaltree();
	//std::cout << "The input population tree is "; poptree->print_poptree();
	//std::cout << "The initial state is " << crrState.transpose() <<"\n";
	//std::cout << "Population size is " << poptree->size() <<"\n";
	//coaltree->print_labelsPopIDs();

	double condlP = 0.0;
	int nPops = poptree->size();

	if(nPops == 1)
	{
		condlP = coaltree->conditionalPr_singlePop(poptree->get_popSize());
	}
	else
	{

		//---- Get the space of lumped states -----//
		Eigen::VectorXi totalN = coaltree->get_totalNumEachKind(poptree->size());
		Eigen::MatrixXi space = getSpace(totalN,poptree->size());

		//--- get the current state id ----//
		int id_crrState = getID(space,crrState);

		// REMOVE
		// std::cout << "the state space is \n" << space <<"\n";

		//---- Get the transition rate matrix ----//
		Eigen::MatrixXd rateMat = poptree->getTransitionRateMat(space, totalN);
		//---- Matrix decomposition -----//
		Eigen::EigenSolver<Eigen::MatrixXd> svd;
		svd.compute(rateMat);
		// eigen value and eigen vector matrices
		Eigen::MatrixXcd D = svd.eigenvalues().asDiagonal();
		Eigen::MatrixXcd V = svd.eigenvectors();
		Eigen::MatrixXcd V_inv = svd.eigenvectors().inverse();
		// Checking eigen values and eigen vectors
		/*
		Eigen::MatrixXcd I = V * V_inv;
		if(abs(I.sum().real() - space.rows()-1) > 3.0e-7)
		  {
			poptree->print_allPopTree();
			std::cout << "Eigen vectors may not be accurate. That is, V*V^(-1) != I.\n"
				<< "The sum of all elements in V*V^(-1) is " << I.sum().real()
				<< "\nThe sum of absolute differences is "<< abs(I.sum().real() - space.rows()+1) <<".\n";
		  }
		Eigen::MatrixXcd rateMat2 = V * D * V_inv;
		if(abs(rateMat2.sum().real() - rateMat.sum()) > 3.0e-7)
		  {
			poptree->print_allPopTree();
			std::cout << "Eigen values and eigen vectors may not be accurate. rateMat != V*D*V^(-1).\n"
				<<"The absolute difference of the sums of elements in rateMat and V*D*V^(-1) is "
					<< abs(rateMat2.sum().real() - rateMat.sum()) <<"\n";
		  }
		*/

		// ---- get the smallest event time (either coalescent or speciation)---//
		double coalT = coaltree->get_minCoalTime(-1.0);
		double splittingT = poptree->get_minSplittingTime(-1.0);
		double eventTime = 0.0;
		int coalesced = 0; // Indicator. 0 if the event is a speciation; 1 if the event is a coalescence.
		if(coalT > splittingT)
		{
			eventTime = splittingT;
			coalesced = 0;
		}
		else
		{
			eventTime = coalT;
			coalesced = 1;
		}

		//REMOVE
		//std::cout << "the smallest coalescent time is " << coalT
		//		<<" and the smallest splitting time is " << splittingT <<"\n";

		//--- Get exp(D*t) ---//
		Eigen::MatrixXcd expD = D*eventTime;
		for(unsigned int i=0; i<D.rows(); i++)
		{
			complex<double> a = expD(i,i);
			expD(i,i) = exp(a);
		}


		//--- Truncate a coalescent tree and assign possible labels to the new tips ---//
		int no_newTips = 0;
		node* shorterTree = new node;
		shorterTree->deepCopy_root(coaltree);
		shorterTree->truncate(no_newTips, eventTime);

		// --- List of possible population assignments --- //
		Eigen::MatrixXi popAssign_nextState;
		popAssign_nextState = getNextPopAssignments(no_newTips, nPops);

		// REMOVE
		//std::cout << "the number of new tips is " << no_newTips <<"\n";
		//std::cout << popAssign_nextState <<"\n";

		//std::cout << "coalesced = " << coalesced <<"; nPops = " << nPops <<"\n";

		// --- Recursively compute the conditional distribution of a tree ----//
		if(coalesced == 0 && nPops == 2)
		{
			double transitionP = 0.0;
			for(int i=0; i<popAssign_nextState.rows(); i++)
			{
				// Assign population ID to the new tips
				Eigen::VectorXi popIDs(no_newTips);
				for(int j=0; j< no_newTips; j++)
					popIDs(j) = popAssign_nextState(i,j);
				Eigen::VectorXi nextState(crrState.size());
				nextState.setZero();
				int idx =0;
				int popID_firstCoal = shorterTree->assignPopulation(popIDs, idx, nextState, nPops);
				// int foundThePop = 0;
				// double popSize_1stCoal = poptree->findPopSize(popID_firstCoal,foundThePop);

				// REMOVE
				//std::cout << "The 1st coal event occurred in pop with ID = "
				//		<< popID_firstCoal << " and size = " << popSize_1stCoal <<"\n";

				// get the next state id
				int id_nextState = getID(space,nextState);

				// compute the size of inverse image of the next lumped state
				int no_nextOriginalStates = getNumNextOrignalStates(nextState,totalN);

				// REMOVE
				// std::cout << "New lumped state is \n";
				// coaltree->print_labelsPopIDs();
				// std::cout <<"The number of the next possible original state is " <<no_nextOriginalStates <<"\n";

				// get the transition probability
				complex<double> prob = V.row(id_crrState) * expD * V_inv.col(id_nextState);
				transitionP += prob.real()/  no_nextOriginalStates;

				// REMOVE
				//std::cout <<"The transition probability is "<< transitionP <<"\n";

				//if(coalesced) // coalescent event
				//transitionP *= 2/popSize_1stCoal;
			}

			// ----  Compute the conditional probability of the truncated tree ----//
			double prob_subtree = 1.0;
			if(no_newTips >= 2)
			{
				// get the remaining coalescent tree
				node* subtree = new node;
				subtree->deepCopy_root(shorterTree);
				int maxid = 0;
				shorterTree->maxTipID(maxid);
				if(coalesced)
					subtree->truncateFirstCoalNode(maxid+1);
				else
					subtree->updatePopIDs(poptree,eventTime);
				// get the new labels on new tips according to their population IDs
				//subtree->newLabels();
				// get the new initial state
				//Eigen::VectorXi newInitialState(nPops*(nPops-1));
				//newInitialState.setZero();
				//subtree->get_tipState(newInitialState,nPops);

				// get the truncated pop tree
				//popTree* subPoptree = poptree->deepCopy_root();
				//subPoptree->truncate_updateIDs(eventTime, coalesced);


				// REMOVE
				//std::cout << "Fully truncated coalescent tree is ";
				//subtree->print_coaltree();
				//std::cout << "Fully truncated poptree tree is ";
				//subPoptree->print_poptree();
				//std::cout << "New initial state is " << newInitialState.transpose() <<"\n";

					// get the conditional probability of the truncated tree
				prob_subtree = subtree->conditionalPr_singlePop(poptree->get_popSize());
				// compute_conditionalProb_recursion(subtree,subPoptree,newInitialState);
				condlP += transitionP * prob_subtree;
				subtree->deleteCoalTree();
				delete subtree;
			}


		}
		else
		{
			for(int i=0; i<popAssign_nextState.rows(); i++)
			{
				// Assign population ID to the new tips
				Eigen::VectorXi popIDs(no_newTips);
				for(int j=0; j< no_newTips; j++)
					popIDs(j) = popAssign_nextState(i,j);
				Eigen::VectorXi nextState(crrState.size());
				nextState.setZero();
				int idx =0;
				int popID_firstCoal = shorterTree->assignPopulation(popIDs, idx, nextState, nPops);
				// Find the population size where the coalescent happened.
				double popSize_1stCoal = 0.0;
				if(coalesced == 1)
				{
					int foundThePop = 0;
					popSize_1stCoal = poptree->findPopSize(popID_firstCoal,foundThePop);
				}


				// get the next state id
				int id_nextState = getID(space,nextState);

				// compute the size of inverse image of the next lumped state
				int no_nextOriginalStates = getNumNextOrignalStates(nextState,totalN);


				// get the transition probability
				complex<double> prob = V.row(id_crrState) * expD * V_inv.col(id_nextState);
				double transitionP =prob.real()/  no_nextOriginalStates;
				if(coalesced == 1) // coalescent event
					transitionP *= 2/popSize_1stCoal;


				// ----  Compute the conditional probability of the truncated tree ----//
				double prob_subtree = 1.0;
				if(no_newTips >= 2)
				{
					// get the remaining coalescent tree
					// FIXME

					//node* subtree = new node;
					//subtree->deepCopy_root(subtree,shorterTree);
					//node* subtree = shorterTree->deepCopy_root();


					node subtree_obj;
					subtree_obj.deepCopy_obj(shorterTree);
					subtree_obj.desc[0]->par = &subtree_obj;
					subtree_obj.desc[1]->par = &subtree_obj;
					node* subtree = &subtree_obj;
					//subtree->print_coaltree();

					int maxid = 0;
					shorterTree->maxTipID(maxid);
					if(coalesced)
						subtree->truncateFirstCoalNode(maxid+1);
					else
						subtree->updatePopIDs(poptree,eventTime);
					// get the new labels on new tips according to their population IDs
					subtree->newLabels();
					// get the new initial state
					Eigen::VectorXi newInitialState(nPops*(nPops-1));
					newInitialState.setZero();
					subtree->get_tipState(newInitialState,nPops);

					// get the truncated pop tree
					popTree* subPoptree = poptree->deepCopy_root();
					subPoptree->truncate_updateIDs(eventTime, coalesced);


					// REMOVE
					// std::cout << "Fully truncated coalescent tree is ";
					// subtree->print_coaltree();
					// std::cout << "Fully truncated poptree tree is ";
					// subPoptree->print_poptree();
					// std::cout << "New initial state is " << newInitialState.transpose() <<"\n";

					// get the conditional probability of the truncated tree
					prob_subtree = compute_conditionalProb_recursion(subtree,subPoptree,newInitialState);
					subtree->deleteCoalTree();
					//delete subtree;
					subPoptree->deletePopTree();
					delete subPoptree;

				}
				condlP += transitionP * prob_subtree;
			}


			// REMOVE
			// std::cout << "population assignment is " << popAssign_nextState.row(i) <<"\n";
			// std::cout << "next state = " << nextState.transpose()
			// 		<< "\n nextStateID = " << id_nextState <<"\n";
			// std::cout <<  "conditional Prob = " << condlP <<"\n";

		}
		shorterTree->deleteCoalTree();
		delete shorterTree;
	}

	// REMOVE
	//std::cout << "\nIn compute_conditionalProb_recursion()\n";
	//coaltree->print_coaltree();
	//std::cout << "condlProb = " << condlP <<"\n\n";

	return condlP;
}




double compute_conditionalProb_recursion_noTreeCopy(node* coaltree, popTree* poptree, Eigen::VectorXi crrState)
{
	// REMOVE
	//std::cout << "\nIn compute_conditionalProb_recursion()\n";
	//std::cout << "The input coalescent tree is "; coaltree->print_coaltree();
	//std::cout << "The input population tree is "; poptree->print_poptree();
	//std::cout << "The initial state is " << crrState.transpose() <<"\n";
	//std::cout << "Population size is " << poptree->size() <<"\n";
	//coaltree->print_labelsPopIDs();

	double condlP = 0.0;
	int nPops = poptree->size();

	if(nPops == 1)
	{
		condlP = coaltree->conditionalPr_singlePop(poptree->get_popSize());
	}
	else
	{

		//---- Get the space of lumped states -----//
		Eigen::VectorXi totalN = coaltree->get_totalNumEachKind(poptree->size());
		Eigen::MatrixXi space = getSpace(totalN,poptree->size());

		//--- get the current state id ----//
		int id_crrState = getID(space,crrState);

		// REMOVE
		// std::cout << "the state space is \n" << space <<"\n";

		//---- Get the transition rate matrix ----//
		Eigen::MatrixXd rateMat = poptree->getTransitionRateMat(space, totalN);
		//---- Matrix decomposition -----//
		Eigen::EigenSolver<Eigen::MatrixXd> svd;
		svd.compute(rateMat);
		// eigen value and eigen vector matrices
		Eigen::MatrixXcd D = svd.eigenvalues().asDiagonal();
		Eigen::MatrixXcd V = svd.eigenvectors();
		Eigen::MatrixXcd V_inv = svd.eigenvectors().inverse();
		// Checking eigen values and eigen vectors
		Eigen::MatrixXcd I = V * V_inv;
		if(abs(I.sum().real() - space.rows()-1) > 3.0e-7)
		  {
		    poptree->print_allPopTree();
		    std::cout << "Eigen vectors may not be accurate. That is, V*V^(-1) != I.\n"
			      << "The sum of all elements in V*V^(-1) is " << I.sum().real()
			      << "\nThe sum of absolute differences is "<< abs(I.sum().real() - space.rows()+1) <<".\n";
		  }
		Eigen::MatrixXcd rateMat2 = V * D * V_inv;
		if(abs(rateMat2.sum().real() - rateMat.sum()) > 3.0e-7)
		  {
		    poptree->print_allPopTree();
		    std::cout << "Eigen values and eigen vectors may not be accurate. rateMat != V*D*V^(-1).\n"
			      <<"The absolute difference of the sums of elements in rateMat and V*D*V^(-1) is "
			      << abs(rateMat2.sum().real() - rateMat.sum()) <<"\n";
		  }

		// ---- get the smallest event time (either coalescent or speciation)---//
		double coalT = coaltree->get_minCoalTime(-1.0);
		double splittingT = poptree->get_minSplittingTime(-1.0);
		double eventTime = 0.0;
		int coalesced = 0; // Indicator. 0 if the event is a speciation; 1 if the event is a coalescence.
		if(coalT > splittingT)
		{
			eventTime = splittingT;
			coalesced = 0;
		}
		else
		{
			eventTime = coalT;
			coalesced = 1;
		}

		//REMOVE
		//std::cout << "the smallest coalescent time is " << coalT
		//		<<" and the smallest splitting time is " << splittingT <<"\n";

		//--- Get exp(D*t) ---//
		Eigen::MatrixXcd expD = D*eventTime;
		for(unsigned int i=0; i<D.rows(); i++)
		{
			complex<double> a = expD(i,i);
			expD(i,i) = exp(a);
		}


		//--- Truncate a coalescent tree and assign possible labels to the new tips ---//
		int no_newTips = 0;
		node* shorterTree = new node;
		shorterTree->deepCopy_root(coaltree);
		shorterTree->truncate(no_newTips, eventTime);

		// --- List of possible population assignments --- //
		Eigen::MatrixXi popAssign_nextState;
		popAssign_nextState = getNextPopAssignments(no_newTips, nPops);

		// REMOVE
		//std::cout << "the number of new tips is " << no_newTips <<"\n";
		//std::cout << popAssign_nextState <<"\n";

		//std::cout << "coalesced = " << coalesced <<"; nPops = " << nPops <<"\n";

		// --- Recursively compute the conditional distribution of a tree ----//
		if(coalesced == 0 && nPops == 2)
		{
			double transitionP = 0.0;
			for(int i=0; i<popAssign_nextState.rows(); i++)
			{
				// Assign population ID to the new tips
				Eigen::VectorXi popIDs(no_newTips);
				for(int j=0; j< no_newTips; j++)
					popIDs(j) = popAssign_nextState(i,j);
				Eigen::VectorXi nextState(crrState.size());
				nextState.setZero();
				int idx =0;
				int popID_firstCoal = shorterTree->assignPopulation(popIDs, idx, nextState, nPops);
				// int foundThePop = 0;
				// double popSize_1stCoal = poptree->findPopSize(popID_firstCoal,foundThePop);

				// REMOVE
				//std::cout << "The 1st coal event occurred in pop with ID = "
				//		<< popID_firstCoal << " and size = " << popSize_1stCoal <<"\n";

				// get the next state id
				int id_nextState = getID(space,nextState);

				// compute the size of inverse image of the next lumped state
				int no_nextOriginalStates = getNumNextOrignalStates(nextState,totalN);

				// REMOVE
				// std::cout << "New lumped state is \n";
				// coaltree->print_labelsPopIDs();
				// std::cout <<"The number of the next possible original state is " <<no_nextOriginalStates <<"\n";

				// get the transition probability
				complex<double> prob = V.row(id_crrState) * expD * V_inv.col(id_nextState);
				transitionP += prob.real()/  no_nextOriginalStates;

				// REMOVE
				//std::cout <<"The transition probability is "<< transitionP <<"\n";

				//if(coalesced) // coalescent event
				//transitionP *= 2/popSize_1stCoal;
			}

			// ----  Compute the conditional probability of the truncated tree ----//
			double prob_subtree = 1.0;
			if(no_newTips >= 2)
			{
				// get the remaining coalescent tree
				node* subtree = new node;
				subtree->deepCopy_root(shorterTree);
				int maxid = 0;
				shorterTree->maxTipID(maxid);
				if(coalesced)
					subtree->truncateFirstCoalNode(maxid+1);
				else
					subtree->updatePopIDs(poptree,eventTime);
				// get the new labels on new tips according to their population IDs
				//subtree->newLabels();
				// get the new initial state
				//Eigen::VectorXi newInitialState(nPops*(nPops-1));
				//newInitialState.setZero();
				//subtree->get_tipState(newInitialState,nPops);

				// get the truncated pop tree
				//popTree* subPoptree = poptree->deepCopy_root();
				//subPoptree->truncate_updateIDs(eventTime, coalesced);


				// REMOVE
				//std::cout << "Fully truncated coalescent tree is ";
				//subtree->print_coaltree();
				//std::cout << "Fully truncated poptree tree is ";
				//subPoptree->print_poptree();
				//std::cout << "New initial state is " << newInitialState.transpose() <<"\n";

					// get the conditional probability of the truncated tree
				prob_subtree = subtree->conditionalPr_singlePop(poptree->get_popSize());
				// compute_conditionalProb_recursion(subtree,subPoptree,newInitialState);
				condlP += transitionP * prob_subtree;
				subtree->deleteCoalTree();
				delete subtree;
			}


		}
		else
		{
			for(int i=0; i<popAssign_nextState.rows(); i++)
			{
				// Assign population ID to the new tips
				Eigen::VectorXi popIDs(no_newTips);
				for(int j=0; j< no_newTips; j++)
					popIDs(j) = popAssign_nextState(i,j);
				Eigen::VectorXi nextState(crrState.size());
				nextState.setZero();
				int idx =0;
				int popID_firstCoal = shorterTree->assignPopulation(popIDs, idx, nextState, nPops);
				// Find the population size where the coalescent happened.
				double popSize_1stCoal = 0.0;
				if(coalesced == 1)
				{
					int foundThePop = 0;
					popSize_1stCoal = poptree->findPopSize(popID_firstCoal,foundThePop);
				}


				// get the next state id
				int id_nextState = getID(space,nextState);

				// compute the size of inverse image of the next lumped state
				int no_nextOriginalStates = getNumNextOrignalStates(nextState,totalN);


				// get the transition probability
				complex<double> prob = V.row(id_crrState) * expD * V_inv.col(id_nextState);
				double transitionP =prob.real()/  no_nextOriginalStates;
				if(coalesced == 1) // coalescent event
					transitionP *= 2/popSize_1stCoal;


				// ----  Compute the conditional probability of the truncated tree ----//
				double prob_subtree = 1.0;
				if(no_newTips >= 2)
				{
					// get the remaining coalescent tree
					// FIXME

					//node* subtree = new node;
					//subtree->deepCopy_root(subtree,shorterTree);
					//node* subtree = shorterTree->deepCopy_root();


					node subtree_obj;
					subtree_obj.deepCopy_obj(shorterTree);
					subtree_obj.desc[0]->par = &subtree_obj;
					subtree_obj.desc[1]->par = &subtree_obj;
					node* subtree = &subtree_obj;
					//subtree->print_coaltree();

					int maxid = 0;
					shorterTree->maxTipID(maxid);
					if(coalesced)
						subtree->truncateFirstCoalNode(maxid+1);
					else
						subtree->updatePopIDs(poptree,eventTime);
					// get the new labels on new tips according to their population IDs
					subtree->newLabels();
					// get the new initial state
					Eigen::VectorXi newInitialState(nPops*(nPops-1));
					newInitialState.setZero();
					subtree->get_tipState(newInitialState,nPops);

					// get the truncated pop tree
					popTree* subPoptree = poptree->deepCopy_root();
					subPoptree->truncate_updateIDs(eventTime, coalesced);


					// REMOVE
					// std::cout << "Fully truncated coalescent tree is ";
					// subtree->print_coaltree();
					// std::cout << "Fully truncated poptree tree is ";
					// subPoptree->print_poptree();
					// std::cout << "New initial state is " << newInitialState.transpose() <<"\n";

					// get the conditional probability of the truncated tree
					prob_subtree = compute_conditionalProb_recursion(subtree,subPoptree,newInitialState);
					subtree->deleteCoalTree();
					//delete subtree;
					subPoptree->deletePopTree();
					delete subPoptree;

				}
				condlP += transitionP * prob_subtree;
			}


			// REMOVE
			// std::cout << "population assignment is " << popAssign_nextState.row(i) <<"\n";
			// std::cout << "next state = " << nextState.transpose()
			// 		<< "\n nextStateID = " << id_nextState <<"\n";
			// std::cout <<  "conditional Prob = " << condlP <<"\n";

		}
		shorterTree->deleteCoalTree();
		delete shorterTree;
	}

	return condlP;
}



// old version? - YC 6/5/2014
/***
 * Computing the conditional distribution of coalescent tree given a population tree
 */
double compute_conditionalProb(node* coaltree, popTree* poptree, locus lc)
{
	//---- Assign populations to tips and relabel tips according to population IDs----//
	coaltree->assignPopulations2Tips(lc);

	//---- Current lumped states ------//
	Eigen::VectorXi crrState(poptree->size()*(poptree->size()-1));
	crrState.setZero();
	coaltree->get_tipState(crrState,poptree->size());


	return compute_conditionalProb_recursion(coaltree, poptree, crrState);
}



void Chain::compute_observedStates_fromSubtrees(unsigned int id_listTopo, unsigned int id_coalEvent,
		std::vector<nodeSimple*> subtrees)
{
	// REMOVE
	//std::cout <<"\nIn Chain::compute_observedStates_fromSubtrees()\n";
	//std::cout << "id_listTopo = " << id_listTopo <<"; id_coalEvent = " << id_coalEvent <<"\n";
	//for(unsigned int i=0; i<subtrees.size(); i++)
	//{
	//	subtrees.at(i)->print_topo();
	//	std::cout <<"\n";
	//}

  unsigned int n_lineages = subtrees.size();
  unsigned int whereTheLargestRank = 0;
  unsigned int maxRank = 0;
  for(unsigned int i=0; i<n_lineages; i++)
    {
      // string state = subtrees.at(i)->convert2Newick_topo_root();
      unsigned int state = 0;
      if(subtrees.at(i)->get_isTip()==1)
	state = subtrees.at(i)->getPopID();
      else
	state = subtrees.at(i)->getRank();

      unsigned int sameState = 0;
      unsigned int insert_loc = 0;
      unsigned int count = 0;
      while(sameState==0 && count < states_observed.at(id_listTopo).at(id_coalEvent).size())
	{
			//if(state.compare(states_observed.at(id_listTopo).at(id_coalEvent).at(count)) == 0) // same string
	  if(state == states_observed.at(id_listTopo).at(id_coalEvent).at(count)) // same string
	    {
	      sameState = 1;
	      states_observed_freq.at(id_listTopo).at(id_coalEvent).at(count)++;
	    }
	  else if(state > states_observed.at(id_listTopo).at(id_coalEvent).at(count))
	    {
	      insert_loc = count+1;
	    }
	  count++;
	}
      if(sameState == 0)
	{
	  if(insert_loc >=  states_observed.at(id_listTopo).at(id_coalEvent).size())
	    {
	      states_observed.at(id_listTopo).at(id_coalEvent).push_back(state);
	      states_observed_freq.at(id_listTopo).at(id_coalEvent).push_back(1);
	    }
	  else
	    {
	      //std::cout << "New state is " << state << " and insert_loc is " << insert_loc <<"\n";
	      
	      std::vector<unsigned int>::iterator iter1 = states_observed.at(id_listTopo).at(id_coalEvent).begin();
	      std::vector<unsigned int>::iterator iter2 = states_observed_freq.at(id_listTopo).at(id_coalEvent).begin();
	      states_observed.at(id_listTopo).at(id_coalEvent).insert(iter1+insert_loc, state);
	      states_observed_freq.at(id_listTopo).at(id_coalEvent).insert(iter2+insert_loc, 1);
	    }
	}
      
      
      if(maxRank < subtrees.at(i)->getRank())
	{
	  maxRank = subtrees.at(i)->getRank();
	  whereTheLargestRank = i;
	}
    }
  
	//REMOVE
	/*
	std::cout <<"\nIn Chain::compute_observedStates_fromSubtrees()\n";
	std::cout << "state freq\n";
	for(unsigned int i=0; i<states_observed.at(id_listTopo).at(id_coalEvent).size(); i++)
	{
		std::cout << states_observed.at(id_listTopo).at(id_coalEvent).at(i)
				<<" " <<states_observed_freq.at(id_listTopo).at(id_coalEvent).at(i)
				<<"\n";
	}
	*/

	// REMOVE
	//std::cout << "Largest rank is " << maxRank <<" and the location is " << whereTheLargestRank <<"\n";

	if(maxRank != 0)
	{
		nodeSimple* tr1 = subtrees.at(whereTheLargestRank)->getFirstChild();
		nodeSimple* tr2 = subtrees.at(whereTheLargestRank)->getSecondChild();
		subtrees.at(whereTheLargestRank) = tr1;
		std::vector<nodeSimple*>::iterator iter = subtrees.begin();
		subtrees.insert(iter+whereTheLargestRank+1,tr2);
		compute_observedStates_fromSubtrees(id_listTopo, id_coalEvent-1, subtrees);
	}

}


void Chain::compute_observedStates_fromTopo()
{
	// REMOVE
  // std::cout << "In Chain::compute_observedStates_fromTopo()\n";

  unsigned int numUniqTopo = list_trees.size();
  states_observed.resize(numUniqTopo);
  states_observed_freq.resize(numUniqTopo);
  nKinds_lineages.resize(numUniqTopo);
  
  for(unsigned int i=0; i<numUniqTopo; i++)
    {
		//REMOVE
		//list_trees.at(i)->print_topo();
      
      unsigned int nGeneCopies = list_trees.at(i)->getSize();
      std::vector<nodeSimple*> subtrees;
      subtrees.push_back(list_trees.at(i)->getFirstChild());
      subtrees.push_back(list_trees.at(i)->getSecondChild());
      states_observed.at(i).resize(nGeneCopies-1);
      states_observed_freq.at(i).resize(nGeneCopies-1);
      nKinds_lineages.at(i).resize(nGeneCopies-1);

      compute_observedStates_fromSubtrees(i,nGeneCopies-2,subtrees);

      for(unsigned int j=0; j<nGeneCopies-1; j++)
	{
	  nKinds_lineages.at(i).at(j) = states_observed.at(i).at(j).size();
	}
    }	
  
  return;
}

/***
 * Compute the number of ways to assign a set of 'freq' elements
 * into 'nPops' populations.
 * Note that this number is different from Stirling numbers, where populations
 * are not distinguished.
 */
unsigned int Chain::compute_size_stateSpaces(unsigned int freq, unsigned int nPops)
{
	unsigned int size = 0;
	if(nPops == 2)
	{
		size = freq+1;
	}
	else if(nPops ==3)
	{
		size = (freq+1)*(freq+2)/2;
	}
	else
	{
		for(unsigned int i=0; i<=freq; i++)
		{
			size += compute_size_stateSpaces(freq-i, nPops-1);
		}
	}

	return size;
}

Eigen::MatrixXd Chain::compute_stateSpaces_recursion(std::vector<unsigned int> freq, unsigned int nPops) //, unsigned int nr, unsigned int nc)
{
	// REMOVE
	//std::cout << "In Chain::compute_stateSpaces_recursion()\n";
	//std::cout << "nPops = " << nPops <<" and freq is \n";
	//for(unsigned int i=0; i< freq.size(); i++)
//		std::cout << freq.at(i) << " ";
	//std::cout << "\n";

  Eigen::MatrixXd mat;
  unsigned int nKinds = freq.size();
  if(nPops == 2)
    {
      if(nKinds == 1)
	{
	  unsigned int sub_size = freq.at(0)+1;
	  mat.resize(sub_size,nPops);
	  for(unsigned int i=0; i<sub_size; i++)
	    {
	      mat(i,0) = freq.at(0)-i;
	      mat(i,1) = i;
	    }
	}
      else if(nKinds > 1)
	{
	  std::vector<unsigned int> sub_freq = freq;
	  sub_freq.erase(sub_freq.begin());
	  Eigen::MatrixXd subMat;
	  subMat = compute_stateSpaces_recursion(sub_freq,nPops);
	  unsigned int sub_nr = subMat.rows();
	  unsigned int sub_size = freq.at(0)+1;
	  mat.resize(sub_size*sub_nr,nKinds*nPops);
	  mat.setZero();
	  for(unsigned int i=0; i<sub_size; i++)
	    {
	      mat.block(i*sub_nr,0,sub_nr,1).fill(freq.at(0)-i);
	      mat.block(i*sub_nr,nKinds,sub_nr,1).fill(i);
	      
	      mat.block(i*sub_nr,1,sub_nr,nKinds-1) = subMat.block(0,0,sub_nr,nKinds-1);
	      mat.block(i*sub_nr,nKinds+1,sub_nr, nKinds-1) = subMat.block(0,nKinds-1,sub_nr,nKinds-1);
	    }
	}
    }
  
  // REMOVE
  //std::cout << mat <<"\n\n";
  
  return mat;
}

void Chain::compute_stateSpaces(unsigned int nPops)
{
  unsigned int numUniqTopo = list_trees.size();
  stateSpaces.resize(numUniqTopo);
  for(unsigned int i=0; i<numUniqTopo; i++)
    {
      unsigned int nGeneCopies = list_trees.at(i)->getSize();
      stateSpaces.at(i).resize(nGeneCopies-1);
      for(unsigned int j=0; j<nGeneCopies-1; j++)
	{
	  if(i > 0 && j==0) /// the state spaces for the first time period (until the first coalescent event) are the same across all the trees
	    stateSpaces.at(i).at(0) = stateSpaces.at(0).at(0);
	  else
	    {
	      stateSpaces.at(i).at(j) = compute_stateSpaces_recursion(states_observed_freq.at(i).at(j), nPops);
	    }
	  
	  // REMOVE
	  // std::cout << stateSpaces.at(i).at(j) <<"\n\n";
	}
    }
  
  
  
  return;
}

/***
 * Find the initial state of a tree topology.
 * Note that the initial state should be the same across tree topologies.
 */
unsigned int Chain::find_initialState(unsigned int nPops)
{
  unsigned int id_initial = 0;

	// REMOVE
	/*
	std::cout << "in Chain::find_initialState()\n";
	//std::cout << stateSpaces.at(0).at(0) <<"\n\n";
	std::cout << "state freq\n";
	for(unsigned int i=0; i<states_observed_freq.at(0).at(0).size(); i++)
	{
		std::cout << states_observed.at(0).at(0).at(i) << " ";
		std::cout << states_observed_freq.at(0).at(0).at(i) << "\n";
	}
	std::cout << "\n";
	*/

  unsigned int nKinds = states_observed_freq.at(0).at(0).size();
  /*
  if(nKinds != nPops)
    {
      std::cout << "\n *** Error in Chain::find_initialState() ***\n"
	"The number of populations (nPops) and the number of kinds of lineages (nKinds)"
	"should be the same, but nPops = " << nPops << " and nKinds = " << nKinds <<".\n\n";
    }
  */
  Eigen::MatrixXd initialState(1,nKinds*nPops);
  initialState.setZero();
  for(unsigned int i=0; i<nKinds; i++)
    initialState(0,i*nKinds+i) = states_observed_freq.at(0).at(0).at(i);
  


  unsigned int found = 0;
  unsigned int count = 0;
  while(count<stateSpaces.at(0).at(0).rows() && found ==0 )
    {
      if(stateSpaces.at(0).at(0).row(count) == initialState)
	{
	  id_initial = count;
	  found = 1;
	  
	  // REMOVE
	  //std::cout << stateSpaces.at(0).at(0).row(count) <<"\n";
	  //std::cout << "id_initial = " << id_initial <<"\n";
	}
      count++;
    }
  
  // REMOVE
  /*
  std::cout << "The initial state is ";
  std::cout << initialState << "\n";
  std::cout << "id_initial = " << id_initial <<"\n";
  */
  
  return id_initial;
}

Eigen::MatrixXd Chain::find_minimumFreq_forNextTimePeriod(unsigned int id_samples,unsigned int  id_period)
{
	unsigned int nKinds =states_observed_freq.at(id_samples).at(id_period).size();
	unsigned int nKinds_next =states_observed_freq.at(id_samples).at(id_period).size();

	Eigen::MatrixXd minimumFreq(1,nKinds);
	minimumFreq.setZero();

	unsigned int count_crr = 0;
	unsigned int count_next = 0;
	unsigned int found = 0;
	while(found < 2  && count_crr < nKinds && count_next < nKinds_next)
	{
		unsigned int lineage =states_observed.at(id_samples).at(id_period).at(count_crr);
		unsigned int lin_next =states_observed.at(id_samples).at(id_period+1).at(count_next);
		unsigned int freq = states_observed_freq.at(id_samples).at(id_period).at(count_crr);
		if(lineage == lin_next) // same kind of lineage
		{
			unsigned int freq_next =states_observed_freq.at(id_samples).at(id_period+1).at(count_next);
			if(freq-1 == freq_next)
			{
				// this lineage, lin, is coalesced
				minimumFreq(0,count_crr) = 1;
				found++;
			}
			else if(freq-2 == freq_next)
			{
				// two of this kinds of lineages, lin, are coalesced
				minimumFreq(0,count_crr) = 2;
				found = 2;
			}
			count_next++;
		}
		else if(freq == 1)
		{
			minimumFreq(0,count_crr) = 1;
			found++;
		}
		else if(freq == 2)
		{
			minimumFreq(0,count_crr) = 2;
			found = 2;
		}
		count_crr++;
	}
	return minimumFreq;
}

std::string Chain::find_nextCoalLineage(unsigned int id_samples,unsigned int  id_period)
{
	std::string lin;

	unsigned int nKinds =states_observed.at(id_samples).at(id_period).size();
	unsigned int next_nKinds =states_observed.at(id_samples).at(id_period+1).size();

	unsigned int found = 0;
	unsigned int count =0;
	while(found ==0 && count < nKinds && count <next_nKinds)
	{
		//std::string next_state = states_observed.at(id_samples).at(id_period+1).at(count);
		//std::string prev_state = states_observed.at(id_samples).at(id_period).at(count);
		unsigned int next_state = states_observed.at(id_samples).at(id_period+1).at(count);
		unsigned int prev_state = states_observed.at(id_samples).at(id_period).at(count);
		//if(next_state.compare(prev_state) == 0)
		if(next_state == prev_state)
		{
			unsigned int next_freq = states_observed_freq.at(id_samples).at(id_period+1).at(count);
			unsigned int prev_freq = states_observed_freq.at(id_samples).at(id_period).at(count);
			if(prev_freq+1 == next_freq)
			{
				found = 1;
				lin = next_state;
			}
		}
		// else if(next_state.compare(prev_state) != 0)
		else
		{
			found = 1;
			lin = next_state;
		}
		count++;
	}

	return lin;
}

void Chain::compute_possiblePaths(unsigned int nPops)
{

	unsigned int numUniqTopo = list_trees.size();
	possiblePaths.resize(numUniqTopo);
	for(unsigned int i=0; i<numUniqTopo; i++)
	{
		unsigned int nGeneCopies = list_trees.at(i)->getSize();
		possiblePaths.at(i).resize(nGeneCopies-1);
		for(unsigned int j=0; j<nGeneCopies-1; j++)
		{
			if(j==0)
				possiblePaths.at(i).at(j).resize(nPops+1);
			else
				possiblePaths.at(i).at(j).resize(nPops*2);
		}
	}

	unsigned int id_initialState = find_initialState(nPops);
	for(unsigned int i=0; i<numUniqTopo; i++)
	{
		// REMOVE
		//std::cout << "tree is ";
		//list_trees.at(i)-> print_topo();
		//std::cout << "id_initialState = " <<id_initialState << "\n";

		unsigned int nGeneCopies = list_trees.at(i)->getSize();
		for(unsigned int j=0; j<nGeneCopies-1; j++)
		{
			// std::cout << "j = " << j <<"\n\n";

			//std::cout << stateSpaces.at(i).at(j) <<"\n\n";


			if(j==0) // initial state of a tree topology
			{
				possiblePaths.at(i).at(j).at(0).push_back(id_initialState);
			}

			unsigned int nKinds = states_observed.at(i).at(j).size();

			// Final state of each time period
			if(j!= nGeneCopies-2) // there are more than two lineages
			{
				unsigned int next_nKinds =states_observed.at(i).at(j+1).size();

				// Find the minimum frequency of each kind of lineages so that the next coalescent event can happen
				Eigen::MatrixXd minimumFreq = find_minimumFreq_forNextTimePeriod(i,j);

				// std::string coalLin = find_nextCoalLineage(i,j);
				// unsigned int coalLin = states_observed.at(i).at(j+1).at(next_nKinds-1);

				// REMOVE
				// std::cout << "minimumFreq  = \n" << minimumFreq <<"\n";

				//Eigen::MatrixXd remainingFreq; // the number of lineages in each kind after remove the lineages to coalesce in the next time period

				unsigned int spaceSize = stateSpaces.at(i).at(j).rows();
				for(unsigned int id_state=0; id_state<spaceSize; id_state++)
				{
					for(unsigned int p=0; p<nPops; p++)
					{
						unsigned int coal_possible = 1;
						unsigned int count = 0;
						while(coal_possible == 1 && count <nKinds)
						{
							if(stateSpaces.at(i).at(j)(id_state,p*nKinds+count) < minimumFreq(count))
							{
								coal_possible = 0;
							}
							count++;
						}

						if(coal_possible == 1)
						{
							if(j==0)
								possiblePaths.at(i).at(j).at(1+p).push_back(id_state);
							else
								possiblePaths.at(i).at(j).at(nPops+p).push_back(id_state);

							// Find next initial state
							Eigen::MatrixXd nextState(1,next_nKinds*nPops);
							nextState.setZero();
							unsigned int count_crr = 0;
							for(unsigned int k=0; k< next_nKinds-1; k++)
							{
								unsigned int state_next = states_observed.at(i).at(j+1).at(k);
								unsigned int found_sameState = 0;
								while(found_sameState == 0 && count_crr < nKinds)
								{
									unsigned int state_crr = states_observed.at(i).at(j).at(count_crr);
									if(state_next == state_crr)
									{
										found_sameState = 1;
										for(unsigned int pp=0; pp<nPops; pp++)
										{
											if(pp == p)
												nextState(0,pp*next_nKinds+k) =
													stateSpaces.at(i).at(j)(id_state,pp*nKinds+ count_crr)
													- minimumFreq(0,count_crr);
											else
												nextState(0,pp*next_nKinds+k) =
													stateSpaces.at(i).at(j)(id_state,pp*nKinds+count_crr);
										}
									}
									count_crr++;
								}
							}
							for(unsigned int pp=0; pp<nPops; pp++)
							{
								if(pp == p)
									nextState(0,pp*next_nKinds+next_nKinds-1) =1;
							}

							// REMOVE
							//std::cout << "crr state is " <<	stateSpaces.at(i).at(j).row(id_state) <<"\n";
							//std::cout << "next state is " << nextState <<"\n";

							unsigned int id_nextState = 0;
							unsigned int found_id_nextState =0;
							while(found_id_nextState == 0 && id_nextState <stateSpaces.at(i).at(j+1).rows())
							{
								if(stateSpaces.at(i).at(j+1).row(id_nextState) == nextState)
								{
									found_id_nextState = 1;
									possiblePaths.at(i).at(j+1).at(p).push_back(id_nextState);
								}
								id_nextState++;
							}
						}
					}
				}

				// REMOVE
				/*
				std::cout << "minmumFreq = \n"<< minimumFreq <<"\n";
				std::cout << stateSpaces.at(i).at(j) <<"\n\n";
				for(unsigned int p=0; p<nPops; p++)
				{
					for(unsigned int k=0; k< possiblePaths.at(i).at(j).at(p).size(); k++)
						std::cout <<"p="<<p <<" k="<< k<< " id = "<< possiblePaths.at(i).at(j).at(p).at(k) <<"\n";
				}
				*/


			}
			else // there are two lineages
			{
				for(unsigned int id_state =0; id_state < stateSpaces.at(i).at(j).rows(); id_state++)
				{
					for(unsigned int p=0; p<nPops; p++)
					{
						if(stateSpaces.at(i).at(j).block(id_state,p*nKinds,1,nKinds).sum() == 2)
						{
							if(j==0)
								possiblePaths.at(i).at(j).at(1+p).push_back(id_state);
							else
								possiblePaths.at(i).at(j).at(nPops+p).push_back(id_state);
						}
					}
				}
			}


		}
	}

	// REMOVE
	std::cout << "*** In Chain::compute_possiblePaths() ***\n";
	for(unsigned int i=0; i<numUniqTopo; i++)
	  { 	    
	    for(unsigned int j=0; j<list_trees.at(i)->getSize()-1; j++)
	      {
		for(unsigned int p=0; p<possiblePaths.at(i).at(j).size(); p++)
		  {
		    for(unsigned int pp=0; pp<possiblePaths.at(i).at(j).at(p).size(); pp++)
		      {
			std::cout << possiblePaths.at(i).at(j).at(p).at(pp) <<" ";
		      }
		    std::cout <<"\n";
		  }
	      }
	  }
	std::cout << "*** End of Chain::compute_possiblePaths() ***\n";

	
	return;
}


void Chain::compute_numTrees(unsigned int nPops)
{
  unsigned int size =list_trees.size();
  numTrees.resize(size);
  for(unsigned int i=0; i<size; i++)
    {
      numTrees.at(i) = list_trees.at(i)->compute_nTrees_sameTopo(nPops);
    }
  
	// REMOVE
	/*
	std::cout << "*** In Chain::compute_numTrees() ***\n";
	for(unsigned int i=0; i<size; i++)
	  {
	    std::cout << numTrees.at(i) << " ";
	  }	
	std::cout << "\n*** End of Chain::compute_numTrees() ***\n";
	*/

	return;
}


void Chain::initializeLmode(IM im, unsigned int crrProcID, unsigned int nProcs)
{
  percentile_4upperBoundOfTMRCA =100; // this option is not used any longer - YC 6/23/2015

  n_MCMCgen = im.get_nSampledTrees(); 
  n_loci = im.get_nLoci();
  nTrees = im.get_nTrees2compute();
  lociInParallel = im.get_lociInParallel();
  TreesWithMaxP = im.get_TreesWithMaxP();
  multiLocusSpecific_mutationRate = im.get_multiLocusSpecific_mutationRate();
  nProcesses = nProcs;
  cpuID = crrProcID;

  if(lociInParallel==1) // CPUs handle different loci
    {
      nSubSample = nTrees;
      sampleID_start = 0;
      if(TreesWithMaxP==1)
	sampleID_end = n_MCMCgen-1;
      else
	sampleID_end = nTrees-1;
      if(nProcs > n_loci)
	{
	  std::cout << "\n\n *** Error in Chain::initializeLmode ***\n";
	  std::cout << "The number of CPUs (nProcs = "<< nProcs <<")"
		    << " should be smaller than or equal to the number "
		    << "of loci (n_loci = " <<n_loci <<")\n\n";
	}

      numSubLoci = static_cast<unsigned int>(n_loci/nProcs);
      unsigned int nRemainder = n_loci - numSubLoci*nProcs;   
      for(unsigned int i=0; i<nProcs; i++)
	{
	  if(crrProcID == i && nProcs-i <= nRemainder)
	    numSubLoci++;
	}   
      if(nProcs - crrProcID > nRemainder)
	locusID_start = numSubLoci*crrProcID;
      else
	locusID_start = (numSubLoci-1)*(nProcs-nRemainder)
	  + numSubLoci*(crrProcID-(nProcs-nRemainder));
      
      locusID_end = locusID_start + numSubLoci-1;
      std::cout << "On process " << crrProcID
      << ": reading " << nSubSample << " sample from locusID = " << locusID_start
		<< " to locusID = " << locusID_end
		<< " (the number of loci = " << numSubLoci << ")"
		<<"\n";
    }
  else // CPUs handle different samples
    {
      numSubLoci = n_loci;
      locusID_start =0; 
      locusID_end = n_loci-1;
      if(nProcs > nTrees)
	{
	  std::cout << "\n *** Error in Chain::initializeLmode ***\n";
	  std::cout << "The number of CPUs (nProcs = "<< nProcs <<")"
		    << " should be smaller than or equal to the number "
		    << "of tree samples per locus (nTrees = " <<nTrees
		    <<") to consider in L mode.\n\n";
	}
      nSubSample = static_cast<unsigned int>(nTrees/nProcs);
      unsigned int nRemainder = nTrees - nSubSample*nProcs;    
      for(unsigned int i=0; i<nProcs; i++)
	{
	  if(crrProcID == i && nProcs-i <= nRemainder)
	    nSubSample++;
	}  
      if(nProcs - crrProcID > nRemainder)
	sampleID_start = nSubSample*crrProcID;
      else
	sampleID_start = (nSubSample-1)*(nProcs-nRemainder)
	  + nSubSample*(crrProcID-(nProcs-nRemainder));
      
      sampleID_end = sampleID_start + nSubSample-1;

      std::cout << "On process " << crrProcID
		<< ": reading sample from id = " << sampleID_start
		<< " to id = " << sampleID_end
		<<"\n";
    }
  
  

  return;
}




void Chain::prepare_Lmode(popTree* poptree)
{
  // std::cout << "poptree->get_age() = " << poptree->get_age()<<"\n";
  compute_observedStates_fromTopo();
  if(poptree->get_age() > 0) // not a single population
    {
      compute_stateSpaces(poptree->size());
      //	compute_possiblePaths(poptree->size());
      compute_possiblePaths_nPossibleCoalEvents(poptree->size());
      compute_coefficientsOfElementsInTransitionRateMat(poptree->size());
    }
  compute_numTrees(poptree->size());
  compute_numSameCoalEvents_inAncPop();
}

// Computing the term of the number of same coalescent events in the ancestral population
void Chain::compute_numSameCoalEvents_inAncPop()
{
  // std::cout << "\nIn Chain::compute_numSameCoalEvents_inAncPop\n";

  unsigned int numUniqTopo = list_trees.size();

  // std::cout << "numUniqTopo = " << numUniqTopo <<"\n";

  numSameCoalEvents_inAncPop.resize(numUniqTopo);
  for(unsigned int i=0; i<numUniqTopo; i++)
    {
      unsigned int nGeneCopies = list_trees.at(i)->getSize();
      numSameCoalEvents_inAncPop.at(i).resize(nGeneCopies-2);
      for(unsigned int j=0; j<nGeneCopies-2; j++)
	{ 
	  numSameCoalEvents_inAncPop.at(i).at(j)=1;
	  unsigned int nKinds = nKinds_lineages.at(i).at(j);
	  unsigned int nKinds_next = nKinds_lineages.at(i).at(j+1);
	  std::vector<unsigned int> states_freqs = states_observed_freq.at(i).at(j);
	  std::vector<unsigned int> states_freqs_next = states_observed_freq.at(i).at(j+1);
	  std::vector<unsigned int> states = states_observed.at(i).at(j);
	  std::vector<unsigned int> states_next = states_observed.at(i).at(j+1);
	  unsigned int nLin_involvedInCoal = 0;
	  for(unsigned int k=0; k<nKinds; k++)
	    {
	      if(states_freqs.at(k) >=2 && nLin_involvedInCoal <2)
		{
		  unsigned int find_sameState=0;
		  unsigned int count=0;
		  while(count<nKinds_next && find_sameState==0)
		    {
		      if(states.at(k)==states_next.at(count))
			{
			  find_sameState=1;
			}
		      else
			count++;			  
		    }
		  if(find_sameState==1)
		    {
		      if(states_freqs.at(k) > states_freqs_next.at(count))
			{
			  unsigned int diff = states_freqs.at(k) - states_freqs_next.at(count);
			  if(diff==1)
			    {	  
			      numSameCoalEvents_inAncPop.at(i).at(j) *= states_freqs.at(k);
			      nLin_involvedInCoal += diff;
			    }
			  else if(diff==2)
			    {			  
			      numSameCoalEvents_inAncPop.at(i).at(j) *= states_freqs.at(k)*(states_freqs.at(k)-1)/2;
			      nLin_involvedInCoal += diff;
			    }
			  else
			    {
			      std::cout << "\n*** Error In Chain::compute_numSameCoalEvents_inAncPop ***\n";
			      std::cout << "states.at(k) = " << states.at(k) 
					<<" states_next.at(count) = " << states_next.at(count) <<"\n"
					<< " states_freqs.at(k) = " << states_freqs.at(k)
					<< " states_freqs_next.at(count) = " << states_freqs_next.at(count) <<"\n\n";
			    }
			}
		      else if(states_freqs.at(k) < states_freqs_next.at(count))
			{
			  std::cout << "\n*** Error In Chain::compute_numSameCoalEvents_inAncPop ***\n";
			  std::cout << "states.at(k) = " << states.at(k) 
				    <<" states_next.at(count) = " << states_next.at(count) <<"\n"
				    << " states_freqs.at(k) = " << states_freqs.at(k)
				    << " states_freqs_next.at(count) = " << states_freqs_next.at(count) <<"\n\n";
			}
		    }
		  else // no same state
		    {
		      if(states_freqs.at(k)==2)
			{	  
			  numSameCoalEvents_inAncPop.at(i).at(j) *= states_freqs.at(k)*(states_freqs.at(k)-1)/2;
			  nLin_involvedInCoal += 2;			  
			}
		      else
			{
			  std::cout << "\n*** Error In Chain::compute_numSameCoalEvents_inAncPop ***\n";
			  std::cout << "k= " << k << " states.at(k) = " << states.at(k) 
				    << " states_freqs.at(k) = " << states_freqs.at(k) <<"\n\n";
			  std::cout << "treeID = " << i << " j(the id for coalescent event)= " << j <<"\n";
			  for(unsigned int kk=0; kk<states.size(); kk++)
			    {
			      std::cout << "kk=" << kk <<"states.at(kk) = " << states.at(kk) 
					<< " states_freqs.at(kk) = " << states_freqs.at(kk) <<"\n"; 
			    }
			}
		    }
		} // End of if(states_freqs.at(k) >=2 && nLin_involvedInCoal <2)
	    }	  
	}      
    }
  
  return;
}

void Chain::compute_possiblePaths_nPossibleCoalEvents(unsigned int nPops)
{
  unsigned int numUniqTopo = list_trees.size();
  possiblePaths.resize(numUniqTopo);
  nPossibleCoalEvents.resize(numUniqTopo);
  for(unsigned int i=0; i<numUniqTopo; i++)
    {
      unsigned int nGeneCopies = list_trees.at(i)->getSize();
      possiblePaths.at(i).resize(nGeneCopies-1);
      nPossibleCoalEvents.at(i).resize(nGeneCopies-1);
      for(unsigned int j=0; j<nGeneCopies-1; j++)
	{
	  if(j==0)
	    possiblePaths.at(i).at(j).resize(nPops+1);
	  else
	    possiblePaths.at(i).at(j).resize(nPops*2);
	  
	  nPossibleCoalEvents.at(i).at(j).resize(nPops+1); // the last case is for the ancestral population.
	}
    }
  
  unsigned int id_initialState = find_initialState(nPops);
  for(unsigned int i=0; i<numUniqTopo; i++)
    {
      // REMOVE
      /*
      std::cout << "tree is ";
      list_trees.at(i)-> print_topo();
      std::cout << "id_initialState = " <<id_initialState << "\n";
      */

      unsigned int nGeneCopies = list_trees.at(i)->getSize();
      for(unsigned int j=0; j<nGeneCopies-1; j++)
	{
	  // std::cout << "j = " << j <<"\n\n";
	  
	  // std::cout << stateSpaces.at(i).at(j) <<"\n\n";
	  
	  
	  if(j==0) // initial state of a tree topology
	    {
	      possiblePaths.at(i).at(j).at(0).push_back(id_initialState);
	    }
	  
	  unsigned int nKinds = states_observed.at(i).at(j).size();
	  
	  // Final state of each time period
	  if(j!= nGeneCopies-2) // there are more than two lineages
	    {
	      unsigned int next_nKinds =states_observed.at(i).at(j+1).size();
	      
	      // Find the minimum frequency of each kind of lineages so that the next coalescent event can happen
	      Eigen::MatrixXd minimumFreq = find_minimumFreq_forNextTimePeriod(i,j);
	      
	      // std::string coalLin = find_nextCoalLineage(i,j);
	      // unsigned int coalLin = states_observed.at(i).at(j+1).at(next_nKinds-1);
	      
	      // REMOVE
	      //std::cout << "minimumFreq  = \n" << minimumFreq <<"\n";
	      
	      //Eigen::MatrixXd remainingFreq; // the number of lineages in each kind after remove the lineages to coalesce in the next time period


	      // Computing the number of possible coalescent events in the ancestral population.
	      unsigned int count_coalEvents = 1;
	      unsigned int coal_possible = 1;
	      unsigned int count = 0;
	      while(coal_possible == 1 && count <nKinds)
		{
		  unsigned int numLin = states_observed_freq.at(i).at(j).at(count);
		  unsigned int requiredNum = minimumFreq(count);
		  // std::cout << "popID = " << nPops << " numLin = " << numLin << " and requiredNum = " << requiredNum <<"\n";
		  if(numLin < requiredNum)
		    {
		      coal_possible = 0;
		    }
		  else if(requiredNum > 0)
		    {
		      count_coalEvents *= choose(numLin,requiredNum);
		    }
		  count++;
		}
	      if(coal_possible == 1)			    
		nPossibleCoalEvents.at(i).at(j).at(nPops).push_back(count_coalEvents);			  		  
	    
	      
	      
	      unsigned int spaceSize = stateSpaces.at(i).at(j).rows();
	      for(unsigned int id_state=0; id_state<spaceSize; id_state++)
		{
		  for(unsigned int p=0; p<nPops; p++)
		    {    
		      unsigned int count_coalEvents = 1;
		      unsigned int coal_possible = 1;
		      unsigned int count = 0;
		      while(coal_possible == 1 && count <nKinds)
			{
			  unsigned int numLin = stateSpaces.at(i).at(j)(id_state,p*nKinds+count);
			  unsigned int requiredNum = minimumFreq(count);
			  // std::cout << "popID = " << p << "numLin = " << numLin << " and requiredNum = " << requiredNum <<"\n";
			  if(numLin < requiredNum)
			    {
			      coal_possible = 0;
			    }
			  else if(requiredNum > 0)
			    {
			      count_coalEvents *= choose(numLin,requiredNum);
			    }
			  count++;
			}
		      
		      if(coal_possible == 1)
			{
			  nPossibleCoalEvents.at(i).at(j).at(p).push_back(count_coalEvents);
			  
			  if(j==0)
			    possiblePaths.at(i).at(j).at(1+p).push_back(id_state);
			  else
			    possiblePaths.at(i).at(j).at(nPops+p).push_back(id_state);
			  
			  // Find next initial state
			  Eigen::MatrixXd nextState(1,next_nKinds*nPops);
			  nextState.setZero();
			  unsigned int count_crr = 0;
			  for(unsigned int k=0; k< next_nKinds-1; k++)
			    {
			      unsigned int state_next = states_observed.at(i).at(j+1).at(k);
			      unsigned int found_sameState = 0;
			      while(found_sameState == 0 && count_crr < nKinds)
				{
				  unsigned int state_crr = states_observed.at(i).at(j).at(count_crr);
				  if(state_next == state_crr)
				    {
				      found_sameState = 1;
				      for(unsigned int pp=0; pp<nPops; pp++)
					{
					  if(pp == p)
					    nextState(0,pp*next_nKinds+k) =
					      stateSpaces.at(i).at(j)(id_state,pp*nKinds+ count_crr) - minimumFreq(0,count_crr);
					  else
					    nextState(0,pp*next_nKinds+k) = stateSpaces.at(i).at(j)(id_state,pp*nKinds+count_crr);
					}
				    }
				  count_crr++;
				}
				}
			      for(unsigned int pp=0; pp<nPops; pp++)
				{
				  if(pp == p)
				    nextState(0,pp*next_nKinds+next_nKinds-1) =1;
				}
			      
			      // REMOVE
			      //std::cout << "crr state is " <<	stateSpaces.at(i).at(j).row(id_state) <<"\n";
			      //std::cout << "next state is " << nextState <<"\n";
			      
			      unsigned int id_nextState = 0;
			      unsigned int found_id_nextState =0;
			      while(found_id_nextState == 0 && id_nextState <stateSpaces.at(i).at(j+1).rows())
				{
				  if(stateSpaces.at(i).at(j+1).row(id_nextState) == nextState)
				    {
				      found_id_nextState = 1;
				      possiblePaths.at(i).at(j+1).at(p).push_back(id_nextState);
				    }
				  id_nextState++;
				}
			    } // END of if(coal_possible == 1)
		    } // END of for(unsigned int p=0; p<nPops+1; p++)
		} // END of for(unsigned int id_state=0; id_state<spaceSize; id_state++)
		  
	      // REMOVE
	      /*
		std::cout << "minmumFreq = \n"<< minimumFreq <<"\n";
		std::cout << stateSpaces.at(i).at(j) <<"\n\n";
		for(unsigned int p=0; p<nPops; p++)
		{
		for(unsigned int k=0; k< possiblePaths.at(i).at(j).at(p).size(); k++)
		std::cout <<"p="<<p <<" k="<< k<< " id = "<< possiblePaths.at(i).at(j).at(p).at(k) <<"\n";
				}
	      */
	      
	      
	    }
	  else // there are two lineages
	    {
	      
	      for(unsigned int p=0; p<nPops+1; p++)
		nPossibleCoalEvents.at(i).at(j).at(p).push_back(1);
	      
	      for(unsigned int id_state =0; id_state < stateSpaces.at(i).at(j).rows(); id_state++)
		{
		  for(unsigned int p=0; p<nPops; p++)
		    {
		      //std::cout <<" stateSpaces.at(i).at(j).block(id_state,p*nKinds,1,nKinds) =" 
		      //	<< stateSpaces.at(i).at(j).block(id_state,p*nKinds,1,nKinds) <<"\n";
		      if(stateSpaces.at(i).at(j).block(id_state,p*nKinds,1,nKinds).sum() == 2)
			{
			  if(j==0)
			    possiblePaths.at(i).at(j).at(1+p).push_back(id_state);
			  else
			    possiblePaths.at(i).at(j).at(nPops+p).push_back(id_state);
			}
		    }
		}
	    }
	  
	}
    }
  
 
  
  return;
}



void Chain::compute_eigenValuesVectors_rateMat(popTree* poptree)
{
	unsigned int nPops = poptree->size();

	// FIXME YC 5/9/2014
	// It works for up to 2 population
	// Note: It might be better to have a migration rate matrix and a vector of population sizes
	// as members of class 'popTree'.
	Eigen::MatrixXd migRate(nPops,nPops);
	Eigen::MatrixXd popSize(1,nPops);
	migRate.setZero();
	popSize.setZero();

	for(unsigned int p1=0; p1< nPops; p1++)
	{
		for(unsigned int p2=0; p2< nPops; p2++)
		{
			if(p1 != p2)
				migRate(p1,p2) = poptree->find_migrationRate(p1+1, p2+1);
		}
		popSize(0,p1) = poptree->find_popSize(p1+1);
	}

	// REMOVE
	//std::cout << "Migration rate matrix is\n";
	//std::cout << migRate <<"\n";
	//std::cout << "Population sizes are \n";
	//std::cout << popSize <<"\n";

	unsigned int numUniqTopo = list_trees.size();
	Qeigen.resize(numUniqTopo);
	// transitionRateMat.resize(numUniqTopo);
	for(unsigned int i=0; i<numUniqTopo; i++)
	{
		unsigned int nGeneCopies = list_trees.at(i)->getSize();
		Qeigen.at(i).resize(nGeneCopies-1);
		// transitionRateMat.at(i).resize(nGeneCopies-1);
		for(unsigned int j=0; j<nGeneCopies-1; j++)
		{
			unsigned int n_states = stateSpaces.at(i).at(j).rows();
			unsigned int nKinds = nKinds_lineages.at(i).at(j); // states_observed.at(i).at(j).size(); // the number of kinds of lineages
			Qeigen.at(i).at(j).resize(3);
			Eigen::MatrixXd matQ(n_states+1,n_states+1);
			matQ.setZero();
			// transitionRateMat.at(i).at(j).resize(n_states+1,n_states+1);
			//transitionRateMat.at(i).at(j).setZero();

			// Transition from state r to c
			for(unsigned int r =0; r<n_states; r++)
			{
				for(unsigned int c = r+1; c<n_states; c++)
				{
					Eigen::MatrixXd diff = stateSpaces.at(i).at(j).row(r) - stateSpaces.at(i).at(j).row(c);
					if(diff.cwiseAbs().sum() == 2 && diff.sum()==0) // if the sum of absolute elements is 2 and the element sum is 0
					{
						// FIXME YC 5/9/2014
						// The following part works for up to 2 populations.
						unsigned int found = 0;
						unsigned int count = 0;
						unsigned int popID_from = 0;
						unsigned int popID_to = 0;
						unsigned int id_kind = 0;
						while(found==0 && count < nPops)
						{
							if(diff.block(0,nKinds*count,1,nKinds).sum() == 1)
							{
								popID_from = count+1;
							}
							else if(diff.block(0,nKinds*count,1,nKinds).sum()  == -1)
							{
								popID_to = count+1;
							}
							if(popID_from > 0 && popID_to >0)
							{
								found = 1;
								unsigned int k=0;
								unsigned int find_kind =0;
								while(k < nKinds && find_kind ==0)
								{
									int kind = diff.block(0,nKinds*count,1,nKinds)(0,k);
									if(kind == 1 || kind == -1)
									{
										id_kind = k;
										find_kind = 1;
									}
									k++;
								}
							}
							count++;
						}
						popID_from--;
						popID_to--;

						// REMOVE
						//std::cout << "r = " << r <<" c = "<<c <<"\n";
						//std::cout << "From: "<<stateSpaces.at(i).at(j).row(r) <<"\n";
						//std::cout << "To: "<< stateSpaces.at(i).at(j).row(c) <<"\n";
						//std::cout << "popID_from: " << popID_from <<", popID_to: " << popID_to << "\n";
						//std::cout << "Kind of lineage is " << id_kind <<"\n";

						// Transition rate from state r to state c
						matQ(r,c) = migRate(popID_from, popID_to) *
								stateSpaces.at(i).at(j)(r, nKinds*(popID_from)+id_kind);

						// Transition rate from state c to state r
						matQ(c,r) = migRate(popID_to, popID_from) *
								stateSpaces.at(i).at(j)(c, nKinds*(popID_to)+id_kind);

					}
				}

				// absorbing states
				for(unsigned int p=0; p < nPops; p++)
				{
					unsigned int totalL = stateSpaces.at(i).at(j).block(r,p*nKinds,1,nKinds).sum();
					if(totalL >= 2)
						matQ(r,n_states) += totalL*(totalL-1)/popSize(0,p);
				}

				// diagonal elements
				matQ(r,r) = -matQ.row(r).sum();
			}

			// REMOVE
			//std::cout << "State space is\n";
			//std::cout << stateSpaces.at(i).at(j)<<"\n\n";
			//std::cout << "Transition rate Matrix is \n";
			//std::cout << matQ <<"\n\n";

			// Matrix decomposition
			Eigen::EigenSolver<Eigen::MatrixXd> svd;
			svd.compute(matQ);
			// eigen value and eigen vector matrices
			Qeigen.at(i).at(j).at(0) = svd.eigenvalues().asDiagonal();
			Qeigen.at(i).at(j).at(1) = svd.eigenvectors();
			Qeigen.at(i).at(j).at(2) = svd.eigenvectors().inverse();
			// Checking eigen values and eigen vectors
			Eigen::MatrixXcd I = Qeigen.at(i).at(j).at(1) * Qeigen.at(i).at(j).at(2);
			if(abs(I.sum().real() - n_states-1) > 3.0e-7)
			  {
			    poptree->print_allPopTree();
			    std::cout << "Eigen vectors may not be accurate. That is, V*V^(-1) != I.\n"
				      << "The sum of all elements in V*V^(-1) is " << I.sum().real()
				      << "\nThe sum of absolute differences is "<< abs(I.sum().real() - n_states-1) <<".\n";
			  }
			Eigen::MatrixXcd rateMat2 = Qeigen.at(i).at(j).at(1) * svd.eigenvalues().asDiagonal() * Qeigen.at(i).at(j).at(2);
			if(abs(rateMat2.sum().real() - matQ.sum()) > 3.0e-7)
			  {
			    poptree->print_allPopTree();
			    std::cout << "Eigen values and eigen vectors may not be accurate. rateMat != V*D*V^(-1).\n"
				      <<"The absolute difference of the sums of elements in rateMat and V*D*V^(-1) is "
				      << abs(rateMat2.sum().real() - matQ.sum()) <<"\n";
			  }

		}
	}

	return;

}





void Chain::compute_coefficientsOfElementsInTransitionRateMat(unsigned int nPops)
{	
  //std::cout << "Chain::compute_coefficientsOfElementsInTransitionRateMat()\n";
  //std::cout << "nPops = "<< nPops <<"\n";

  // FIXME YC 5/9/2014
  // It works for up to 2 population
  // Note: It might be better to have a migration rate matrix and a vector of population sizes
  // as members of class 'popTree'.
  /*
  Eigen::MatrixXd migRate(nPops,nPops);
  Eigen::MatrixXd popSize(1,nPops);
  migRate.setZero();
  popSize.setZero();


  for(unsigned int p1=0; p1< nPops; p1++)
    {
      for(unsigned int p2=0; p2< nPops; p2++)
	{
	  if(p1 != p2)
	    migRate(p1,p2) = poptree->find_migrationRate(p1+1, p2+1);
	}
      popSize(0,p1) = poptree->find_popSize(p1+1);
    }
  */

  unsigned int numUniqTopo = stateSpaces.size(); // list_trees.size();
  
  // Qeigen.resize(numUniqTopo);
  coeff4TransitionRateMat.resize(numUniqTopo);
  // transitionRateMat.resize(numUniqTopo);
  for(unsigned int i=0; i<numUniqTopo; i++)
    {
      unsigned int nGeneCopies = stateSpaces.at(i).size()+1; //list_trees.at(i)->getSize();
      // Qeigen.at(i).resize(nGeneCopies-1);
      coeff4TransitionRateMat.at(i).resize(nGeneCopies-1);
      for(unsigned int j=0; j<nGeneCopies-1; j++)
	{
	  unsigned int n_states = stateSpaces.at(i).at(j).rows();
	  unsigned int nKinds = nKinds_lineages.at(i).at(j); // states_observed.at(i).at(j).size(); // the number of kinds of lineages
	  
	  coeff4TransitionRateMat.at(i).at(j).resize(4);
	  
	  // Qeigen.at(i).at(j).resize(3);
	  /*
	  Eigen::MatrixXd mat_pop1(n_states+1,n_states+1);
	  mat_pop1.setZero();
	  Eigen::MatrixXd mat_pop2(n_states+1,n_states+1);
	  mat_pop2.setZero();
	  Eigen::MatrixXd mat_mig1(n_states+1,n_states+1);
	  mat_mig1.setZero();
	  Eigen::MatrixXd mat_mig2(n_states+1,n_states+1);
	  mat_mig22.setZero();
	  */
	  for(unsigned int r=0; r< 4; r++)
	    {
	      coeff4TransitionRateMat.at(i).at(j).at(r).resize(n_states+1,n_states+1);
	      coeff4TransitionRateMat.at(i).at(j).at(r).setZero();
	    }
	  
	  // Transition from state r to c
	  for(unsigned int r =0; r<n_states; r++)
	    {
	      for(unsigned int c = r+1; c<n_states; c++)
		{
		  Eigen::MatrixXd diff = stateSpaces.at(i).at(j).row(r) - stateSpaces.at(i).at(j).row(c);
		  if(diff.cwiseAbs().sum() == 2 && diff.sum()==0) // if the sum of absolute elements is 2 and the element sum is 0
		    {
		      // FIXME YC 5/9/2014
		      // The following part works for up to 2 populations.
		      unsigned int found = 0;
		      unsigned int count = 0;
		      unsigned int popID_from = 0;
		      unsigned int popID_to = 0;
		      unsigned int id_kind = 0;
		      while(found==0 && count < nPops)
			{
			  if(diff.block(0,nKinds*count,1,nKinds).sum() == 1)
			    {
			      popID_from = count+1;
			    }
			  else if(diff.block(0,nKinds*count,1,nKinds).sum()  == -1)
			    {
			      popID_to = count+1;
			    }
			  if(popID_from > 0 && popID_to >0)
			    {
			      found = 1;
			      unsigned int k=0;
			      unsigned int find_kind =0;
			      while(k < nKinds && find_kind ==0)
				{
				  int kind = diff.block(0,nKinds*count,1,nKinds)(0,k);
				  if(kind == 1 || kind == -1)
				    {
				      id_kind = k;
				      find_kind = 1;
				    }
				  k++;
				}
			    }
			  count++;
			}
		      popID_from--;
		      popID_to--;
		      
		      
		      // Transition rate from state r to state c	      
		      // matQ(r,c) = migRate(popID_from, popID_to) *stateSpaces.at(i).at(j)(r, nKinds*(popID_from)+id_kind);
		      if(popID_from < popID_to) // migration rate 1
			coeff4TransitionRateMat.at(i).at(j).at(nPops)(r,c) = stateSpaces.at(i).at(j)(r, nKinds*(popID_from)+id_kind);
		      else // migration rate 2
			coeff4TransitionRateMat.at(i).at(j).at(nPops+1)(r,c) = stateSpaces.at(i).at(j)(r, nKinds*(popID_from)+id_kind);
		      
		      
		      // Transition rate from state c to state r
		      // matQ(c,r) = migRate(popID_to, popID_from) *stateSpaces.at(i).at(j)(c, nKinds*(popID_to)+id_kind);
		      if(popID_to < popID_from) // migration rate 1
			coeff4TransitionRateMat.at(i).at(j).at(nPops)(c,r) = stateSpaces.at(i).at(j)(c, nKinds*(popID_to)+id_kind);
		      else // migration rate 2
			coeff4TransitionRateMat.at(i).at(j).at(nPops+1)(c,r) = stateSpaces.at(i).at(j)(c, nKinds*(popID_to)+id_kind);
		      
		    }
		}
	      
	      // absorbing states
	      for(unsigned int p=0; p < nPops; p++)
		{
		  unsigned int totalL = stateSpaces.at(i).at(j).block(r,p*nKinds,1,nKinds).sum();
		  if(totalL >= 2)
		    {
		      // matQ(r,n_states) += totalL*(totalL-1)/popSize(0,p);
		      coeff4TransitionRateMat.at(i).at(j).at(p)(r,n_states) = totalL*(totalL-1);
		    }
		}
	      
	      // diagonal elements
	      // matQ(r,r) = -matQ.row(r).sum();
	    }
	  
	 
	  
	}
    }
  
  return;

}



// YC 3/25/2015
void Chain::compute_eigenValuesVectors_subMatOfEigenVectors(popTree* poptree)
{
  std::chrono::high_resolution_clock::time_point start_t, end_t;
  start_t= std::chrono::high_resolution_clock::now();

  unsigned int nPops = poptree->size();

  // FIXME YC 5/9/2014
  // It works for up to 2 population
  // Note: It might be better to have a migration rate matrix and a vector of population sizes
  // as members of class 'popTree'.
  Eigen::MatrixXd migRate(nPops,nPops);
  Eigen::MatrixXd popSize(1,nPops);
  migRate.setZero();
  popSize.setZero();
  
  for(unsigned int p1=0; p1< nPops; p1++)
    {
      for(unsigned int p2=0; p2< nPops; p2++)
	{
	  if(p1 != p2)
	    migRate(p1,p2) = poptree->find_migrationRate(p1+1, p2+1);
	}
      popSize(0,p1) = poptree->find_popSize(p1+1);
    }

  unsigned int numUniqTopo = stateSpaces.size(); // list_trees.size();
  
  Qeigen.resize(numUniqTopo);
  for(unsigned int i=0; i<numUniqTopo; i++)
    {
      unsigned int nGeneCopies = stateSpaces.at(i).size()+1; //list_trees.at(i)->getSize();
      Qeigen.at(i).resize(nGeneCopies-1);
      // transitionRateMat.at(i).resize(nGeneCopies-1);
      for(unsigned int j=0; j<nGeneCopies-1; j++)
	{
	  Qeigen.at(i).at(j).resize(3);
	  

	  unsigned int n_states = stateSpaces.at(i).at(j).rows();


	  // unsigned int nKinds = nKinds_lineages.at(i).at(j); // states_observed.at(i).at(j).size(); // the number of kinds of lineages
	  Qeigen.at(i).at(j).resize(3);
	  Eigen::MatrixXd matQ(n_states+1,n_states+1);
	  matQ.setZero();

	  for(unsigned int p=0; p<nPops; p++)
	    matQ = matQ + coeff4TransitionRateMat.at(i).at(j).at(p)/popSize(0,p);

	  unsigned int count = nPops;
	  for(unsigned int p1=0; p1< nPops; p1++)
	    {
	      for(unsigned int p2=0; p2< nPops; p2++)
		{
		  if(p1 != p2)
		    {
		      matQ = matQ + coeff4TransitionRateMat.at(i).at(j).at(count) * migRate(p1,p2);
		      count++;
		    }
		}
	    }

	  // diagonal elements
	  for(unsigned int r =0; r<n_states; r++)
	    matQ(r,r) = -matQ.row(r).sum();

	  // Matrix decomposition
	  Eigen::EigenSolver<Eigen::MatrixXd> svd;
	  svd.compute(matQ);
	  // eigen value and eigen vector matrices
	  Qeigen.at(i).at(j).at(0) = svd.eigenvalues();
	  //Qeigen.at(i).at(j).at(0) = svd.eigenvalues().asDiagonal();
	  Qeigen.at(i).at(j).at(1) = svd.eigenvectors();
	  Qeigen.at(i).at(j).at(2) = svd.eigenvectors().inverse();
	  // Checking eigen values and eigen vectors
	  /*
	  Eigen::MatrixXcd I = Qeigen.at(i).at(j).at(1) * Qeigen.at(i).at(j).at(2);
	  if(abs(I.sum().real() - n_states-1) > pow(10,-4))
	    {
	      poptree->print_allPopTree();
	      std::cout << "*** Warning *** In Chain::compute_eigenValuesVectors()"
			<<" Eigen vectors may not be accurate. That is, V*V^(-1) != I.\n"
			<< "The sum of all elements in V*V^(-1) is " << I.sum().real()
			<< "\nThe sum of absolute differences is "<< abs(I.sum().real() - n_states-1) <<".\n";
	      std::cout << "matQ = " << matQ <<"\n";
	    }
	  Eigen::MatrixXcd rateMat2 = Qeigen.at(i).at(j).at(1) * svd.eigenvalues().asDiagonal() * Qeigen.at(i).at(j).at(2);
	  if(abs(rateMat2.sum().real() - matQ.sum()) > pow(10,-5))
	    {
	      poptree->print_allPopTree();
	      std::cout << "*** Warning *** In Chain::compute_eigenValuesVectors()"
			<< "Eigen values and eigen vectors may not be accurate. rateMat != V*D*V^(-1).\n"
			<<"The absolute difference of the sums of elements in rateMat and V*D*V^(-1) is "
			<< abs(rateMat2.sum().real() - matQ.sum()) <<"\n";
	      std::cout << "matQ = " << matQ <<"\n";
	      std::cout << "rateMat2 = " << rateMat2 <<"\n";
	    }
	  */
	   

	}
    }

  //--- END of computing eigen values and eigen vectors ---//
  end_t= std::chrono::high_resolution_clock::now();
  eachComputingTime_eigen = std::chrono::duration_cast<std::chrono::microseconds>(end_t - start_t);


  //-- Getting sub-matrices of eigen vector matrices ---//
  start_t= std::chrono::high_resolution_clock::now();

  subMatV.resize(numUniqTopo);
  eigenvalues.resize(numUniqTopo);
  for(unsigned int i=0; i<numUniqTopo; i++)
    {
      unsigned int nGeneCopies = stateSpaces.at(i).size()+1; 
      subMatV.at(i).resize(nGeneCopies-1);
      eigenvalues.at(i).resize(nGeneCopies-1);
      for(unsigned int j=0; j<nGeneCopies-1; j++)
	{
	  eigenvalues.at(i).at(j) = Qeigen.at(i).at(j).at(0);
	  subMatV.at(i).at(j).resize(3);
	  // Each 'subMatV.at(i).at(j)' has 4 matrices:
	  //  case 0: submatrix of eigen vector matrix for a coalescent event (possible paths of a coalescent event)
	  //  case 1: submatrix of inverse eigenvector matrix for a coalescent event (possible paths of a coalescent event)
	  //  case 2: submstrix of inverse eigenvector matrix for population merging (all transient paths)
	  

	  // case 0 
	  if(j==0)
	    { 
	      subMatV.at(i).at(j).at(0) =  Qeigen.at(i).at(j).at(1).row(possiblePaths.at(i).at(j).at(0).at(0));	      
	    }
	  else
	    {
	      unsigned int nr = 0;
	      for(unsigned int p=0; p<nPops; p++)
		{
		  nr += possiblePaths.at(i).at(j).at(p).size();
		}
	      subMatV.at(i).at(j).at(0).resize(nr, Qeigen.at(i).at(j).at(1).cols());
	      unsigned int count_r = 0;
	      for(unsigned int p=0; p<nPops; p++)
		{
		  for(unsigned int k=0; k<possiblePaths.at(i).at(j).at(p).size(); k++)
		    {
		      unsigned int id_row = possiblePaths.at(i).at(j).at(p).at(k);
		      subMatV.at(i).at(j).at(0).row(count_r) =Qeigen.at(i).at(j).at(1).row(id_row);
		      count_r++;
		    }
		}    
	    }

	  // case 1
	  unsigned int id_next = 0;
	  if(j==0)
	    id_next =1;
	  else
	    id_next = nPops;
	  unsigned int nc = 0;
	  for(unsigned int p=0; p<nPops; p++)
	    nc += possiblePaths.at(i).at(j).at(p+id_next).size();
	  
	  subMatV.at(i).at(j).at(1).resize(Qeigen.at(i).at(j).at(2).rows(), nc);
	  unsigned int count = 0;
	  for(unsigned int p=0; p<nPops; p++)
	    {
	      for(unsigned int k=0; k<possiblePaths.at(i).at(j).at(p+id_next).size(); k++)
		{
		  unsigned int id_col = possiblePaths.at(i).at(j).at(p+id_next).at(k);
		  unsigned int nCoalEvents = nPossibleCoalEvents.at(i).at(j).at(p).at(k);
		  subMatV.at(i).at(j).at(1).col(count) =Qeigen.at(i).at(j).at(2).col(id_col)* nCoalEvents *2/popSize(0,p);
		  count++;
		}
	    }	 

	  // case 2
	  unsigned int nTransientStates = Qeigen.at(i).at(j).at(2).cols()-1;
	  subMatV.at(i).at(j).at(2) = Qeigen.at(i).at(j).at(2).leftCols(nTransientStates); // keep the first 'nTransientStates' columns
	}
    }
  

  //-- End of getting sub-matrices of eigen vector matrices ---//
  end_t= std::chrono::high_resolution_clock::now();
  eachComputingTime_eigen_subMatOfEigenVectors = std::chrono::duration_cast<std::chrono::microseconds>(end_t - start_t);

  Qeigen.resize(0);

  return;
  
}



// YC 2/24/2016
void Chain::compute_eigenValuesVectors_subMatOfEigenVectors_MPI(popTree* poptree)
{
  std::chrono::high_resolution_clock::time_point start_t, end_t;
  start_t= std::chrono::high_resolution_clock::now();

  unsigned int nPops = poptree->size();

  // Assigning topologies to CPUs
  unsigned int numUniqTopo = stateSpaces.size();
  unsigned int subNumTopo = static_cast<unsigned int>(numUniqTopo/nProcesses);
  unsigned int smallerSubNumTopo = subNumTopo;
  int nRemainder = numUniqTopo - subNumTopo*nProcesses;   
  for(unsigned int i=0; i<nProcesses; i++)
    {
      if(cpuID == i && nProcesses-i <= nRemainder)
	subNumTopo++;
    }   
  if(nProcesses - cpuID > nRemainder)
    topoID_start = subNumTopo*cpuID;
  else
    topoID_start = (subNumTopo-1)*(nProcesses-nRemainder)
      + subNumTopo*(cpuID-(nProcesses-nRemainder));
      
  topoID_end = topoID_start + subNumTopo-1;

  /*  
  std::cout << "cpuID = " << cpuID << " topoID_start = " << topoID_start <<" topoID_end = " << topoID_end 
	    << " subNumTopo = "  << subNumTopo
	    <<"\n";
  */
  

  // FIXME YC 5/9/2014
  // It works for up to 2 population
  // Note: It might be better to have a migration rate matrix and a vector of population sizes
  // as members of class 'popTree'.
  Eigen::MatrixXd migRate(nPops,nPops);
  Eigen::MatrixXd popSize(1,nPops);
  migRate.setZero();
  popSize.setZero();
  
  for(unsigned int p1=0; p1< nPops; p1++)
    {
      for(unsigned int p2=0; p2< nPops; p2++)
	{
	  if(p1 != p2)
	    migRate(p1,p2) = poptree->find_migrationRate(p1+1, p2+1);
	}
      popSize(0,p1) = poptree->find_popSize(p1+1);
    }

  
  Qeigen.resize(subNumTopo);
  for(int i=0; i<numUniqTopo; i++)
    {
      unsigned int nGeneCopies = stateSpaces.at(i).size()+1;
      if(subNumTopo>0 && i >= topoID_start && i <= topoID_end)
	{
	  Qeigen.at(i-topoID_start).resize(nGeneCopies-1);
	  for(unsigned int j=0; j<nGeneCopies-1; j++)
	    {
	      Qeigen.at(i-topoID_start).at(j).resize(3);
	  

	      unsigned int n_states = stateSpaces.at(i).at(j).rows();

	      // Qeigen.at(i-topoID_start).at(j).resize(3);
	      Eigen::MatrixXd matQ(n_states+1,n_states+1);
	      matQ.setZero();

	      for(unsigned int p=0; p<nPops; p++)
		matQ = matQ + coeff4TransitionRateMat.at(i).at(j).at(p)/popSize(0,p);

	      unsigned int count = nPops;
	      for(unsigned int p1=0; p1< nPops; p1++)
		{
		  for(unsigned int p2=0; p2< nPops; p2++)
		    {
		      if(p1 != p2)
			{
			  matQ = matQ + coeff4TransitionRateMat.at(i).at(j).at(count) * migRate(p1,p2);
			  count++;
			}
		    }
		}
	      
	      // diagonal elements
	      for(unsigned int r =0; r<n_states; r++)
		matQ(r,r) = -matQ.row(r).sum();
	      
	      // Matrix decomposition
	      Eigen::EigenSolver<Eigen::MatrixXd> svd;
	      svd.compute(matQ);
	      // eigen value and eigen vector matrices
	      Qeigen.at(i-topoID_start).at(j).at(0) = svd.eigenvalues();
	      Qeigen.at(i-topoID_start).at(j).at(1) = svd.eigenvectors();
	      Qeigen.at(i-topoID_start).at(j).at(2) = svd.eigenvectors().inverse();
	      // Checking eigen values and eigen vectors
	      /*
		Eigen::MatrixXcd I = Qeigen.at(i).at(j).at(1) * Qeigen.at(i).at(j).at(2);
		if(abs(I.sum().real() - n_states-1) > pow(10,-4))
		{
		poptree->print_allPopTree();
		std::cout << "*** Warning *** In Chain::compute_eigenValuesVectors()"
		<<" Eigen vectors may not be accurate. That is, V*V^(-1) != I.\n"
		<< "The sum of all elements in V*V^(-1) is " << I.sum().real()
		<< "\nThe sum of absolute differences is "<< abs(I.sum().real() - n_states-1) <<".\n";
		std::cout << "matQ = " << matQ <<"\n";
		}
		Eigen::MatrixXcd rateMat2 = Qeigen.at(i).at(j).at(1) * svd.eigenvalues().asDiagonal() * Qeigen.at(i).at(j).at(2);
		if(abs(rateMat2.sum().real() - matQ.sum()) > pow(10,-5))
		{
		poptree->print_allPopTree();
		std::cout << "*** Warning *** In Chain::compute_eigenValuesVectors()"
		<< "Eigen values and eigen vectors may not be accurate. rateMat != V*D*V^(-1).\n"
		<<"The absolute difference of the sums of elements in rateMat and V*D*V^(-1) is "
		<< abs(rateMat2.sum().real() - matQ.sum()) <<"\n";
		std::cout << "matQ = " << matQ <<"\n";
		std::cout << "rateMat2 = " << rateMat2 <<"\n";
		}
	      */    
	    }
	}
  
      //--- END of computing eigen values and eigen vectors ---//
      end_t= std::chrono::high_resolution_clock::now();
      eachComputingTime_eigen = std::chrono::duration_cast<std::chrono::microseconds>(end_t - start_t);
    }
  //-- Getting sub-matrices of eigen vector matrices ---//
  start_t= std::chrono::high_resolution_clock::now();

  // std::cout << "cpuID = " << cpuID <<" here\n";

  subMatV.resize(numUniqTopo);
  eigenvalues.resize(numUniqTopo);
  for(int i=0; i<numUniqTopo; i++)
    {
      unsigned int nGeneCopies = stateSpaces.at(i).size()+1; 
      subMatV.at(i).resize(nGeneCopies-1);
      eigenvalues.at(i).resize(nGeneCopies-1);
      // std::cout << "cpuID = " << cpuID <<" here3\n";
      for(unsigned int j=0; j<nGeneCopies-1; j++)
	{
	  subMatV.at(i).at(j).resize(3);
	  if(subNumTopo>0 && i >= topoID_start && i <= topoID_end)
	    {
	      eigenvalues.at(i).at(j) = Qeigen.at(i-topoID_start).at(j).at(0);
	      // Each 'subMatV.at(i).at(j)' has 4 matrices:
	      //  case 0: submatrix of eigen vector matrix for a coalescent event (possible paths of a coalescent event)
	      //  case 1: submatrix of inverse eigenvector matrix for a coalescent event (possible paths of a coalescent event)
	      //  case 2: submstrix of inverse eigenvector matrix for population merging (all transient paths)
	      
	      
	      // case 0 
	      if(j==0)
		{ 
		  subMatV.at(i).at(j).at(0) =  Qeigen.at(i-topoID_start).at(j).at(1).row(possiblePaths.at(i).at(j).at(0).at(0));	      
		}
	      else
		{
		  unsigned int nr = 0;
		  for(unsigned int p=0; p<nPops; p++)
		    {
		      nr += possiblePaths.at(i).at(j).at(p).size();
		    }
		  subMatV.at(i).at(j).at(0).resize(nr, Qeigen.at(i-topoID_start).at(j).at(1).cols());
		  unsigned int count_r = 0;
		  for(unsigned int p=0; p<nPops; p++)
		    {
		      for(unsigned int k=0; k<possiblePaths.at(i).at(j).at(p).size(); k++)
			{
			  unsigned int id_row = possiblePaths.at(i).at(j).at(p).at(k);
			  subMatV.at(i).at(j).at(0).row(count_r) =Qeigen.at(i-topoID_start).at(j).at(1).row(id_row);
			  count_r++;
			}
		    }    
		}
	      
	      // case 1
	      unsigned int id_next = 0;
	      if(j==0)
		id_next =1;
	      else
		id_next = nPops;
	      unsigned int nc = 0;
	      for(unsigned int p=0; p<nPops; p++)
		nc += possiblePaths.at(i).at(j).at(p+id_next).size();
	      
	      subMatV.at(i).at(j).at(1).resize(Qeigen.at(i-topoID_start).at(j).at(2).rows(), nc);
	      unsigned int count = 0;
	      for(unsigned int p=0; p<nPops; p++)
		{
		  for(unsigned int k=0; k<possiblePaths.at(i).at(j).at(p+id_next).size(); k++)
		    {
		      unsigned int id_col = possiblePaths.at(i).at(j).at(p+id_next).at(k);
		      unsigned int nCoalEvents = nPossibleCoalEvents.at(i).at(j).at(p).at(k);
		      subMatV.at(i).at(j).at(1).col(count) =Qeigen.at(i-topoID_start).at(j).at(2).col(id_col)* nCoalEvents *2/popSize(0,p);
		      count++;
		    }
		}	 
	      
	      // case 2
	      unsigned int nTransientStates = Qeigen.at(i-topoID_start).at(j).at(2).cols()-1;
	      subMatV.at(i).at(j).at(2) = Qeigen.at(i-topoID_start).at(j).at(2).leftCols(nTransientStates); // keep the first 'nTransientStates' columns
	    }

	} // END of if(subNumTopo>0 && i >= topoID_start && i <= topoID_end)

      // std::cout << "cpuID = " << cpuID <<" i=" << i  << " here2\n";
      
      int senderCPUID = nProcesses-nRemainder+i;
      // std::cout << "cpuID = " << cpuID <<" i=" << i  <<" senderCPUID="<<senderCPUID<<"\n";
      if(smallerSubNumTopo !=0)
	{
	  senderCPUID = static_cast<unsigned int>(i/smallerSubNumTopo);   
	  if(nProcesses-nRemainder <=  senderCPUID)
	    senderCPUID = nProcesses-nRemainder + static_cast<unsigned int>((i-(nProcesses-nRemainder)*smallerSubNumTopo)/(smallerSubNumTopo+1));
	}
      /*
      std::cout << "cpuID = " << cpuID <<" i=" << i  
		<<" nProcesses="<<nProcesses <<" nRemainder"<<nRemainder
		<<" senderCPUID = "<< senderCPUID
		<< " \n";
      */
     
      for(unsigned int j=0; j< nGeneCopies-1; j++)
	{	  
	  unsigned int size = 0;
	  // eigen values
	  if(subNumTopo>0 && i >= topoID_start && i <= topoID_end)
	    {
	      size = eigenvalues.at(i).at(j).size();
	      // std::cout << "Sender's ID = " << cpuID << " size="<<size <<"\n";
	    }
	  MPI::COMM_WORLD.Bcast(&size, 1, MPI::UNSIGNED, senderCPUID);
	  MPI::COMM_WORLD.Barrier();	  
	  eigenvalues.at(i).at(j).resize(size);
	  for(unsigned int rr=0; rr<size; rr++)
	    {
	      std::complex<double> evalue = 0;
	      if(subNumTopo>0 && i >= topoID_start && i <= topoID_end)
		{
		  evalue = eigenvalues.at(i).at(j)(rr);
		}
	      MPI::COMM_WORLD.Bcast(&evalue, 1, MPI::DOUBLE_COMPLEX, senderCPUID);
	      MPI::COMM_WORLD.Barrier();
	      eigenvalues.at(i).at(j)(rr) = evalue;
	    }
 	  
	  unsigned int sizeR =0;
	  unsigned int sizeC =0;

	  // eigen vectors
	  if(subNumTopo>0 && i >= topoID_start && i <= topoID_end)
	    {
	      sizeR =subMatV.at(i).at(j).at(0).rows(); 
	      sizeC =subMatV.at(i).at(j).at(0).cols();
	    }	
	  MPI::COMM_WORLD.Bcast(&sizeR, 1, MPI::UNSIGNED, senderCPUID);
	  MPI::COMM_WORLD.Barrier();	  
	  MPI::COMM_WORLD.Bcast(&sizeC, 1, MPI::UNSIGNED, senderCPUID);
	  MPI::COMM_WORLD.Barrier();	  
	  subMatV.at(i).at(j).at(0).resize(sizeR, sizeC);
	  
	  for(unsigned int rr=0; rr<sizeR; rr++)
	    {
	      for(unsigned int cc=0; cc<sizeC; cc++)
		{
            	  std::complex<double> evalue = 0;
		  if(subNumTopo>0 && i >= topoID_start && i <= topoID_end)
		    {
	              evalue = subMatV.at(i).at(j).at(0)(rr,cc);
		    }
		  MPI::COMM_WORLD.Bcast(&evalue, 1, MPI::DOUBLE_COMPLEX, senderCPUID);
		  MPI::COMM_WORLD.Barrier();
		  subMatV.at(i).at(j).at(0)(rr,cc) = evalue;		  
		}	      
	    }


	  // inverse eigenvectors for a coalescent
	  sizeR=0; sizeC=0;
	  if(subNumTopo>0 && i >= topoID_start && i <= topoID_end)
	    {
	      sizeR =subMatV.at(i).at(j).at(1).rows(); 
	      sizeC =subMatV.at(i).at(j).at(1).cols();   
	    }	
	  MPI::COMM_WORLD.Bcast(&sizeR, 1, MPI::UNSIGNED, senderCPUID);
	  MPI::COMM_WORLD.Barrier();	  
	  MPI::COMM_WORLD.Bcast(&sizeC, 1, MPI::UNSIGNED, senderCPUID);
	  MPI::COMM_WORLD.Barrier();	  
	  subMatV.at(i).at(j).at(1).resize(sizeR, sizeC);
	  
	  for(unsigned int rr=0; rr<sizeR; rr++)
	    {
	      for(unsigned int cc=0; cc<sizeC; cc++)
		{
            	  std::complex<double> evalue = 0;
		  if(subNumTopo>0 && i >= topoID_start && i <= topoID_end)
		    {
	              evalue = subMatV.at(i).at(j).at(1)(rr,cc);
		    }
		  MPI::COMM_WORLD.Bcast(&evalue, 1, MPI::DOUBLE_COMPLEX, senderCPUID);
		  MPI::COMM_WORLD.Barrier();
		  subMatV.at(i).at(j).at(1)(rr,cc) = evalue;		  
		}	      
	    }

	  // inverse eigenvectors for population merging
	  sizeR=0; sizeC=0;
	  if(subNumTopo>0 && i >= topoID_start && i <= topoID_end)
	    {
	      sizeR =subMatV.at(i).at(j).at(2).rows(); 
	      sizeC =subMatV.at(i).at(j).at(2).cols();   
	    }	
	  MPI::COMM_WORLD.Bcast(&sizeR, 1, MPI::UNSIGNED, senderCPUID);
	  MPI::COMM_WORLD.Barrier();	  
	  MPI::COMM_WORLD.Bcast(&sizeC, 1, MPI::UNSIGNED, senderCPUID);
	  MPI::COMM_WORLD.Barrier();	  
	  subMatV.at(i).at(j).at(2).resize(sizeR, sizeC);
	  
	  for(unsigned int rr=0; rr<sizeR; rr++)
	    {
	      for(unsigned int cc=0; cc<sizeC; cc++)
		{
            	  std::complex<double> evalue = 0;
		  if(subNumTopo>0 && i >= topoID_start && i <= topoID_end)
		    {
	              evalue = subMatV.at(i).at(j).at(2)(rr,cc);
		    }
		  MPI::COMM_WORLD.Bcast(&evalue, 1, MPI::DOUBLE_COMPLEX, senderCPUID);
		  MPI::COMM_WORLD.Barrier();
		  subMatV.at(i).at(j).at(2)(rr,cc) = evalue;		  
		}	      
	    }

	}
    }

  //-- End of getting sub-matrices of eigen vector matrices ---//
  end_t= std::chrono::high_resolution_clock::now();
  eachComputingTime_eigen_subMatOfEigenVectors = std::chrono::duration_cast<std::chrono::microseconds>(end_t - start_t);

  Qeigen.resize(0);

  return;
  
}




// old version
// YC 2/16/2014
void Chain::compute_eigenValuesVectors(popTree* poptree)
{

  //clock_start = clock();

  unsigned int nPops = poptree->size();

  // FIXME YC 5/9/2014
  // It works for up to 2 population
  // Note: It might be better to have a migration rate matrix and a vector of population sizes
  // as members of class 'popTree'.
  Eigen::MatrixXd migRate(nPops,nPops);
  Eigen::MatrixXd popSize(1,nPops);
  migRate.setZero();
  popSize.setZero();
  
  for(unsigned int p1=0; p1< nPops; p1++)
    {
      for(unsigned int p2=0; p2< nPops; p2++)
	{
	  if(p1 != p2)
	    migRate(p1,p2) = poptree->find_migrationRate(p1+1, p2+1);
	}
      popSize(0,p1) = poptree->find_popSize(p1+1);
    }

  unsigned int numUniqTopo = stateSpaces.size(); // list_trees.size();
  
  Qeigen.resize(numUniqTopo);
  // transitionRateMat.resize(numUniqTopo);
  for(unsigned int i=0; i<numUniqTopo; i++)
    {
      unsigned int nGeneCopies = stateSpaces.at(i).size()+1; //list_trees.at(i)->getSize();
      Qeigen.at(i).resize(nGeneCopies-1);
      // transitionRateMat.at(i).resize(nGeneCopies-1);
      for(unsigned int j=0; j<nGeneCopies-1; j++)
	{
	  Qeigen.at(i).at(j).resize(3);
	  

	  unsigned int n_states = stateSpaces.at(i).at(j).rows();


	  // unsigned int nKinds = nKinds_lineages.at(i).at(j); // states_observed.at(i).at(j).size(); // the number of kinds of lineages
	  Qeigen.at(i).at(j).resize(3);
	  Eigen::MatrixXd matQ(n_states+1,n_states+1);
	  matQ.setZero();

	  for(unsigned int p=0; p<nPops; p++)
	    matQ = matQ + coeff4TransitionRateMat.at(i).at(j).at(p)/popSize(0,p);

	  unsigned int count = nPops;
	  for(unsigned int p1=0; p1< nPops; p1++)
	    {
	      for(unsigned int p2=0; p2< nPops; p2++)
		{
		  if(p1 != p2)
		    {
		      matQ = matQ + coeff4TransitionRateMat.at(i).at(j).at(count) * migRate(p1,p2);
		      count++;
		    }
		}
	    }

	  // diagonal elements
	  for(unsigned int r =0; r<n_states; r++)
	    matQ(r,r) = -matQ.row(r).sum();

	  // Matrix decomposition
	  Eigen::EigenSolver<Eigen::MatrixXd> svd;
	  svd.compute(matQ);
	  // eigen value and eigen vector matrices
	  Qeigen.at(i).at(j).at(0) = svd.eigenvalues();
	  //Qeigen.at(i).at(j).at(0) = svd.eigenvalues().asDiagonal();
	  Qeigen.at(i).at(j).at(1) = svd.eigenvectors();
	  Qeigen.at(i).at(j).at(2) = svd.eigenvectors().inverse();
	  // Checking eigen values and eigen vectors
	  Eigen::MatrixXcd I = Qeigen.at(i).at(j).at(1) * Qeigen.at(i).at(j).at(2);
	  if(abs(I.sum().real() - n_states-1) > pow(10,-4))
	    {
	      poptree->print_allPopTree();
	      std::cout << "*** Warning *** In Chain::compute_eigenValuesVectors()"
			<<" Eigen vectors may not be accurate. That is, V*V^(-1) != I.\n"
			<< "The sum of all elements in V*V^(-1) is " << I.sum().real()
			<< "\nThe sum of absolute differences is "<< abs(I.sum().real() - n_states-1) <<".\n";
	      std::cout << "matQ = " << matQ <<"\n";
	    }
	  Eigen::MatrixXcd rateMat2 = Qeigen.at(i).at(j).at(1) * svd.eigenvalues().asDiagonal() * Qeigen.at(i).at(j).at(2);
	  if(abs(rateMat2.sum().real() - matQ.sum()) > pow(10,-5))
	    {
	      poptree->print_allPopTree();
	      std::cout << "*** Warning *** In Chain::compute_eigenValuesVectors()"
			<< "Eigen values and eigen vectors may not be accurate. rateMat != V*D*V^(-1).\n"
			<<"The absolute difference of the sums of elements in rateMat and V*D*V^(-1) is "
			<< abs(rateMat2.sum().real() - matQ.sum()) <<"\n";
	      std::cout << "matQ = " << matQ <<"\n";
	      std::cout << "rateMat2 = " << rateMat2 <<"\n";
	    }
	  
	   

	}
    }

  return;
  
}











void Chain::share_treeIDs_coalTimes_logPriorTrees(unsigned int crrProc_ID)
{
  unsigned int nSample = 0;
  if(crrProc_ID ==0)
    {
      nSample = treeIDs.size();
    }
  MPI::COMM_WORLD.Bcast(&nSample, 1,  MPI::UNSIGNED, 0);
  if(crrProc_ID !=0 )
    {
      treeIDs.resize(nSample);
      coalTimes.resize(nSample);
      logPrior_trees.resize(nSample);
    }
  for(unsigned int ngen =0; ngen<nSample; ngen++)
    {
      MPI::COMM_WORLD.Bcast(&logPrior_trees.at(ngen), 1,  MPI::DOUBLE, 0);
      
      unsigned int nloci = 0;
      if(crrProc_ID ==0)
	{
	  nloci = coalTimes.at(ngen).size();
	}
      MPI::COMM_WORLD.Bcast(&nloci, 1,  MPI::UNSIGNED, 0);
      if(crrProc_ID != 0)
	{
	  treeIDs.at(ngen).resize(nloci);
	  coalTimes.at(ngen).resize(nloci);
	}
      for(unsigned int l=0; l<nloci; l++)
	{
	  MPI::COMM_WORLD.Bcast(&treeIDs.at(ngen).at(l), 1,  MPI::UNSIGNED, 0);
	  unsigned int nGeneCopies = 0;
	  if(crrProc_ID ==0)
	    {
	      nGeneCopies = coalTimes.at(ngen).at(l).size()+1;
	    }
	  MPI::COMM_WORLD.Bcast(&nGeneCopies, 1,  MPI::UNSIGNED, 0);
	  if(crrProc_ID !=0)
	    {
	      coalTimes.at(ngen).at(l).resize(nGeneCopies-1);
	    }
	  std::list<double>::iterator iter = coalTimes.at(ngen).at(l).begin();
	  for(unsigned int g=0; g<nGeneCopies-1; g++)
	    {
	      double ctime = 0.0;
	      if(crrProc_ID == 0)
		{
		  ctime = *iter;		      
		}
	      MPI::COMM_WORLD.Bcast(&ctime, 1,  MPI::DOUBLE, 0);
	      if(crrProc_ID != 0)
		{
		  *iter = ctime;
		}
	      ++iter;
	    }	    
	}
      
    }
  return;    
}


void Chain::share_stateSpaces_nKindsLineages_possiblePaths_numTrees(unsigned int crrProc_ID)
{
  unsigned int numUniqTopo = 0;
  if(crrProc_ID ==0)
    {
      numUniqTopo = stateSpaces.size();
    }
  MPI::COMM_WORLD.Bcast(&numUniqTopo, 1,  MPI::UNSIGNED, 0);
  if(crrProc_ID !=0)
    {
      nKinds_lineages.resize(numUniqTopo);
      stateSpaces.resize(numUniqTopo);
      possiblePaths.resize(numUniqTopo);
      numTrees.resize(numUniqTopo);
    }
  for(unsigned int tp = 0; tp<numUniqTopo; tp++)
    {
      MPI::COMM_WORLD.Bcast(&numTrees.at(tp), 1,  MPI::UNSIGNED, 0);
      
      unsigned int nGeneCopies = 0;
      if(crrProc_ID == 0)
	{
	  nGeneCopies = stateSpaces.at(tp).size()+1;
	}
      MPI::COMM_WORLD.Bcast(&nGeneCopies, 1,  MPI::UNSIGNED, 0);
      if(crrProc_ID !=0)
	{
	  nKinds_lineages.at(tp).resize(nGeneCopies-1);
	  stateSpaces.at(tp).resize(nGeneCopies-1);
	  possiblePaths.at(tp).resize(nGeneCopies-1);
	}
      for(unsigned int g=0; g<nGeneCopies-1; g++)
	{
	  MPI::COMM_WORLD.Bcast(&nKinds_lineages.at(tp).at(g), 1,  MPI::UNSIGNED, 0);
	  unsigned int nr = 0;
	  unsigned int nc = 0;
	  if(crrProc_ID == 0)
	    {
	      nr = stateSpaces.at(tp).at(g).rows();
	      nc = stateSpaces.at(tp).at(g).cols();
	    }
	  MPI::COMM_WORLD.Bcast(&nr, 1,  MPI::UNSIGNED, 0);
	  MPI::COMM_WORLD.Bcast(&nc, 1,  MPI::UNSIGNED, 0);
	  if(crrProc_ID !=0)
	    {
	      stateSpaces.at(tp).at(g).resize(nr,nc);
	      stateSpaces.at(tp).at(g).setZero();
	    }
	  for(unsigned int r = 0; r<nr;r++)
	    {
	      for(unsigned int c =0; c<nc; c++)
		{
		  MPI::COMM_WORLD.Bcast(&stateSpaces.at(tp).at(g)(r,c), 1,  MPI::DOUBLE, 0);
		  
		}
	    }
	  unsigned int sizeI=0;
	  if(crrProc_ID == 0)
	    {
	      sizeI = possiblePaths.at(tp).at(g).size();
	    }
	  MPI::COMM_WORLD.Bcast(&sizeI, 1,  MPI::UNSIGNED, 0);
	  if(crrProc_ID !=0)
	    {
	      possiblePaths.at(tp).at(g).resize(sizeI);
	    }
	  for(unsigned int i=0; i<sizeI; i++)
	    {
	      unsigned int sizeJ =0;
	      if(crrProc_ID == 0)
		{
		  sizeJ = possiblePaths.at(tp).at(g).at(i).size();
		}
	      MPI::COMM_WORLD.Bcast(&sizeJ, 1,  MPI::UNSIGNED, 0);
	      if(crrProc_ID !=0)
		{
		  possiblePaths.at(tp).at(g).at(i).resize(sizeJ);
		}
	      for(unsigned int j=0; j<sizeJ; j++)
		{
		  MPI::COMM_WORLD.Bcast(&possiblePaths.at(tp).at(g).at(i).at(j), 1,  MPI::UNSIGNED, 0);		  
		}   
	    }	  
	} 
    }
  return;
}


void Chain::share_stateSpaces_nKindsLineages_possiblePaths_numTrees_nPossibleCoalEvents(unsigned int crrProc_ID)
{
  unsigned int numUniqTopo = 0;
  if(crrProc_ID ==0)
    {
      numUniqTopo = stateSpaces.size();
    }
  MPI::COMM_WORLD.Bcast(&numUniqTopo, 1,  MPI::UNSIGNED, 0);

  if(crrProc_ID !=0)
    {
      nKinds_lineages.resize(numUniqTopo);
      stateSpaces.resize(numUniqTopo);
      possiblePaths.resize(numUniqTopo);
      numTrees.resize(numUniqTopo);
      nPossibleCoalEvents.resize(numUniqTopo);
    }
  for(unsigned int tp = 0; tp<numUniqTopo; tp++)
    {
      // Sharing numTrees
      MPI::COMM_WORLD.Bcast(&numTrees.at(tp), 1,  MPI::UNSIGNED, 0);
      
      unsigned int nGeneCopies = 0;
      if(crrProc_ID == 0)
	{
	  nGeneCopies = stateSpaces.at(tp).size()+1;
	}
      MPI::COMM_WORLD.Bcast(&nGeneCopies, 1,  MPI::UNSIGNED, 0);
      if(crrProc_ID !=0)
	{
	  nKinds_lineages.at(tp).resize(nGeneCopies-1);
	  stateSpaces.at(tp).resize(nGeneCopies-1);
	  possiblePaths.at(tp).resize(nGeneCopies-1);
	  nPossibleCoalEvents.at(tp).resize(nGeneCopies-1);
	}
      for(unsigned int g=0; g<nGeneCopies-1; g++)
	{
	  // Sharing nKinds_lineages
	  MPI::COMM_WORLD.Bcast(&nKinds_lineages.at(tp).at(g), 1,  MPI::UNSIGNED, 0);

	  // Sharing stateSpaces
	  unsigned int nr = 0;
	  unsigned int nc = 0;
	  if(crrProc_ID == 0)
	    {
	      nr = stateSpaces.at(tp).at(g).rows();
	      nc = stateSpaces.at(tp).at(g).cols();
	    }
	  MPI::COMM_WORLD.Bcast(&nr, 1,  MPI::UNSIGNED, 0);
	  MPI::COMM_WORLD.Bcast(&nc, 1,  MPI::UNSIGNED, 0);
	  if(crrProc_ID !=0)
	    {
	      stateSpaces.at(tp).at(g).resize(nr,nc);
	      stateSpaces.at(tp).at(g).setZero();
	    }
	  for(unsigned int r = 0; r<nr;r++)
	    {
	      for(unsigned int c =0; c<nc; c++)
		{
		  MPI::COMM_WORLD.Bcast(&stateSpaces.at(tp).at(g)(r,c), 1,  MPI::DOUBLE, 0);
		  
		}
	    }

	  // Sharing possiblePaths
	  unsigned int sizeI=0;
	  if(crrProc_ID == 0)
	    {
	      sizeI = possiblePaths.at(tp).at(g).size();
	    }
	  MPI::COMM_WORLD.Bcast(&sizeI, 1,  MPI::UNSIGNED, 0);
	  if(crrProc_ID !=0)
	    {
	      possiblePaths.at(tp).at(g).resize(sizeI);
	    }
	  for(unsigned int i=0; i<sizeI; i++)
	    {
	      unsigned int sizeJ =0;
	      if(crrProc_ID == 0)
		{
		  sizeJ = possiblePaths.at(tp).at(g).at(i).size();
		}
	      MPI::COMM_WORLD.Bcast(&sizeJ, 1,  MPI::UNSIGNED, 0);
	      if(crrProc_ID !=0)
		{
		  possiblePaths.at(tp).at(g).at(i).resize(sizeJ);
		}
	      for(unsigned int j=0; j<sizeJ; j++)
		{
		  MPI::COMM_WORLD.Bcast(&possiblePaths.at(tp).at(g).at(i).at(j), 1,  MPI::UNSIGNED, 0);		  
		}   
	    }
	  
	  // Sharing nPossibleCoalEvents
	  unsigned int nAllPops = 0;
	  if(crrProc_ID == 0)
	    {
	      nAllPops = nPossibleCoalEvents.at(tp).at(g).size();
	    }
	  MPI::COMM_WORLD.Bcast(&nAllPops, 1,  MPI::UNSIGNED, 0);
	  if(crrProc_ID !=0)
	    {
	      nPossibleCoalEvents.at(tp).at(g).resize(nAllPops);
	    }
	  for(unsigned int p=0; p<nAllPops; p++)
	    {
	      unsigned int sizeC = 0;
	      if(crrProc_ID == 0)
		{
		  sizeC = nPossibleCoalEvents.at(tp).at(g).at(p).size();
		}
	      MPI::COMM_WORLD.Bcast(&sizeC, 1,  MPI::UNSIGNED, 0);	      
	      if(crrProc_ID !=0)
		{
		  nPossibleCoalEvents.at(tp).at(g).at(p).resize(sizeC);
		}
	      for(unsigned int j=0; j<sizeC; j++)
		{
		  MPI::COMM_WORLD.Bcast(&nPossibleCoalEvents.at(tp).at(g).at(p).at(j), 1,  MPI::UNSIGNED, 0);		  
		}   	      
	    }
	  // END of sharing nPossibleCoalEvents         
	} // END of for(unsigned int g=0; g<nGeneCopies-1; g++)
    } // END of for(unsigned int tp = 0; tp<numUniqTopo; tp++)

  return;
}



void Chain::share_eigenValuesVectors(unsigned int crrProc_ID)
{
  unsigned int numUniqTopo = 0;
  if(crrProc_ID == 0)
    {
      numUniqTopo = Qeigen.size();
    }
  MPI::COMM_WORLD.Bcast(&numUniqTopo, 1,  MPI::UNSIGNED, 0);
  if(crrProc_ID !=0)
    {
      Qeigen.resize(numUniqTopo);
    }
  for(unsigned int i=0; i<numUniqTopo; i++)
    {
      unsigned int nGeneCopies = 0;
      if(crrProc_ID == 0)
	{	  
	  nGeneCopies = Qeigen.at(i).size()+1;
	}
      MPI::COMM_WORLD.Bcast(&nGeneCopies, 1,  MPI::UNSIGNED, 0);
      if(crrProc_ID != 0)
	{
	  Qeigen.at(i).resize(nGeneCopies-1);
	}
      for(unsigned int j=0; j<nGeneCopies-1; j++)
	{
	  if(crrProc_ID != 0)
	    {
	      Qeigen.at(i).at(j).resize(3);
	    }
	  for(unsigned int k=0; k<3; k++)
	    {
	      unsigned int nrow = 0;
	      unsigned int ncol = 0;
	      if(crrProc_ID == 0)
		{
		  nrow = Qeigen.at(i).at(j).at(k).rows();
		  ncol = Qeigen.at(i).at(j).at(k).cols();
		}
	      MPI::COMM_WORLD.Bcast(&nrow, 1,  MPI::UNSIGNED, 0);
	      MPI::COMM_WORLD.Bcast(&ncol, 1,  MPI::UNSIGNED, 0);
	      if(crrProc_ID != 0)
		{
		  Qeigen.at(i).at(j).at(k).resize(nrow,ncol);
		  Qeigen.at(i).at(j).at(k).setZero();
		}
	      for(unsigned int r=0; r<nrow; r++)
		{
		  for(unsigned int c =0; c<ncol; c++)
		    {
		      MPI::COMM_WORLD.Bcast(&Qeigen.at(i).at(j).at(k)(r,c), 1,  MPI::DOUBLE, 0);
		    }
		}
	    }
	}
    }
  return;
}



// old version
double Chain::compute_conditionalProb(unsigned int id_sample, unsigned int id_locus, popTree* poptree)
{
  unsigned int nPops = poptree->size();

  // FIXME YC 5/9/2014
  // It works for up to 2 population
  // Note: It might be better to have a migration rate matrix and a vector of population sizes
  // as members of class 'popTree'.
  Eigen::MatrixXd popSize(1,nPops);
  popSize.setZero();
  for(unsigned int p1=0; p1< nPops; p1++)
    {
      popSize(0,p1) = poptree->find_popSize(p1+1);
    }
  double splittingTime = poptree->get_age();
  
  double prob = 1;
  unsigned int trID = treeIDs.at(id_sample).at(id_locus);
  std::list<double> eventT = coalTimes.at(id_sample).at(id_locus);
  eventT.push_back(splittingTime);
  eventT.sort();
  
  
  std::list<double>::iterator iter = eventT.begin();
  
  // REMOVE
  // std::cout << "id_split = " << id_split <<"\n";
  //for(iter = eventT.begin(); iter !=eventT.end(); ++iter)
  //	std::cout << *iter <<" ";
  //std::cout <<"\n";
  
  unsigned int nGeneCopies = coalTimes.at(id_sample).at(id_locus).size()+1; //nKinds_lineages.at(trID).size()+1; //list_trees.at(trID)->getSize();
  unsigned int ancestralPop = 0;
  Eigen::MatrixXcd probMat;
  unsigned int i=0;
  iter = eventT.begin();
  while(i<nGeneCopies-1 && iter!=eventT.end() && ancestralPop == 0 )
    {
      if(*iter > splittingTime)
	{
	  
	  double popSize = poptree->get_popSize();
	  unsigned int nLineages = nGeneCopies - i;
	  double expterm = 0.0;
	  prob  *= pow(2/popSize,nLineages-1);
			for(unsigned int l=nLineages; l>=2; l--, ++iter)
			{
				expterm += l*(l-1)* (*iter)*2/popSize;
			}
			prob *= exp(-1*expterm);
			ancestralPop = 1;
		}
		else if(*iter == splittingTime)
		{
			double waitingTime = *iter;
			Eigen::MatrixXcd V;
			Eigen::MatrixXcd V_inv;
			if(i==0)
			{
				V = Qeigen.at(trID).at(i).at(1).row(possiblePaths.at(trID).at(i).at(0).at(0));
				V_inv = Qeigen.at(trID).at(i).at(2);

				//--- Get exp(D*t) ---//
				Eigen::MatrixXcd expD = Qeigen.at(trID).at(i).at(0) *waitingTime;
				for(unsigned int di=0; di<expD.rows(); di++)
				{
					complex<double> val = expD(di,di);
					expD(di,di) = exp(val);
				}

				probMat = V*expD*V_inv;
			}
			else
			{
				V_inv = Qeigen.at(trID).at(i).at(2);

				unsigned int nr = 0;
				for(unsigned int p=0; p<nPops; p++)
				{
					nr += possiblePaths.at(trID).at(i).at(p).size();
				}
				V.resize(nr, Qeigen.at(trID).at(i).at(1).rows());
				unsigned int count_r = 0;
				for(unsigned int p=0; p<nPops; p++)
				{
					for(unsigned int k=0; k<possiblePaths.at(trID).at(i).at(p).size(); k++)
					{
						unsigned int id_row = possiblePaths.at(trID).at(i).at(p).at(k);
						V.row(count_r) =Qeigen.at(trID).at(i).at(1).row(id_row);
						count_r++;
					}
				}

				//--- Get exp(D*t) ---//
				Eigen::MatrixXcd expD = Qeigen.at(trID).at(i).at(0) *waitingTime;
				for(unsigned int di=0; di<expD.rows(); di++)
				{
					complex<double> val = expD(di,di);
					expD(di,di) = exp(val);
				}
				probMat *= V*expD*V_inv;
			}
			complex<double> p_tmp = probMat.sum();
			prob = p_tmp.real();
		}
		else
		{
			double waitingTime = *iter;
			Eigen::MatrixXcd V;
			Eigen::MatrixXcd V_inv;
			if(i==0)
			{
				unsigned int nc = 0;
				for(unsigned int p=0; p<nPops; p++)
					nc += possiblePaths.at(trID).at(i).at(p+1).size();

				V = Qeigen.at(trID).at(i).at(1).row(possiblePaths.at(trID).at(i).at(0).at(0));
				V_inv.resize(Qeigen.at(trID).at(i).at(1).rows(), nc);
				unsigned int count = 0;
				for(unsigned int p=0; p<nPops; p++)
				{
					for(unsigned int k=0; k<possiblePaths.at(trID).at(i).at(p+1).size(); k++)
					{
						unsigned int id_col = possiblePaths.at(trID).at(i).at(p+1).at(k);
						V_inv.col(count) =Qeigen.at(trID).at(i).at(2).col(id_col)* 2/popSize(0,p);
						count++;
					}
				}

				//--- Get exp(D*t) ---//
				Eigen::MatrixXcd expD = Qeigen.at(trID).at(i).at(0) *waitingTime;
				for(unsigned int di=0; di<expD.rows(); di++)
				{
					complex<double> val = expD(di,di);
					expD(di,di) = exp(val);
				}

				probMat = V*expD*V_inv;

				// REMOVE:
				/*
			std::cout << "Coal time is " << *iter <<"\n";
			std::cout << "The full transition prob mat is\n" <<Qeigen.at(trID).at(i).at(1)*expD*Qeigen.at(trID).at(i).at(2) <<"\n";
			std::cout << "V = \n" << V <<"\n"
					<< "V_inv = \n"<< V_inv <<"\n"
					<< "Full V =\n" << Qeigen.at(trID).at(i).at(1) <<"\n"
					<< "Full V_inv =\n" << Qeigen.at(trID).at(i).at(2) <<"\n";
 			std::cout << "V*expD*V_inv is\n" <<  V*expD*V_inv <<"\n";
			std::cout << "probMat is" <<"\n";
			std::cout << probMat <<"\n\n";
			*/

			}
			else
			{
				unsigned int nr = 0;
				unsigned int nc = 0;
				for(unsigned int p=0; p<nPops; p++)
				{
					nr += possiblePaths.at(trID).at(i).at(p).size();
					nc += possiblePaths.at(trID).at(i).at(nPops+p).size();
				}

				V.resize(nr, Qeigen.at(trID).at(i).at(1).rows());
				V_inv.resize(Qeigen.at(trID).at(i).at(1).rows(), nc);
				unsigned int count_r = 0;
				unsigned int count_c = 0;
				for(unsigned int p=0; p<nPops; p++)
				{
					for(unsigned int k=0; k<possiblePaths.at(trID).at(i).at(p).size(); k++)
					{
						unsigned int id_row = possiblePaths.at(trID).at(i).at(p).at(k);
						V.row(count_r) =Qeigen.at(trID).at(i).at(1).row(id_row);
						count_r++;
					}
					for(unsigned int k=0; k<possiblePaths.at(trID).at(i).at(nPops+p).size(); k++)
					{
						unsigned int id_col = possiblePaths.at(trID).at(i).at(nPops+p).at(k);
						V_inv.col(count_c) =Qeigen.at(trID).at(i).at(2).col(id_col) * 2/popSize(0,p);
						count_c++;
					}
				}

				//--- Get exp(D*t) ---//
				Eigen::MatrixXcd expD = Qeigen.at(trID).at(i).at(0) *waitingTime;
				for(unsigned int di=0; di<expD.rows(); di++)
				{
					complex<double> val = expD(di,di);
					expD(di,di) = exp(val);
				}
				probMat *= V*expD*V_inv;


			// REMOVE:
			//std::cout << "Coal time is " << *iter <<"\n";
			//std::cout << "The full transition prob mat is\n" <<Qeigen.at(trID).at(i).at(1)*expD*Qeigen.at(trID).at(i).at(2) <<"\n";
			//std::cout << "V = \n" << V <<"\n"
								//<< "V_inv = \n"<< V_inv <<"\n"
								//<< "Full V =\n" << Qeigen.at(trID).at(i).at(1) <<"\n"
								//<< "Full V_inv =\n" << Qeigen.at(trID).at(i).at(2) <<"\n";
			//std::cout << "V*expD*V_inv is\n" <<  V*expD*V_inv <<"\n";
			//std::cout << "probMat is" <<"\n";
			//std::cout << probMat <<"\n\n";
			}
			i++;
		}
		++iter;
	} // End of while loop


	if(ancestralPop ==0)
	{
		complex<double> p_tmp = probMat(0,0)+probMat(0,1);
		prob = p_tmp.real();
	}

	return prob;
}


// YC 3/25/2015
// 'crrProcID' is required for debugging only. This function is called on each process.
double Chain::compute_logConditionalProb(unsigned int id_sample, unsigned int id_locus, popTree* poptree, unsigned int crrProcID)
{
  std::chrono::high_resolution_clock::time_point start_t, end_t;
  std::chrono::high_resolution_clock::time_point start_t_case, end_t_case;
  start_t= std::chrono::high_resolution_clock::now();


  unsigned int nPops = poptree->size();
  // std::cout << "nPops = " << nPops <<"\n";
  
  // FIXME YC 5/9/2014
  // It works for up to 2 population
  // Note: It might be better to have a migration rate matrix and a vector of population sizes
  // as members of class 'popTree'.
  Eigen::MatrixXd popSize(1,nPops);
  popSize.setZero();
  for(unsigned int p1=0; p1< nPops; p1++)
    {
      popSize(0,p1) = poptree->find_popSize(p1+1);
    }
  double splittingTime = poptree->get_age(); // Sorts the elements in ascending order

  // std::cout << "popSize = " << popSize <<"\n";
  // std::cout << "splittingTime = " << splittingTime <<"\n";

  double logProb = 0;

  unsigned int trID = treeIDs.at(id_sample).at(id_locus);
  // std::cout << "trID = " << trID <<"\n";
  std::list<double> eventT = coalTimes.at(id_sample).at(id_locus);


  eventT.sort();

  // std::cout <<"nGeneCopies = " << coalTimes.at(id_sample).at(id_locus).size()+1 <<"\n";
  
  std::list<double>::iterator iter = eventT.begin();
  
  unsigned int nGeneCopies = coalTimes.at(id_sample).at(id_locus).size()+1; 
  unsigned int ancestralPop = 0;
  Eigen::MatrixXcd probMat;
  unsigned int count_events=0;
  iter = eventT.begin();
  while(count_events <nGeneCopies-1 && iter!=eventT.end())
    {  
      // std::cout << "*iter = " << *iter << "\n";
      if(*iter <=splittingTime)
	{	  
	  //--- Case 3: coalescent events in sampling populations ---//
	  // std::cout << "case 3\n";

	  Eigen::MatrixXcd V = subMatV.at(trID).at(count_events).at(0);
	  Eigen::ArrayXcd eigenval = eigenvalues.at(trID).at(count_events);
	  std::list<double>::iterator iter_prev = iter;
	  
	  double waitingTime = 0.0;
	  if(count_events==0)
	    {
	      waitingTime = *iter;
	    }
	  else
	    {		      
	      iter_prev--;		 
	      waitingTime = *iter - *iter_prev;
	    }
	  
	  Eigen::MatrixXcd V_inv = subMatV.at(trID).at(count_events).at(1);
	  Eigen::MatrixXcd probMat_each;
	  
	  if(count_events==0) // the first coalescent event
	    {
	      //std::cout << "Case3-1 : count_events = "<< count_events<< "\n";
	      count_case3_1++;	      
	      start_t_case= std::chrono::high_resolution_clock::now();
	      
	      //--- Computing V*exp(D*t) ---//
	      V = V.cwiseProduct((eigenval*waitingTime).exp().matrix().transpose());
	      probMat = V*V_inv;
	      
	      end_t_case= std::chrono::high_resolution_clock::now();	      
	      eachComputingTime_condiProb_case3_1 += std::chrono::duration_cast<std::chrono::microseconds>(end_t_case - start_t_case);
	    }
	  else // non-first coalescent event
	    {
	      //std::cout << "Case3-2 : count_events = "<< count_events<< "\n";
	      count_case3_2++;
	      start_t_case= std::chrono::high_resolution_clock::now();

	      //Eigen::MatrixXcd expD = (eigenval*waitingTime).exp().matrix().asDiagonal();
	      if(probMat.rows() != 1 || probMat.cols() != V.rows())
		{
		  std::cout << "\n ***** Error in Chain::compute_conditionalProb_MPI *****\n";	
		  std::cout << "\n ***** Case 3-2: different dimensions *****\n";
		  std::cout << "\n ***** id_sample = "<< id_sample
			    << ", id_locus = "<< id_locus <<", crrProcID = "<< crrProcID << " and poptree is ";
		  poptree->print_allPopTree();
		}
	      probMat = probMat*V; // row-vector
	      probMat = probMat.cwiseProduct((eigenval*waitingTime).exp().matrix().transpose());
	      probMat = probMat*V_inv;	      
	      
	      end_t_case= std::chrono::high_resolution_clock::now();	      
	      eachComputingTime_condiProb_case3_2 += std::chrono::duration_cast<std::chrono::microseconds>(end_t_case - start_t_case);
	    }

	  //std::cout << "Case 3: *iter = " << *iter << " splittingtime = " << splittingTime
	  //	    << " waitingTime= " << waitingTime << " and probMat=" << probMat <<"\n";

	  count_events++; // counting the number of events that happened before the splitting time.
	  	  
	  // Added by YC 3/8/2015
	  // During optimization (differential evolution), proposed splitting times
	  // can be same as coalescent times. 
	  // If the splitting time and a coalescent time are the same, then
	  // we assume the coalescent event happened in a sampling population 
	  // before the splitting event backward in time.
	  // If the splitting time is same as a coalescent time smaller than
	  // the TMRCA (tree height) (i.e., count_events <nGeneCopies-1), then
	  // we have a flag of 'ancestralPop = 1'. 
	  // If ancestralPop =1, then we do NOT need to compute the probability of
	  // no coalescent events between the last observed coalescent event time in 
	  // in sampling populations and the splitting time and 'Case 1' is called. 
	  if(*iter == splittingTime  && count_events <nGeneCopies-1)
	    {
	      ancestralPop = 1;
	    }
	  
	  ++iter;

	}
      else if(*iter > splittingTime && ancestralPop ==0 && splittingTime!=0)
	{
	  //--- Case 2: No coalescent event between the last coalescent ---//
	  //--- event in the sampling populations and the splitting time --//
	  
	  // std::cout << "Case 2\n";
	  
	  Eigen::MatrixXcd V = subMatV.at(trID).at(count_events).at(0);
	  Eigen::ArrayXcd eigenval = eigenvalues.at(trID).at(count_events);

	  std::list<double>::iterator iter_prev = iter;
	  
	  double waitingTime = 0.0;
	  if(count_events==0)
	    {
	      waitingTime = splittingTime;
	    }
	  else
	    {		      
	      iter_prev--; 
	      waitingTime = splittingTime - *iter_prev;
	    }	  
	  
	  Eigen::MatrixXcd V_inv = subMatV.at(trID).at(count_events).at(2);
	  Eigen::MatrixXcd probMat_each;
	  if(count_events==0)
	    {
	      count_case2_1++;
	      start_t_case= std::chrono::high_resolution_clock::now();
	      	      
	      //--- Get exp(D*t) ---//
	      V=V.cwiseProduct((eigenval*waitingTime).exp().matrix().transpose());

	      if(V.cols() != V_inv.rows())
		{
		  std::cout << "\n ***** Error in Chain::compute_conditionalProb_MPI *****\n";	
		  std::cout << "\n ***** Case 2-1: different dimensions *****\n";
		  std::cout << "trID = " << trID <<"\n";
		  std::cout << "possiblePaths.at(trID).at(count_events).at(0).at(0) ="<<possiblePaths.at(trID).at(count_events).at(0).at(0)<<"\n";
		  std::cout << "V = " << V <<"\n"
			    <<"V_inv = " << V_inv<<"\n"
			    <<"eigenval = " <<eigenval <<"\n";
		}
	      else
		{
		  probMat = V*V_inv;
		}

	      end_t_case= std::chrono::high_resolution_clock::now();	      
	      eachComputingTime_condiProb_case2_1 += std::chrono::duration_cast<std::chrono::microseconds>(end_t_case - start_t_case);
	    }
	  else // if(count_events != 0 && waitingTime > 0.0)
	    {
	      // std::cout << "Case2-2 : count_events = "<< count_events<< "\n";
	      count_case2_2++;
	      start_t_case= std::chrono::high_resolution_clock::now();

	      
	      if(probMat.rows()!=1 || probMat.cols() != V.rows())
		{
		  std::cout << "\n ***** Error in Chain::compute_conditionalProb_MPI *****\n";	
		  std::cout << "\n ***** Case 2-2: different dimensions *****\n";
		}
	      probMat *= V; // row-vector
	      probMat = probMat.cwiseProduct((eigenval*waitingTime).exp().matrix().transpose());
	      probMat *= V_inv;	   

	      end_t_case= std::chrono::high_resolution_clock::now();	      
	      eachComputingTime_condiProb_case2_2 += std::chrono::duration_cast<std::chrono::microseconds>(end_t_case - start_t_case);
	    }
	  
	  double prob_noAbsorbing = probMat.sum().real();
	  /*
	  for(unsigned int id_No_AbsorbingState = 0; id_No_AbsorbingState < probMat.cols()-1; id_No_AbsorbingState++)
	    prob_noAbsorbing += probMat(0,id_No_AbsorbingState).real();
	  */

	  logProb += log(prob_noAbsorbing);
	  
	  // The flag that 'Case 1' should be called.
	  ancestralPop =1;
	}
      else if(splittingTime==0 || *iter > splittingTime && ancestralPop == 1)
	{	  
	  //--- Case 1: coalescents in the ancestral population. ---//
	  //  Some or all coalescent events happened in the ancestral population.
	  
	  // std::cout << "case 1\n";

	  count_case1++;
	  start_t_case= std::chrono::high_resolution_clock::now();
  
	  double popSize = poptree->get_popSize();
	  unsigned int nLineages = nGeneCopies - count_events; // the number of remaining lineages
	  double expterm = 0.0;
	  for(unsigned int l=nLineages; l>=2; l--, ++iter)
	    {
	      std::list<double>::iterator iter_prev = iter;
	      double waitingTime =0.0;
	      if(l == nLineages) // the first coalescent event in the ancestral population
		waitingTime = *iter - splittingTime;
	      else
		{
		  iter_prev--;
		  waitingTime = *iter - *iter_prev;
		}

	      expterm += l*(l-1)* waitingTime / popSize;
	    }
	  logProb += (nLineages-1)*log(2/popSize) - expterm;

	  // YC 4/29/2015
	  // Computing the term of the number of same coalescent events
	  unsigned int trID = treeIDs.at(id_sample).at(id_locus);
	  double log_numSameCoalEvents =0;
	  for(unsigned int c=count_events; c < nGeneCopies-2; c++)
	    {
	      log_numSameCoalEvents += log(static_cast<double> (numSameCoalEvents_inAncPop.at(trID).at(c)) );
	    }
	  logProb += log_numSameCoalEvents;

	  count_events = nGeneCopies-1;
	  ancestralPop = 1;

	  end_t_case= std::chrono::high_resolution_clock::now();	      
	  eachComputingTime_condiProb_case1 += std::chrono::duration_cast<std::chrono::microseconds>(end_t_case - start_t_case);
	}
    } // End of while loop
  
  
  if(ancestralPop ==0)
    {

      // std::cout << "Case4 : count_events = "<< count_events<< "\n";

      complex<double> p_tmp;
      if(probMat.rows() != 1 || probMat.cols() != 2)
	{
	  std::cout << "\n ***** Error in Chain::compute_conditionalProb_MPI *****\n";	
	  std::cout << "\n ***** Case 4: different dimensions *****\n";
	  std::cout << "\n ***** id_sample = "<< id_sample
		    << ", id_locus = "<< id_locus <<", crrProcID = "<< crrProcID << " and poptree is ";
	  poptree->print_allPopTree();
	  std::cout << "The probMat should be 1x2 matrix, but ";
	  std::cout << "probMat = " << probMat <<"\n";
	}
      else
	{
	  p_tmp = probMat(0,0)+probMat(0,1);
	}
      double p_real = p_tmp.real();
      if(p_real <0)
	{
	  p_real *= -1;
	  if(p_real < pow(10,-10))
	    p_real = 0;
	  else
	    {
	      /*
	      std::cout << "\n***Error in Chain::compute_logConditionalProb()\n";
	      std::cout << "p_tmp.real() = " << p_tmp.real() <<"\n";
	      poptree->print_allPopTree();
	      */
	    }
	  p_real = pow(10,-10);
	}

      logProb = log(p_real);
	 
    }

  end_t= std::chrono::high_resolution_clock::now();
  eachComputingTime_condiProb += std::chrono::duration_cast<std::chrono::microseconds>(end_t - start_t);
  count_condiProbFunctionCalls++;

  // std::cout << "logProb = " << logProb <<"\n";

  return logProb;
}


// YC 1/14/2015
// 'crrProcID' is required for debugging only. This function is called on each process.
/*
double Chain::compute_logConditionalProb_old(unsigned int id_sample, unsigned int id_locus, popTree* poptree, unsigned int crrProcID)
{
  unsigned int nPops = poptree->size();

  // std::cout << "nPops = " << nPops <<"\n";
  
  // FIXME YC 5/9/2014
  // It works for up to 2 population
  // Note: It might be better to have a migration rate matrix and a vector of population sizes
  // as members of class 'popTree'.
  Eigen::MatrixXd popSize(1,nPops);
  popSize.setZero();
  for(unsigned int p1=0; p1< nPops; p1++)
    {
      popSize(0,p1) = poptree->find_popSize(p1+1);
    }
  double splittingTime = poptree->get_age(); // Sorts the elements in ascending order


  double logProb = 0;

  unsigned int trID = treeIDs.at(id_sample).at(id_locus);
  std::list<double> eventT = coalTimes.at(id_sample).at(id_locus);


  eventT.push_back(splittingTime);
  eventT.sort();

  
  std::list<double>::iterator iter = eventT.begin();
  
  unsigned int nGeneCopies = coalTimes.at(id_sample).at(id_locus).size()+1; //nKinds_lineages.at(trID).size()+1; //list_trees.at(trID)->getSize();
  unsigned int ancestralPop = 0;
  Eigen::MatrixXcd probMat;
  unsigned int count_events=0;
  iter = eventT.begin();
  while(count_events <nGeneCopies-1 && iter!=eventT.end() && ancestralPop == 0 )
    {

      // std::cout << "*iter = " << *iter <<"\n";
      
      if(*iter > splittingTime)
	{
	  // case 1: coalescent in the ancestral population
	  // std::cout << "case 1\n";
	  
	  double popSize = poptree->get_popSize();
	  unsigned int nLineages = nGeneCopies - count_events; // remaining number of lineages
	  double expterm = 0.0;
	  for(unsigned int l=nLineages; l>=2; l--, ++iter)
	    {
	      // Added by YC 6/11/2014
	      std::list<double>::iterator iter_prev = iter;
	      iter_prev--;
	      double waitingTime = *iter - *iter_prev;
	      
	      // Modified by YC 6/10/2014
	      expterm += l*(l-1)* waitingTime / popSize;
	      // unsigned int nCoalEvents = nPossibleCoalEvents.at(trID).at(nGeneCopies-l).at(nPops).at(0);
	      // prob  *= nCoalEvents * 2/popSize;
	      // logProb += log(nCoalEvents * 2/popSize);
	      //logProb += nCoalEvents * log(2/popSize);
	    }
	  logProb += (nLineages-1)*log(2/popSize) - expterm;
	  // double prob_beforeExp = prob;
	  //prob *= exp(-1*expterm);
	  // REMOVE
	  //if(prob == 0)
	  //  {
	  //    std::cout << "prob = " << prob << " coal. in the ancetral pop on procID "<< crrProcID<< "(id_sample ="<< id_sample <<", id_locus="<<id_locus<<")"<<"expterm = "<<expterm <<"prob_beforeExp="<< prob_beforeExp<<"\n";
	  //  }

	  // END of REMOVE
	  ancestralPop = 1;
	}
      else if(*iter == splittingTime) 
	{
	  // -----------------------------------------------------------------//
	  // If all the coalescent events happened before the splitting time (backward in time),
	  // this part won't be called. 
	  // This part is called if and only if some or all coalescent events
	  // happened after the splitting time or at the same time of the split.
	  // Therefore, in the next round, the case of (*iter > splittingTime) 
	  // should be called and "ancestralPop = 1" will be set in that case.
	  // ----------------------------------------------------------------//

	  //std::cout << "Case2 : count_events = "<< count_events<< "\n";

	  // Added by YC 6/11/2014
	  std::list<double>::iterator iter_prev = iter;
	  
	  double waitingTime = 0.0;
	  if(count_events==0)
	    {
	      waitingTime = *iter;
	    }
	  else
	    {		      
	      iter_prev--; 
	      waitingTime = *iter - *iter_prev;
	    }
	  
	  
	  if(waitingTime > 0.0)
	    {
	      Eigen::MatrixXcd V;
	      Eigen::MatrixXcd V_inv;
	      if(count_events==0)
		{

		  V = Qeigen.at(trID).at(count_events).at(1).row(possiblePaths.at(trID).at(count_events).at(0).at(0));
		  V_inv = Qeigen.at(trID).at(count_events).at(2);

		  //std::cout << "V = " << V <<"\n"
		  //	    << "V_inv = " << V_inv <<"\n";

		  //--- Get exp(D*t) ---//
		  Eigen::MatrixXcd expD = Qeigen.at(trID).at(count_events).at(0) *waitingTime;
		  for(unsigned int di=0; di<expD.rows(); di++)
		    {
		      complex<double> val = expD(di,di);
		      expD(di,di) = exp(val);
		    }
		  if(V.cols() != expD.rows() || expD.cols() != V_inv.rows())
		    {
		      std::cout << "\n ***** Error in Chain::compute_conditionalProb_MPI *****\n";	
		      std::cout << "\n ***** Case 2-1: different dimensions *****\n";
		      std::cout << "trID = " << trID <<"\n";
		      std::cout << "possiblePaths.at(trID).at(count_events).at(0).at(0) ="<<possiblePaths.at(trID).at(count_events).at(0).at(0)<<"\n";
		      std::cout << "V = " << V <<"\n"
				<<"V_inv = " << V_inv<<"\n"
				<<"expD = " <<expD <<"\n";
		    }
		  else
		    {
		      probMat = V*expD*V_inv;
		    }
		  // std::cout << "probMat = " << probMat <<"\n";
		}
	      else if(waitingTime > 0.0) // if(count_events != 0 && waitingTime > 0.0)
		{
		  //std::cout << "Case2-2 : count_events = "<< count_events<< "\n";
		  
		  V_inv = Qeigen.at(trID).at(count_events).at(2);
		  
		  unsigned int nr = 0;
		  for(unsigned int p=0; p<nPops; p++)
		    {
		      nr += possiblePaths.at(trID).at(count_events).at(p).size();
		    }
		  V.resize(nr, Qeigen.at(trID).at(count_events).at(1).rows());
		  unsigned int count_r = 0;
		  for(unsigned int p=0; p<nPops; p++)
		    {
		      for(unsigned int k=0; k<possiblePaths.at(trID).at(count_events).at(p).size(); k++)
			{
			  unsigned int id_row = possiblePaths.at(trID).at(count_events).at(p).at(k);
			  V.row(count_r) =Qeigen.at(trID).at(count_events).at(1).row(id_row);
			  count_r++;
			}
		    }
		  
		  //--- Get exp(D*t) ---//
		  Eigen::MatrixXcd expD = Qeigen.at(trID).at(count_events).at(0) *waitingTime;
		  for(unsigned int di=0; di<expD.rows(); di++)
		    {
		      complex<double> val = expD(di,di);
		      expD(di,di) = exp(val);
		    }
		  try
		    {
		      if(probMat.cols() != V.rows() || V.cols() != expD.rows() || expD.cols() != V_inv.rows())
			{
			  std::cout << "\n ***** Error in Chain::compute_conditionalProb_MPI *****\n";	
			  std::cout << "\n ***** Case 2-2: different dimensions *****\n";
			}
		      else
			{
			  probMat *= V*expD*V_inv;
			}
		    }
		  catch(std::exception &e)
		    {
		      std::cout << "\n ***** Error in Chain::compute_conditionalProb_MPI *****\n";		  
		    }
		}

	      double prob_noAbsorbing = 0.0;
	      for(unsigned int id_No_AbsorbingState = 0; id_No_AbsorbingState < probMat.cols()-1; id_No_AbsorbingState++)
		prob_noAbsorbing += probMat(0,id_No_AbsorbingState).real();
	      // prob *= prob_noAbsorbing;
	      logProb += log(prob_noAbsorbing);


	  
	    }
	  else if(waitingTime == 0.0 && count_events!=0)
	    {
	      //std::cout << "case 2-3 : count_events = "<< count_events<< "\n";

	      // ---------------------------------------------------------//
	      // YC 7/14/2014
	      // When a coalescent event and splitting happened at the same time,
	      // we assume the coalescent event happened in the ancesteral population.
	      // ---------------------------------------------------------//

	      
	      logProb += log(2/poptree->get_popSize());
	      ancestralPop =1;
	      count_events++;
	    }
	 
	
	}
      else // if (*iter < splittingTime)
	{	  
	  // Added by YC 6/11/2014
	  std::list<double>::iterator iter_prev = iter;
	  
	  double waitingTime = 0.0;
	  if(count_events==0)
	    {
	      waitingTime = *iter;
	    }
	  else
	    {		      
	      iter_prev--;		 
	      waitingTime = *iter - *iter_prev;
	    }

	  Eigen::MatrixXcd V;
	  Eigen::MatrixXcd V_inv;
	  //Eigen::MatrixXcd V_inv_test;
	  //Eigen::MatrixXcd Coal_test;
	  
	  if(count_events==0)
	    {
	      //std::cout << "Case3-1 : count_events = "<< count_events<< "\n";

	      unsigned int nc = 0;
	      for(unsigned int p=0; p<nPops; p++)
		nc += possiblePaths.at(trID).at(count_events).at(p+1).size();
	      
	      V = Qeigen.at(trID).at(count_events).at(1).row(possiblePaths.at(trID).at(count_events).at(0).at(0));
	      V_inv.resize(Qeigen.at(trID).at(count_events).at(1).rows(), nc);
	      unsigned int count = 0;
	      for(unsigned int p=0; p<nPops; p++)
		{
		  for(unsigned int k=0; k<possiblePaths.at(trID).at(count_events).at(p+1).size(); k++)
		    {
		      unsigned int id_col = possiblePaths.at(trID).at(count_events).at(p+1).at(k);
		      unsigned int nCoalEvents = nPossibleCoalEvents.at(trID).at(count_events).at(p).at(k);
		      //std::cout << "nPossibleCoalEvents = " << nCoalEvents <<"\n";
		      V_inv.col(count) =Qeigen.at(trID).at(count_events).at(2).col(id_col)* nCoalEvents *2/popSize(0,p);
		      count++;
		    }
		}
	      
	      //--- Get exp(D*t) ---//
	      Eigen::MatrixXcd expD = Qeigen.at(trID).at(count_events).at(0) *waitingTime;
	      for(unsigned int di=0; di<expD.rows(); di++)
		{
		  complex<double> val = expD(di,di);
		  expD(di,di) = exp(val);
		}
	      
	      if(V.cols() != expD.rows() || expD.cols() != V_inv.rows())
		{
		  std::cout << "\n ***** Error in Chain::compute_conditionalProb_MPI *****\n";	
		  std::cout << "\n ***** Case 3-1: different dimensions *****\n";		  
		  std::cout << "\n ***** id_sample = "<< id_sample
			    << ", id_locus = "<< id_locus <<", crrProcID = "<< crrProcID << " and poptree is ";
		  poptree->print_allPopTree();
		}
	      else
		{
		  probMat = V*expD*V_inv;
		}
	      
		      
	    }
	  else
	    {
	      //std::cout << "Case3-2 : count_events = "<< count_events<< "\n";

	      unsigned int nr = 0;
	      unsigned int nc = 0;
	      for(unsigned int p=0; p<nPops; p++)
		{
		  nr += possiblePaths.at(trID).at(count_events).at(p).size();
		  nc += possiblePaths.at(trID).at(count_events).at(nPops+p).size();
		}
	      
	      V.resize(nr, Qeigen.at(trID).at(count_events).at(1).rows());
	      V_inv.resize(Qeigen.at(trID).at(count_events).at(1).rows(), nc);
	      unsigned int count_r = 0;
	      unsigned int count_c = 0;
	      for(unsigned int p=0; p<nPops; p++)
		{
		  for(unsigned int k=0; k<possiblePaths.at(trID).at(count_events).at(p).size(); k++)
		    {
		      unsigned int id_row = possiblePaths.at(trID).at(count_events).at(p).at(k);
		      V.row(count_r) =Qeigen.at(trID).at(count_events).at(1).row(id_row);
		      count_r++;
		    }
		  for(unsigned int k=0; k<possiblePaths.at(trID).at(count_events).at(nPops+p).size(); k++)
		    {
		      unsigned int id_col = possiblePaths.at(trID).at(count_events).at(nPops+p).at(k);
		      unsigned int nCoalEvents = nPossibleCoalEvents.at(trID).at(count_events).at(p).at(k);
		      //std::cout << "nPossibleCoalEvents = " << nCoalEvents <<"\n";
		      V_inv.col(count_c) =Qeigen.at(trID).at(count_events).at(2).col(id_col) * nCoalEvents * 2/popSize(0,p);
		      count_c++;
		    }
		}
	      
	      //--- Get exp(D*t) ---//
	      Eigen::MatrixXcd expD = Qeigen.at(trID).at(count_events).at(0) *waitingTime;
	      for(unsigned int di=0; di<expD.rows(); di++)
		{
		  complex<double> val = expD(di,di);
		  expD(di,di) = exp(val);
		}

	      
	      if(probMat.cols() != V.rows() || V.cols() != expD.rows() || expD.cols() != V_inv.rows())
		{
		  std::cout << "\n ***** Error in Chain::compute_conditionalProb_MPI *****\n";	
		  std::cout << "\n ***** Case 3-2: different dimensions *****\n";
		  std::cout << "\n ***** id_sample = "<< id_sample
			    << ", id_locus = "<< id_locus <<", crrProcID = "<< crrProcID << " and poptree is ";
		  poptree->print_allPopTree();
		}
	      else
		{
		  probMat *= V*expD*V_inv;
		}
	      
	    }
	  count_events++; // counting the number of events that happened before the splitting time.
	}
      ++iter;
    } // End of while loop
  
  
  if(ancestralPop ==0)
    {

      //std::cout << "Case4 : count_events = "<< count_events<< "\n";

      complex<double> p_tmp;
      if(probMat.rows() != 1 || probMat.cols() != 2)
	{
	  std::cout << "\n ***** Error in Chain::compute_conditionalProb_MPI *****\n";	
	  std::cout << "\n ***** Case 4: different dimensions *****\n";
	  std::cout << "\n ***** id_sample = "<< id_sample
		    << ", id_locus = "<< id_locus <<", crrProcID = "<< crrProcID << " and poptree is ";
	  poptree->print_allPopTree();
	  std::cout << "The probMat should be 1x2 matrix, but ";
	  std::cout << "probMat = " << probMat <<"\n";
	}
      else
	{
	  p_tmp = probMat(0,0)+probMat(0,1);
	}
      // prob *= p_tmp.real();
      logProb = log(p_tmp.real());
	  // REMOVE
	  //if(id_sample==98 && id_locus ==9)
	  //  std::cout << "prob = " << prob << " coal. before split on procID "<< crrProcID<< "(id_sample ="<< id_sample <<", id_locus="<<id_locus<<")\n";
	  // END of REMOVE
    }
  
  // std::cout << "On process " << crrProcID <<" and prob(return) = "<< prob <<"\n";
  // std::cout << "log(conditional Prob) = " << logProb << "\n";
  

  // std::cout << "Thus, logProb = " << logProb << "\n";
  
  return logProb;
}
*/


// YC 3/8/2015
// 'crrProcID' is required for debugging only. This function is called on each process.
double Chain::compute_logConditionalProb_old(unsigned int id_sample, unsigned int id_locus, popTree* poptree, unsigned int crrProcID)
{
  std::chrono::high_resolution_clock::time_point start_t, end_t;
  std::chrono::high_resolution_clock::time_point start_t_case, end_t_case;
  start_t= std::chrono::high_resolution_clock::now();

  // REMOVE
  /*
  std::cout << "In  Chain::compute_logConditionalProb_MPI()\n";  
  poptree->print_allPopTree();
  std::cout << "id_sample = " << id_sample << " id_locus = " << id_locus <<"\n";
  std::cout << "coaltime  = ";
  for(std::list<double>::iterator ct_iter=coalTimes.at(id_sample).at(id_locus).begin(); ct_iter != coalTimes.at(id_sample).at(id_locus).end(); ct_iter++)
    {
      std::cout <<  *ct_iter <<" ";
    }
  std::cout <<"\n";
  */


  unsigned int nPops = poptree->size();

  // std::cout << "nPops = " << nPops <<"\n";
  
  // FIXME YC 5/9/2014
  // It works for up to 2 population
  // Note: It might be better to have a migration rate matrix and a vector of population sizes
  // as members of class 'popTree'.
  Eigen::MatrixXd popSize(1,nPops);
  popSize.setZero();
  for(unsigned int p1=0; p1< nPops; p1++)
    {
      popSize(0,p1) = poptree->find_popSize(p1+1);
    }
  double splittingTime = poptree->get_age(); // Sorts the elements in ascending order


  double logProb = 0;

  unsigned int trID = treeIDs.at(id_sample).at(id_locus);
  std::list<double> eventT = coalTimes.at(id_sample).at(id_locus);


  // eventT.push_back(splittingTime);
  eventT.sort();

  
  std::list<double>::iterator iter = eventT.begin();
  
  unsigned int nGeneCopies = coalTimes.at(id_sample).at(id_locus).size()+1; //nKinds_lineages.at(trID).size()+1; //list_trees.at(trID)->getSize();
  unsigned int ancestralPop = 0;
  Eigen::MatrixXcd probMat;
  unsigned int count_events=0;
  iter = eventT.begin();
  while(count_events <nGeneCopies-1 && iter!=eventT.end())
    {      
      if(*iter <=splittingTime)
	{	  
	  //--- Case 3: coalescent events in sampling populations ---//


	  // Added by YC 6/11/2014
	  std::list<double>::iterator iter_prev = iter;
	  
	  double waitingTime = 0.0;
	  if(count_events==0)
	    {
	      waitingTime = *iter;
	    }
	  else
	    {		      
	      iter_prev--;		 
	      waitingTime = *iter - *iter_prev;
	    }
	  
	  Eigen::MatrixXcd V;
	  Eigen::MatrixXcd V_inv;
	  Eigen::MatrixXcd probMat_each;
	  
	  if(count_events==0) // the first coalescent event
	    {
	      //std::cout << "Case3-1 : count_events = "<< count_events<< "\n";
	      count_case3_1++;	      
	      start_t_case= std::chrono::high_resolution_clock::now();

	      //--- Get D*t ---//
	      Eigen::ArrayXcd eigenvalues = Qeigen.at(trID).at(count_events).at(0)*waitingTime;

	      unsigned int nc = 0;
	      for(unsigned int p=0; p<nPops; p++)
		nc += possiblePaths.at(trID).at(count_events).at(p+1).size();
	      
	      //--- Computing V*exp(D*t) ---//
	      V = Qeigen.at(trID).at(count_events).at(1).row(possiblePaths.at(trID).at(count_events).at(0).at(0)).cwiseProduct(eigenvalues.exp().matrix().transpose());
	      //V_inv.resize(Qeigen.at(trID).at(count_events).at(1).rows(), nc);
	      probMat_each.resize(1,nc);
	      unsigned int count = 0;
	      for(unsigned int p=0; p<nPops; p++)
		{
		  for(unsigned int k=0; k<possiblePaths.at(trID).at(count_events).at(p+1).size(); k++)
		    {
		      unsigned int id_col = possiblePaths.at(trID).at(count_events).at(p+1).at(k);
		      unsigned int nCoalEvents = nPossibleCoalEvents.at(trID).at(count_events).at(p).at(k);
		      V_inv =Qeigen.at(trID).at(count_events).at(2).col(id_col)* nCoalEvents *2/popSize(0,p);
		      //V_inv.col(count) =Qeigen.at(trID).at(count_events).at(2).col(id_col)* nCoalEvents *2/popSize(0,p);
		      probMat_each.col(count) = V*V_inv;
		      count++;
		    }
		}
	      
	      //--- Get exp(D*t) ---//
	      //Eigen::ArrayXcd eigenvalues = Qeigen.at(trID).at(count_events).at(0)*waitingTime;
	      //V=V.cwiseProduct(eigenvalues.exp().matrix().transpose());
	      // V = (V.array()*eigenvalues.exp()).matrix();
	      //Eigen::MatrixXd expD = eigenvalues.exp().matrix();
	      /*
	      for(unsigned int di=0; di<expD.rows(); di++)
		{
		  complex<double> val = expD(di,di);
		  expD(di,di) = exp(val*waitingTime);
		}
	      */
	      
	      // if(V.cols() != expD.rows() || expD.cols() != V_inv.rows())
	      /*
	      if(V.cols() != V_inv.rows())
		{
		  std::cout << "\n ***** Error in Chain::compute_conditionalProb_MPI *****\n";	
		  std::cout << "\n ***** Case 3-1: different dimensions *****\n";		  
		  std::cout << "\n ***** id_sample = "<< id_sample
			    << ", id_locus = "<< id_locus <<", crrProcID = "<< crrProcID << " and poptree is ";
		  poptree->print_allPopTree();
		}
	      else
		{
		  // probMat = V*expD*V_inv;
		  //probMat = V*V_inv;
		}
	      */
	      probMat = probMat_each;
	      
	      end_t_case= std::chrono::high_resolution_clock::now();	      
	      eachComputingTime_condiProb_case3_1 += std::chrono::duration_cast<std::chrono::microseconds>(end_t_case - start_t_case);
	    }
	  else // non-first coalescent event
	    {
	      //std::cout << "Case3-2 : count_events = "<< count_events<< "\n";
	      count_case3_2++;
	      start_t_case= std::chrono::high_resolution_clock::now();

	      //--- Get exp(D*t) ---//
	      Eigen::ArrayXcd eigenvalues = Qeigen.at(trID).at(count_events).at(0)*waitingTime;
	      //Eigen::MatrixXd expD = (eigenvalues.exp().matrix()).asDiagonal(); 
	      //Eigen::MatrixXd expD = ((Qeigen.at(trID).at(count_events).at(0).real() *waitingTime).exp()).asDiagonal();
	      /*
	      for(unsigned int di=0; di<expD.rows(); di++)
		{
		  complex<double> val = expD(di,di);
		  expD(di,di) = exp(val);
		}
	      */

	      unsigned int nr = 0;
	      unsigned int nc = 0;
	      for(unsigned int p=0; p<nPops; p++)
		{
		  nr += possiblePaths.at(trID).at(count_events).at(p).size();
		  nc += possiblePaths.at(trID).at(count_events).at(nPops+p).size();
		}
	      
	      V.resize(nr, Qeigen.at(trID).at(count_events).at(1).rows());
	      V_inv.resize(Qeigen.at(trID).at(count_events).at(1).rows(), nc);
	      unsigned int count_r = 0;
	      unsigned int count_c = 0;
	      for(unsigned int p=0; p<nPops; p++)
		{
		  for(unsigned int k=0; k<possiblePaths.at(trID).at(count_events).at(p).size(); k++)
		    {
		      unsigned int id_row = possiblePaths.at(trID).at(count_events).at(p).at(k);
		      V.row(count_r) =Qeigen.at(trID).at(count_events).at(1).row(id_row).cwiseProduct(eigenvalues.exp().matrix().transpose());
		      count_r++;
		    }
		  for(unsigned int k=0; k<possiblePaths.at(trID).at(count_events).at(nPops+p).size(); k++)
		    {
		      unsigned int id_col = possiblePaths.at(trID).at(count_events).at(nPops+p).at(k);
		      unsigned int nCoalEvents = nPossibleCoalEvents.at(trID).at(count_events).at(p).at(k);
		      //std::cout << "nPossibleCoalEvents = " << nCoalEvents <<"\n";
		      V_inv.col(count_c) =Qeigen.at(trID).at(count_events).at(2).col(id_col) * nCoalEvents * 2/popSize(0,p);
		      count_c++;
		    }
		}
	      
	      
	      //if(probMat.cols() != V.rows() || V.cols() != expD.rows() || expD.cols() != V_inv.rows())
	      if(probMat.cols() != V.rows() || V.cols() !=  V_inv.rows())
		{
		  std::cout << "\n ***** Error in Chain::compute_conditionalProb_MPI *****\n";	
		  std::cout << "\n ***** Case 3-2: different dimensions *****\n";
		  std::cout << "\n ***** id_sample = "<< id_sample
			    << ", id_locus = "<< id_locus <<", crrProcID = "<< crrProcID << " and poptree is ";
		  poptree->print_allPopTree();
		}
	      else
		{
		  probMat *= V*V_inv;
		  // probMat *= V*expD*V_inv;
		}
	      
	      end_t_case= std::chrono::high_resolution_clock::now();	      
	      eachComputingTime_condiProb_case3_2 += std::chrono::duration_cast<std::chrono::microseconds>(end_t_case - start_t_case);
	    }

	  //std::cout << "Case 3: *iter = " << *iter << " splittingtime = " << splittingTime
	  //	    << " waitingTime= " << waitingTime << " and probMat=" << probMat <<"\n";

	  count_events++; // counting the number of events that happened before the splitting time.
	  	  
	  // Added by YC 3/8/2015
	  // During optimization (differential evolution), proposed splitting times
	  // can be same as coalescent times. 
	  // If the splitting time and a coalescent time are the same, then
	  // we assume the coalescent event happened in a sampling population 
	  // before the splitting event backward in time.
	  // If the splitting time is same as a coalescent time smaller than
	  // the TMRCA (tree height) (i.e., count_events <nGeneCopies-1), then
	  // we have a flag of 'ancestralPop = 1'. 
	  // If ancestralPop =1, then we do NOT need to compute the probability of
	  // no coalescent events between the last observed coalescent event time in 
	  // in sampling populations and the splitting time and 'Case 1' is called. 
	  if(*iter == splittingTime  && count_events <nGeneCopies-1)
	    {
	      ancestralPop = 1;
	    }
	  
	  ++iter;

	}
      else if(*iter > splittingTime && ancestralPop ==0)
	{
	  //--- Case 2: No coalescent event between the last coalescent ---//
	  //--- event in the sampling populations and the splitting time --//
	  

	  std::list<double>::iterator iter_prev = iter;
	  
	  double waitingTime = 0.0;
	  if(count_events==0)
	    {
	      waitingTime = splittingTime;
	    }
	  else
	    {		      
	      iter_prev--; 
	      waitingTime = splittingTime - *iter_prev;
	    }	  
	  
	  Eigen::MatrixXcd V;
	  Eigen::MatrixXcd V_inv;
	  if(count_events==0)
	    {
	      count_case2_1++;
	      start_t_case= std::chrono::high_resolution_clock::now();

	      V = Qeigen.at(trID).at(count_events).at(1).row(possiblePaths.at(trID).at(count_events).at(0).at(0));
	      V_inv = Qeigen.at(trID).at(count_events).at(2);
	      	      
	      //--- Get exp(D*t) ---//
	      Eigen::ArrayXcd eigenvalues = Qeigen.at(trID).at(count_events).at(0)*waitingTime;
	      V=V.cwiseProduct(eigenvalues.exp().matrix().transpose());
	      //Eigen::MatrixXd expD = (eigenvalues.exp().matrix()).asDiagonal(); 
	      // Eigen::MatrixXd expD = ((Qeigen.at(trID).at(count_events).at(0).real()*waitingTime).exp()).asDiagonal();
	      /*
	      for(unsigned int di=0; di<expD.rows(); di++)
		{
		  complex<double> val = expD(di,di);
		  expD(di,di) = exp(val *waitingTime);
		}
	      */
	      // if(V.cols() != expD.rows() || expD.cols() != V_inv.rows())
	      if(V.cols() != V_inv.rows())
		{
		  std::cout << "\n ***** Error in Chain::compute_conditionalProb_MPI *****\n";	
		  std::cout << "\n ***** Case 2-1: different dimensions *****\n";
		  std::cout << "trID = " << trID <<"\n";
		  std::cout << "possiblePaths.at(trID).at(count_events).at(0).at(0) ="<<possiblePaths.at(trID).at(count_events).at(0).at(0)<<"\n";
		  std::cout << "V = " << V <<"\n"
			    <<"V_inv = " << V_inv<<"\n"
			    <<"eigenvalues = " <<eigenvalues <<"\n";
		}
	      else
		{
		  probMat = V*V_inv;
		  // probMat = V*expD*V_inv;
		}

	      end_t_case= std::chrono::high_resolution_clock::now();	      
	      eachComputingTime_condiProb_case2_1 += std::chrono::duration_cast<std::chrono::microseconds>(end_t_case - start_t_case);
	    }
	  else // if(count_events != 0 && waitingTime > 0.0)
	    {
	      //std::cout << "Case2-2 : count_events = "<< count_events<< "\n";
	      count_case2_2++;
	      start_t_case= std::chrono::high_resolution_clock::now();

	      V_inv = Qeigen.at(trID).at(count_events).at(2);
	      
	      //--- Get exp(D*t) ---//
	      Eigen::ArrayXcd eigenvalues = Qeigen.at(trID).at(count_events).at(0)*waitingTime;
	      //Eigen::MatrixXd expD = (eigenvalues.exp().matrix()).asDiagonal(); 
	      // Eigen::MatrixXd expD = ((Qeigen.at(trID).at(count_events).at(0).real()*waitingTime).exp()).asDiagonal();
	      /*
	      for(unsigned int di=0; di<expD.rows(); di++)
		{
		  complex<double> val = expD(di,di);
		  expD(di,di) = exp(val);
		}
	      */
	      unsigned int nr = 0;
	      for(unsigned int p=0; p<nPops; p++)
		{
		  nr += possiblePaths.at(trID).at(count_events).at(p).size();
		}
	      V.resize(nr, Qeigen.at(trID).at(count_events).at(1).rows());
	      unsigned int count_r = 0;
	      for(unsigned int p=0; p<nPops; p++)
		{
		  for(unsigned int k=0; k<possiblePaths.at(trID).at(count_events).at(p).size(); k++)
		    {
		      unsigned int id_row = possiblePaths.at(trID).at(count_events).at(p).at(k);
		      //V.row(count_r) =Qeigen.at(trID).at(count_events).at(1).row(id_row);
		      V.row(count_r) =Qeigen.at(trID).at(count_events).at(1).row(id_row).cwiseProduct(eigenvalues.exp().matrix().transpose());
		      count_r++;
		    }
		}
	      
	      //if(probMat.cols() != V.rows() || V.cols() != expD.rows() || expD.cols() != V_inv.rows())
	      if(probMat.cols() != V.rows() || V.cols() != V_inv.rows())
		{
		  std::cout << "\n ***** Error in Chain::compute_conditionalProb_MPI *****\n";	
		  std::cout << "\n ***** Case 2-2: different dimensions *****\n";
		}
	      else
		{
		  probMat *= V*V_inv;
		  //probMat *= V*expD*V_inv;
		}
	      end_t_case= std::chrono::high_resolution_clock::now();	      
	      eachComputingTime_condiProb_case2_2 += std::chrono::duration_cast<std::chrono::microseconds>(end_t_case - start_t_case);
	    }
	  
	  double prob_noAbsorbing = 0.0;
	  for(unsigned int id_No_AbsorbingState = 0; id_No_AbsorbingState < probMat.cols()-1; id_No_AbsorbingState++)
	    prob_noAbsorbing += probMat(0,id_No_AbsorbingState).real();

	  logProb += log(prob_noAbsorbing);
	  
	  // The flag that 'Case 1' should be called.
	  ancestralPop =1;
	}
      else if(*iter > splittingTime && ancestralPop == 1)
	{	  
	  //--- Case 1: coalescent in the ancestral population. ---//
	  //  Some or all coalescent events happened in the ancestral population.
	  
	  // std::cout << "case 1\n";
	  count_case1++;
	  start_t_case= std::chrono::high_resolution_clock::now();

	  double popSize = poptree->get_popSize();
	  unsigned int nLineages = nGeneCopies - count_events; // remaining number of lineages
	  double expterm = 0.0;
	  for(unsigned int l=nLineages; l>=2; l--, ++iter)
	    {
	      std::list<double>::iterator iter_prev = iter;
	      double waitingTime =0.0;
	      if(l == nLineages) // the first coalescent event in the ancestral population
		waitingTime = *iter - splittingTime;
	      else
		{
		  iter_prev--;
		  waitingTime = *iter - *iter_prev;
		}

	      expterm += l*(l-1)* waitingTime / popSize;
	    }
	  logProb += (nLineages-1)*log(2/popSize) - expterm;

	  count_events = nGeneCopies-1;
	  ancestralPop = 1;

	  end_t_case= std::chrono::high_resolution_clock::now();	      
	  eachComputingTime_condiProb_case1 += std::chrono::duration_cast<std::chrono::microseconds>(end_t_case - start_t_case);
	}
    } // End of while loop
  
  
  if(ancestralPop ==0)
    {

      //std::cout << "Case4 : count_events = "<< count_events<< "\n";

      complex<double> p_tmp;
      if(probMat.rows() != 1 || probMat.cols() != 2)
	{
	  std::cout << "\n ***** Error in Chain::compute_conditionalProb_MPI *****\n";	
	  std::cout << "\n ***** Case 4: different dimensions *****\n";
	  std::cout << "\n ***** id_sample = "<< id_sample
		    << ", id_locus = "<< id_locus <<", crrProcID = "<< crrProcID << " and poptree is ";
	  poptree->print_allPopTree();
	  std::cout << "The probMat should be 1x2 matrix, but ";
	  std::cout << "probMat = " << probMat <<"\n";
	}
      else
	{
	  p_tmp = probMat(0,0)+probMat(0,1);
	}
      double p_real = p_tmp.real();
      if(p_real <0)
	{
	  p_real *= -1;
	  if(p_real < pow(10,-10))
	    p_real = 0;
	  else
	    {
	      std::cout << "\n***Error in Chain::compute_logConditionalProb()\n";
	      std::cout << "p_tmp.real() = " << p_tmp.real() <<"\n";
	    }
	}

      logProb = log(p_real);
	  // REMOVE
	  //if(id_sample==98 && id_locus ==9)
	  //  std::cout << "prob = " << prob << " coal. before split on procID "<< crrProcID<< "(id_sample ="<< id_sample <<", id_locus="<<id_locus<<")\n";
	  // END of REMOVE
    }
  
  // std::cout << "On process " << crrProcID <<" and prob(return) = "<< prob <<"\n";
  // std::cout << "log(conditional Prob) = " << logProb << "\n";
  /*
  if(prob <0 || prob ==0)
    {
      std::cout << "Error in Chain::compute_conditionalProb_MPI(). id_sample = "<< id_sample
		<< ", id_locus = "<< id_locus <<", crrProcID = "<< crrProcID << " and poptree is ";
      poptree->print_allPopTree();
      std::cout << "prob = " << prob <<"\n";
    }
  */

  // std::cout << "Thus, logProb = " << logProb << "\n";
  end_t= std::chrono::high_resolution_clock::now();
  eachComputingTime_condiProb += std::chrono::duration_cast<std::chrono::microseconds>(end_t - start_t);
  count_condiProbFunctionCalls++;
  
  return logProb;
}



double Chain::compute_conditionalProb_MPI(unsigned int id_sample, unsigned int id_locus, popTree* poptree, unsigned int crrProcID)
{
  unsigned int nPops = poptree->size();
  
  // FIXME YC 5/9/2014
  // It works for up to 2 population
  // Note: It might be better to have a migration rate matrix and a vector of population sizes
  // as members of class 'popTree'.
  Eigen::MatrixXd popSize(1,nPops);
  popSize.setZero();
  for(unsigned int p1=0; p1< nPops; p1++)
    {
      popSize(0,p1) = poptree->find_popSize(p1+1);
    }
  double splittingTime = poptree->get_age();

  double scaler = 100000;
  double prob = 1.0*scaler;
  unsigned int trID = treeIDs.at(id_sample).at(id_locus);
  std::list<double> eventT = coalTimes.at(id_sample).at(id_locus);
  eventT.push_back(splittingTime);
  eventT.sort();
  
  
  std::list<double>::iterator iter = eventT.begin();
  
  unsigned int nGeneCopies = coalTimes.at(id_sample).at(id_locus).size()+1; //nKinds_lineages.at(trID).size()+1; //list_trees.at(trID)->getSize();
  unsigned int ancestralPop = 0;
  Eigen::MatrixXcd probMat;
  unsigned int count_events=0;
  iter = eventT.begin();
  while(count_events <nGeneCopies-1 && iter!=eventT.end() && ancestralPop == 0 )
    {
      
      if(*iter > splittingTime)
	{
	  // case 1: coalescent in the ancestral population
	  double popSize = poptree->get_popSize();
	  unsigned int nLineages = nGeneCopies - count_events; // remaining number of lineages
	  double expterm = 0.0;
	  for(unsigned int l=nLineages; l>=2; l--, ++iter)
	    {
	      // Added by YC 6/11/2014
	      std::list<double>::iterator iter_prev = iter;
	      iter_prev--;
	      double waitingTime = *iter - *iter_prev;
	      
	      // Modified by YC 6/10/2014
	      expterm += l*(l-1)* waitingTime / popSize;
	      unsigned int nCoalEvents = nPossibleCoalEvents.at(trID).at(nGeneCopies-l).at(nPops).at(0);
	      prob  *= nCoalEvents * 2/popSize;
	    }
	  double prob_beforeExp = prob;
	  prob *= exp(-1*expterm);
	  // REMOVE
	  if(prob == 0)
	    {
	      std::cout << "prob = " << prob << " coal. in the ancetral pop on procID "<< crrProcID<< "(id_sample ="<< id_sample <<", id_locus="<<id_locus<<")"<<"expterm = "<<expterm <<"prob_beforeExp="<< prob_beforeExp<<"\n";
	    }
	  //if(id_sample==98 && id_locus ==9)
	  //  std::cout << "prob = " << prob << " coal. in the ancetral pop on procID "<< crrProcID<< "(id_sample ="<< id_sample <<", id_locus="<<id_locus<<")\n";
	  // END of REMOVE
	  ancestralPop = 1;
	}
      else if(*iter == splittingTime) 
	{
	  // -----------------------------------------------------------------//
	  // If all the coalescent events happened before the splitting time (backward in time),
	  // this part won't be called. 
	  // This part is called if and only if some or all coalescent events
	  // happened after the splitting time or at the same time of the split.
	  // Therefore, in the next round, the case of (*iter > splittingTime) 
	  // should be called and "ancestralPop = 1" will be set in that case.
	  // ----------------------------------------------------------------//

	  //std::cout << "Case2 : count_events = "<< count_events<< "\n";

	  // Added by YC 6/11/2014
	  std::list<double>::iterator iter_prev = iter;
	  
	  double waitingTime = 0.0;
	  if(count_events==0)
	    {
	      waitingTime = *iter;
	    }
	  else
	    {		      
	      iter_prev--; 
	      waitingTime = *iter - *iter_prev;
	    }
	  
	  
	  if(waitingTime > 0.0)
	    {
	      Eigen::MatrixXcd V;
	      Eigen::MatrixXcd V_inv;
	      if(count_events==0)
		{
		  // std::cout << "id_sample = " << id_sample <<" id_locus = " << id_locus <<"\n";

		  V = Qeigen.at(trID).at(count_events).at(1).row(possiblePaths.at(trID).at(count_events).at(0).at(0));
		  V_inv = Qeigen.at(trID).at(count_events).at(2);
		  
		  //--- Get exp(D*t) ---//
		  Eigen::MatrixXcd expD = Qeigen.at(trID).at(count_events).at(0) *waitingTime;
		  for(unsigned int di=0; di<expD.rows(); di++)
		    {
		      complex<double> val = expD(di,di);
		      expD(di,di) = exp(val);
		    }
		  if(V.cols() != expD.rows() || expD.cols() != V_inv.rows())
		    {
		      std::cout << "\n ***** Error in Chain::compute_conditionalProb_MPI *****\n";	
		      std::cout << "\n ***** Case 2-1: different dimensions *****\n";
		      std::cout << "trID = " << trID <<"\n";
		      std::cout << "possiblePaths.at(trID).at(count_events).at(0).at(0) ="<<possiblePaths.at(trID).at(count_events).at(0).at(0)<<"\n";
		      std::cout << "V = " << V <<"\n"
				<<"V_inv = " << V_inv<<"\n"
				<<"expD = " <<expD <<"\n";
		    }
		  else
		    {
		      probMat = V*expD*V_inv;
		    }
		  std::cout << "probMat = " << probMat <<"\n";
		}
	      else if(waitingTime > 0.0) // if(count_events != 0 && waitingTime > 0.0)
		{
		  //std::cout << "Case2-2 : count_events = "<< count_events<< "\n";
		  
		  V_inv = Qeigen.at(trID).at(count_events).at(2);
		  
		  unsigned int nr = 0;
		  for(unsigned int p=0; p<nPops; p++)
		    {
		      nr += possiblePaths.at(trID).at(count_events).at(p).size();
		    }
		  V.resize(nr, Qeigen.at(trID).at(count_events).at(1).rows());
		  unsigned int count_r = 0;
		  for(unsigned int p=0; p<nPops; p++)
		    {
		      for(unsigned int k=0; k<possiblePaths.at(trID).at(count_events).at(p).size(); k++)
			{
			  unsigned int id_row = possiblePaths.at(trID).at(count_events).at(p).at(k);
			  V.row(count_r) =Qeigen.at(trID).at(count_events).at(1).row(id_row);
			  count_r++;
			}
		    }
		  
		  //--- Get exp(D*t) ---//
		  Eigen::MatrixXcd expD = Qeigen.at(trID).at(count_events).at(0) *waitingTime;
		  for(unsigned int di=0; di<expD.rows(); di++)
		    {
		      complex<double> val = expD(di,di);
		      expD(di,di) = exp(val);
		    }
		  try
		    {
		      if(probMat.cols() != V.rows() || V.cols() != expD.rows() || expD.cols() != V_inv.rows())
			{
			  std::cout << "\n ***** Error in Chain::compute_conditionalProb_MPI *****\n";	
			  std::cout << "\n ***** Case 2-2: different dimensions *****\n";
			}
		      else
			{
			  probMat *= V*expD*V_inv;
			}
		    }
		  catch(std::exception &e)
		    {
		      std::cout << "\n ***** Error in Chain::compute_conditionalProb_MPI *****\n";		  
		    }
		}

	      double prob_noAbsorbing = 0.0;
	      for(unsigned int id_No_AbsorbingState = 0; id_No_AbsorbingState < probMat.cols()-1; id_No_AbsorbingState++)
		prob_noAbsorbing += probMat(0,id_No_AbsorbingState).real();
	      prob *= prob_noAbsorbing;

	      //std::cout << "prob = " << prob <<"\n";
	    }
	  else if(waitingTime == 0.0)
	    {
	      // std::cout << "case 2-3 : count_events = "<< count_events<< "\n";

	      // ---------------------------------------------------------//
	      // YC 7/14/2014
	      // When a coalescent event and splitting happened at the same time,
	      // we assume the coalescent event happened in the ancesteral population.
	      // ---------------------------------------------------------//

	      prob *= 2 /poptree->get_popSize();
	      ancestralPop =1;
	      count_events++;
	    }

	
	}
      else // if (*iter < splittingTime)
	{	  
	  // Added by YC 6/11/2014
	  std::list<double>::iterator iter_prev = iter;
	  
	  double waitingTime = 0.0;
	  if(count_events==0)
	    {
	      waitingTime = *iter;
	    }
	  else
	    {		      
	      iter_prev--;		 
	      waitingTime = *iter - *iter_prev;
	    }

	  Eigen::MatrixXcd V;
	  Eigen::MatrixXcd V_inv;
	  //Eigen::MatrixXcd V_inv_test;
	  //Eigen::MatrixXcd Coal_test;
	  
	  if(count_events==0)
	    {
	      //std::cout << "Case3-1 : count_events = "<< count_events<< "\n";

	      unsigned int nc = 0;
	      for(unsigned int p=0; p<nPops; p++)
		nc += possiblePaths.at(trID).at(count_events).at(p+1).size();
	      
	      V = Qeigen.at(trID).at(count_events).at(1).row(possiblePaths.at(trID).at(count_events).at(0).at(0));
	      V_inv.resize(Qeigen.at(trID).at(count_events).at(1).rows(), nc);
	      unsigned int count = 0;
	      for(unsigned int p=0; p<nPops; p++)
		{
		  for(unsigned int k=0; k<possiblePaths.at(trID).at(count_events).at(p+1).size(); k++)
		    {
		      unsigned int id_col = possiblePaths.at(trID).at(count_events).at(p+1).at(k);
		      unsigned int nCoalEvents = nPossibleCoalEvents.at(trID).at(count_events).at(p).at(k);
		      //std::cout << "nPossibleCoalEvents = " << nCoalEvents <<"\n";
		      V_inv.col(count) =Qeigen.at(trID).at(count_events).at(2).col(id_col)* nCoalEvents *2/popSize(0,p);
		      count++;
		    }
		}
	      
	      //--- Get exp(D*t) ---//
	      Eigen::MatrixXcd expD = Qeigen.at(trID).at(count_events).at(0) *waitingTime;
	      for(unsigned int di=0; di<expD.rows(); di++)
		{
		  complex<double> val = expD(di,di);
		  expD(di,di) = exp(val);
		}
	      
	      if(V.cols() != expD.rows() || expD.cols() != V_inv.rows())
		{
		  std::cout << "\n ***** Error in Chain::compute_conditionalProb_MPI *****\n";	
		  std::cout << "\n ***** Case 3-1: different dimensions *****\n";		  
		  std::cout << "\n ***** id_sample = "<< id_sample
			    << ", id_locus = "<< id_locus <<", crrProcID = "<< crrProcID << " and poptree is ";
		  poptree->print_allPopTree();
		}
	      else
		{
		  probMat = V*expD*V_inv;
		}
	      
		      
	    }
	  else
	    {
	      //std::cout << "Case3-2 : count_events = "<< count_events<< "\n";

	      unsigned int nr = 0;
	      unsigned int nc = 0;
	      for(unsigned int p=0; p<nPops; p++)
		{
		  nr += possiblePaths.at(trID).at(count_events).at(p).size();
		  nc += possiblePaths.at(trID).at(count_events).at(nPops+p).size();
		}
	      
	      V.resize(nr, Qeigen.at(trID).at(count_events).at(1).rows());
	      V_inv.resize(Qeigen.at(trID).at(count_events).at(1).rows(), nc);
	      unsigned int count_r = 0;
	      unsigned int count_c = 0;
	      for(unsigned int p=0; p<nPops; p++)
		{
		  for(unsigned int k=0; k<possiblePaths.at(trID).at(count_events).at(p).size(); k++)
		    {
		      unsigned int id_row = possiblePaths.at(trID).at(count_events).at(p).at(k);
		      V.row(count_r) =Qeigen.at(trID).at(count_events).at(1).row(id_row);
		      count_r++;
		    }
		  for(unsigned int k=0; k<possiblePaths.at(trID).at(count_events).at(nPops+p).size(); k++)
		    {
		      unsigned int id_col = possiblePaths.at(trID).at(count_events).at(nPops+p).at(k);
		      unsigned int nCoalEvents = nPossibleCoalEvents.at(trID).at(count_events).at(p).at(k);
		      //std::cout << "nPossibleCoalEvents = " << nCoalEvents <<"\n";
		      V_inv.col(count_c) =Qeigen.at(trID).at(count_events).at(2).col(id_col) * nCoalEvents * 2/popSize(0,p);
		      count_c++;
		    }
		}
	      
	      //--- Get exp(D*t) ---//
	      Eigen::MatrixXcd expD = Qeigen.at(trID).at(count_events).at(0) *waitingTime;
	      for(unsigned int di=0; di<expD.rows(); di++)
		{
		  complex<double> val = expD(di,di);
		  expD(di,di) = exp(val);
		}

	      
	      if(probMat.cols() != V.rows() || V.cols() != expD.rows() || expD.cols() != V_inv.rows())
		{
		  std::cout << "\n ***** Error in Chain::compute_conditionalProb_MPI *****\n";	
		  std::cout << "\n ***** Case 3-2: different dimensions *****\n";
		  std::cout << "\n ***** id_sample = "<< id_sample
			    << ", id_locus = "<< id_locus <<", crrProcID = "<< crrProcID << " and poptree is ";
		  poptree->print_allPopTree();
		}
	      else
		{
		  probMat *= V*expD*V_inv;
		}
	      
	    }
	  count_events++; // counting the number of events that happened before the splitting time.
	}
      ++iter;
    } // End of while loop
  
  
  if(ancestralPop ==0)
    {

      //std::cout << "Case4 : count_events = "<< count_events<< "\n";

      complex<double> p_tmp;
      if(probMat.rows() != 1 || probMat.cols() != 2)
	{
	  std::cout << "\n ***** Error in Chain::compute_conditionalProb_MPI *****\n";	
	  std::cout << "\n ***** Case 4: different dimensions *****\n";
	  std::cout << "\n ***** id_sample = "<< id_sample
		    << ", id_locus = "<< id_locus <<", crrProcID = "<< crrProcID << " and poptree is ";
	  poptree->print_allPopTree();
	}
      else
	{
	  p_tmp = probMat(0,0)+probMat(0,1);
	}
      prob *= p_tmp.real();
	  // REMOVE
	  if(id_sample==98 && id_locus ==9)
	    std::cout << "prob = " << prob << " coal. before split on procID "<< crrProcID<< "(id_sample ="<< id_sample <<", id_locus="<<id_locus<<")\n";
	  // END of REMOVE
    }
  
  // std::cout << "On process " << crrProcID <<" and prob(return) = "<< prob <<"\n";
  //std::cout << "conditional Prob = " << prob << "\n";
  if(prob <0 || prob ==0)
    {
      std::cout << "Error in Chain::compute_conditionalProb_MPI(). id_sample = "<< id_sample
		<< ", id_locus = "<< id_locus <<", crrProcID = "<< crrProcID << " and poptree is ";
      poptree->print_allPopTree();
      std::cout << "prob = " << prob <<"\n";
    }
  
  return prob;
}



double Chain::compute_jointPosteriorDensity(popTree* poptree, IM im)
{
  // REMOVE
  //  std::cout << "Computing Eigen values and Eigen vectors..\n";
  compute_eigenValuesVectors_rateMat(poptree);
  //std::cout << "DONE with computing Eigen values and Eigen vectors\n";

	double posteriorDensity = 0.0;

	unsigned int nSample = im.get_nSampledTrees();
	unsigned int nloci = im.get_nLoci();

	for(unsigned int i=0; i< nSample ; i++)
	{
		double eachterm = 1.0;
		for(unsigned int j=0; j<nloci; j++)
		{
			double eachProb = compute_conditionalProb(i,j,poptree) /numTrees.at(treeIDs.at(i).at(j));
			eachterm *= eachProb;

		}
		posteriorDensity += eachterm/ exp(logPrior_trees.at(i));
	}

	Eigen::Vector3d paraMax = im.get_paraMax();
	double priorPopTree = poptree->computeJointPrior(paraMax);
	posteriorDensity *= priorPopTree/(nSample);

	return posteriorDensity;
}


// old version
double Chain::compute_partialJointPosteriorDensity_MPI_overSample(popTree* poptree, IM im, unsigned int crrProcID, unsigned int nProcs)
{
  
  if(crrProcID==0)
    {
      std::cout << "\n In Chain::compute_partialJointPosteriorDensity_MPI_overSample()\n";
      poptree->print_allPopTree();
    }
  
  // compute_eigenValuesVectors_rateMat_MPI(poptree);
  compute_eigenValuesVectors(poptree);

  double posteriorDensity_log = 0.0;
  double posteriorDensity_partial = 1;
  double posteriorDensity_partial_local = 1;
  double posteriorDensity_log_partial = 0.0;
  double posteriorDensity = 0.0;

  unsigned int nSample = im.get_nSampledTrees();
  unsigned int nloci = im.get_nLoci();

  if(percentile_4upperBoundOfTMRCA != 100)
    {
      // By YC 12/15/2014
      // The upper bound of tmrca for tree included in the computation is
      //      (spliitng time) + 99th percentile of 1st coal. time + that of 2nd coal. time + ...
      
      double partialTerm = 0.0;
      unsigned int nLin = coalTimes.at(0).at(0).size()+1;
      for(unsigned int i=nLin; i>=2; i--)
	{
	  partialTerm += 2/(i*(i-1));
	}
      upperBoundOfTMRCA = poptree->get_age() - poptree->get_popSize()/2 * log(1-percentile_4upperBoundOfTMRCA/100) * partialTerm;

    }

  nSamples_belowTheUpperBoundOfTMRCA =0;
  for(unsigned int i=0; i< nSample ; i++)
    {

      unsigned int largerThanTheUpperBoundOfTMRCA = 0;
      if(percentile_4upperBoundOfTMRCA != 100)
	{
	  unsigned int count = 0;
	  if(n_loci == 0)
	    {
	      n_loci = coalTimes.at(i).size();
	    }
	  while(largerThanTheUpperBoundOfTMRCA==0 && count < n_loci)
	    {
	      std::list<double> ctimes = coalTimes.at(i).at(count);
	      ctimes.sort(std::greater<double>()); // sort in descending order
	      std::list<double>::iterator iter_ctimes = ctimes.begin();
	      if(*iter_ctimes > upperBoundOfTMRCA)
		{
		  largerThanTheUpperBoundOfTMRCA= 1;
		}
	      count++;
	    
	    }
	}

      unsigned int procID2compute = i-nProcs*static_cast<int>(i/nProcs);
      if(procID2compute == crrProcID && largerThanTheUpperBoundOfTMRCA ==0)
	{
	  nSamples_belowTheUpperBoundOfTMRCA++;

	  double eachterm_log = 0.0;
	  
	  for(unsigned int j=0; j<nloci; j++)
	    {
	      //std::cout << "i (upto nSample) = " << i << " j(<= nloci) = "<< j <<"\n";
	      // double eachProb = compute_conditionalProb_MPI(i,j,poptree,crrProcID);
	      double eachLogProb = 0;
	      eachLogProb = compute_logConditionalProb(i,j,poptree,crrProcID);
	      double ntrees = numTrees.at(treeIDs.at(i).at(j)); // the number of coalescent trees that have the give tree topology with the special labels
	      
	      eachterm_log += eachLogProb-log(ntrees);
	      //eachterm_log += log(eachProb)-log(ntrees);
	      //std::cout << "\neachterm = " << eachterm_log << "\n";
	      
	      if(eachterm_log == numeric_limits<double>::infinity() || eachterm_log == -1*numeric_limits<double>::infinity())
		{
		  std::cout << "*** Error *** in Chain::compute_jointPosteriorDensity_MPI()\n";
		  std::cout << "\t on " << i << "th sample and " << j<<"th loci: eachterm_log = " << eachterm_log <<"\n";
		  std::cout << "\t eachLogProb =  " << eachLogProb <<" and log(ntrees) = " << log(ntrees) <<"\n";
		}
	    }
	  //std::cout << "i= " << i <<"\n";
	  //std::cout << "logPrior_trees.at(i) = " << logPrior_trees.at(i) <<"\n";
	  if(posteriorDensity_log_partial == 0)
	    {
	      posteriorDensity_log_partial = eachterm_log -logPrior_trees.at(i);
	    }
	  else if(exp(eachterm_log) == 0 || exp(posteriorDensity_log_partial) == 0)
	    {
	      double logLarge, logSmall;
	      if(posteriorDensity_log_partial < eachterm_log)
		{
		  logLarge = eachterm_log-logPrior_trees.at(i);
		  logSmall = posteriorDensity_log_partial;
		}
	      else
		{
		  logLarge = posteriorDensity_log_partial;
		  logSmall = eachterm_log-logPrior_trees.at(i);
		}
	      posteriorDensity_log_partial = logLarge+log(1+exp(logSmall-logLarge));
	      if(posteriorDensity_log_partial == numeric_limits<double>::infinity())
		std::cout << "case1: posteriorDensity_log_partial = " << posteriorDensity_log_partial << " on procID = " << crrProcID <<"\n";
	    }
	  else
	    {
	      double factor2sum = exp(eachterm_log -logPrior_trees.at(i));
	      posteriorDensity_log_partial = log(exp(posteriorDensity_log_partial)+factor2sum);
	      if(posteriorDensity_log_partial == numeric_limits<double>::infinity())
		std::cout << "case2: posteriorDensity_log_partial = " << posteriorDensity_log_partial << "on procID = " << crrProcID <<"\n";
	    }

	} 
    }

  posteriorDensity_partial_local = exp(posteriorDensity_log_partial);
  MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Reduce(&posteriorDensity_partial_local, &posteriorDensity_partial, 1, MPI_DOUBLE,  MPI_SUM, 0);
  MPI::COMM_WORLD.Barrier();
  posteriorDensity_log = log(posteriorDensity_partial);

  Eigen::Vector3d paraMax = im.get_paraMax();
  double priorPopTree = poptree->computeJointPrior(paraMax);
  //posteriorDensity *= priorPopTree/nSample;
  posteriorDensity_log += log(priorPopTree/nSample);

  posteriorDensity = exp(posteriorDensity_log);
  
  return posteriorDensity;	
}


// YC 4/15/2015
// Computing "summaryOfCoalTree" which contains
// the number of lineages in each population
// before each coalescent event without migration
// If lineages sampled from different population
// should be in the same population, the population 
// should be the ancestral population. Therefore, in
// this case, the number lineages in the sampling populations
// should be zero.
void Chain::compute_summaryOfCoalTrees(popTree* poptree)
{
  unsigned int nCoalEvents = coalTimes.at(0).at(0).size();
  unsigned int nUniqTopo = list_trees.size();
  unsigned int isThereAncPop = 0;
  unsigned int nSampledPops = poptree->size();
  if(poptree->get_age()!=0)
     isThereAncPop = 1;
  unsigned int nPops = nSampledPops+isThereAncPop;
  nodeSimple *topo = new nodeSimple;
  summaryOfCoalTree.resize(nUniqTopo);
  for(unsigned int id=0; id <nUniqTopo; id++)
    {
      summaryOfCoalTree.at(id).resize(nCoalEvents);	 
      topo = list_trees.at(id);
      topo->compute_areTipsFromSamePop();
      std::vector<unsigned int> nodePopIDs;
      unsigned int foundAll = 0;
      unsigned int Rank = 2*nCoalEvents+1;
      while(foundAll==0)
	{
	  if(Rank <= nCoalEvents+1)
	    foundAll = 1;
	  else
	    {
	      std::vector<unsigned int> ID= topo->get_nodePopID_withRank(Rank);
	      if(ID.at(1)==1)
		nodePopIDs.push_back(ID.at(0));
	      else
		{
		  std::cout << "\nError in Chain::compute_summaryOfCoalTrees()\n";
		  std::cout << "There is no node whose rank is same as Rank=" << Rank <<"\n";
		}
	    }
	  Rank--;
	}

      // computing the number of lineages in each population
      foundAll = 0;
      std::vector<unsigned int> nlin;
      nlin.resize(nSampledPops);
      std::vector<nodeSimple*> subtrees;
      subtrees.push_back(topo);
      unsigned int count_lin =0;
      while(foundAll==0)
	{
	  for(unsigned int i=0; i<subtrees.size(); i++)
	    {
	      if(subtrees.at(i)->get_isTip()==1)
		{
		  unsigned int sampledPopID = subtrees.at(i)->getPopID();
		  nlin.at(sampledPopID-1)++;
		  subtrees.erase(subtrees.begin()+i);
		  i--;
		  count_lin++;
		}
	      else
		{
		  subtrees.push_back(subtrees.at(i)->getFirstChild());
		  subtrees.push_back(subtrees.at(i)->getSecondChild());
		  subtrees.erase(subtrees.begin()+i);
		  i--;
		}
	      if(count_lin==nCoalEvents+1)
		foundAll=1;	      
	    }
	}
      
      summaryOfCoalTree.at(id).at(0).resize(nPops);
      for(unsigned int i=0; i<nSampledPops; i++)
	{
	  summaryOfCoalTree.at(id).at(0).at(i) = nlin.at(i);
	}
      if(poptree->get_age()!=0)
	summaryOfCoalTree.at(id).at(0).at(nPops-1) = nCoalEvents+1;

      for(unsigned int nc=1; nc<=nCoalEvents; nc++)
	{
	  summaryOfCoalTree.at(id).at(nc).resize(nPops);
	  summaryOfCoalTree.at(id).at(nc).at(nPops-1) =nCoalEvents+1-nc;
	  unsigned int nodeID = nodePopIDs.at(nc);
	  if(nodeID==0) // coalescent event between samples from different population
	    {
	      for(unsigned int i=0; i<nSampledPops; i++)
		summaryOfCoalTree.at(id).at(nc).at(i) = 0;
	    }
	  else
	    {
	      for(unsigned int i=0; i<nSampledPops; i++)
		{
		  if(nodeID == i+1)
		    summaryOfCoalTree.at(id).at(nc).at(i) = summaryOfCoalTree.at(id).at(nc).at(i-1)-1;
		  else		    
		    summaryOfCoalTree.at(id).at(nc).at(i) = summaryOfCoalTree.at(id).at(nc).at(i-1);
		}	      
	    }
	}
    } 
  
  delete topo;
  return;
}



// YC 4/30/2015
void Chain::compute_partialJointPosteriorDensity_overSubSample(popTree* poptree, IM im, unsigned int crrProcID, unsigned int nProcs)
{ 
  
  std::chrono::high_resolution_clock::time_point start_t, end_t;
  
  double migRateMax = im.get_migRateMax();

  // prepare_Lmode(poptree);
  if(migRateMax !=0) // not a single population
    {
      if(poptree->get_age()>0)
	compute_eigenValuesVectors_subMatOfEigenVectors(poptree);
      // compute_eigenValuesVectors_subMatOfEigenVectors_MPI(poptree);
    }
  else if(migRateMax==0)
    {
      compute_summaryOfCoalTrees(poptree);
    }
  

  std::chrono::microseconds initial(0);
  eachComputingTime_condiProb = initial;
  eachComputingTime_condiProb_case1= initial;
  eachComputingTime_condiProb_case2_1 = initial;
  eachComputingTime_condiProb_case2_2 = initial;
  eachComputingTime_condiProb_case3_1 = initial;
  eachComputingTime_condiProb_case3_2 = initial;

  expectationOfCoalProb.resize(n_loci);

  for(unsigned int j=0; j<n_loci; j++)
    {
      expectationOfCoalProb.at(j) = 0.0;
      
      for(unsigned int i=0; i< nSubSample; i++)
	{
	  
	  double eachterm_log = 0.0;
	  double eachLogProb = 0;
	  
	  eachLogProb = compute_logConditionalProb(i,j,poptree,crrProcID);

	  // the number of coalescent trees that have the give tree topology with the special labels
	  double ntrees = numTrees.at(treeIDs.at(i).at(j)); 
	  //std::cout << "ntrees = " << ntrees <<"\n";
	  
	  eachterm_log = eachLogProb-log(ntrees);
	  expectationOfCoalProb.at(j) += exp(eachterm_log);
	} 
      expectationOfCoalProb.at(j) /= n_MCMCgen;
    } 
  
  
 
  return;	
}



// YC 3/16/2016
void Chain::compute_partialJointPosteriorDensity_overSubLoci_ESS(popTree* poptree, IM im, unsigned int crrProcID, unsigned int nProcs)
{ 
  
  std::chrono::high_resolution_clock::time_point start_t, end_t;
  
  double migRateMax = im.get_migRateMax();

  // prepare_Lmode(poptree);
  if(migRateMax !=0) // not a single population
    {
      if(poptree->get_age()>0)
	compute_eigenValuesVectors_subMatOfEigenVectors_MPI(poptree);
    }
  else if(migRateMax==0)
    {
      compute_summaryOfCoalTrees(poptree);
    }
  

  std::chrono::microseconds initial(0);
  eachComputingTime_condiProb = initial;
  eachComputingTime_condiProb_case1= initial;
  eachComputingTime_condiProb_case2_1 = initial;
  eachComputingTime_condiProb_case2_2 = initial;
  eachComputingTime_condiProb_case3_1 = initial;
  eachComputingTime_condiProb_case3_2 = initial;

  logExpectationOfCoalProb.resize(numSubLoci);
  logExpectationOfCoalProbSquared.resize(numSubLoci);

  for(unsigned int j=0; j<numSubLoci; j++)
    {
      logExpectationOfCoalProb.at(j) = 0.0;
      logExpectationOfCoalProbSquared.at(j) = 0.0;
      
      if(TreesWithMaxP==1)
	{
	  unsigned int sampleID = MaxSampleID.at(j);
	  long double eachLogProb = (long double) compute_logConditionalProb(sampleID,j,poptree,crrProcID);
	  long double ntrees = (long double) numTrees.at(treeIDs.at(sampleID).at(j)); 
	  logExpectationOfCoalProb.at(j) = eachLogProb-log(ntrees);
	}
      else if(TreesWithMaxP==0)
	{
	  long double maxLogCondPr = -1*numeric_limits<long double>::infinity();
	  long double maxLogCondPrSquared = -1*numeric_limits<long double>::infinity();
	  std::vector<long double> logCondPr;
	  std::vector<long double> logCondPrSquared;
	  logCondPr.resize(nSubSample);
	  logCondPrSquared.resize(nSubSample);
	  for(unsigned int i=0; i< nSubSample; i++)
	    {
	      long double logPriorTree = 0.0;
	      if(im.get_priorType()==1)
		logPriorTree = (long double) logPrior_each.at(i).at(j);
	      long double eachLogProb = 0;
	      eachLogProb = (long double) compute_logConditionalProb(i,j,poptree,crrProcID);	      
	      // the number of coalescent trees that have the give tree topology with the special labels
	      long double ntrees = numTrees.at(treeIDs.at(i).at(j)); 
	      // eachterm_log = eachLogProb-log(ntrees);

	      logCondPr.at(i) = eachLogProb-log(ntrees)-logPriorTree;
	      logCondPrSquared.at(i) = 2*logCondPr.at(i);
	      // std::cout <<"crrProcID="<<crrProcID<<" sampleID i="<<i<<" logCondPr.at(i)="<<logCondPr.at(i)
	      //		<<" logCondPrSquared.at(i)="<<logCondPrSquared.at(i)<<"\n";

	      if(maxLogCondPr < logCondPr.at(i))
		maxLogCondPr = logCondPr.at(i);
	      if(maxLogCondPrSquared < logCondPrSquared.at(i))
		maxLogCondPrSquared = logCondPrSquared.at(i);
	    }
	  long double ratios =0.0;
	  long double ratiosSquared =0.0;
	  for(unsigned int i=0; i< nSubSample; i++)
	    {
	      ratios += exp(logCondPr.at(i)-maxLogCondPr);
	      ratiosSquared += exp(logCondPrSquared.at(i)-maxLogCondPrSquared);
	    } 
	  logExpectationOfCoalProb.at(j) = maxLogCondPr + log(ratios)-log(n_MCMCgen);
	  logExpectationOfCoalProbSquared.at(j) = maxLogCondPrSquared + log(ratiosSquared)-log(n_MCMCgen);
	}
      else
	{
	  std::cout << "\n*** Error in Chain::compute_partialJointPosteriorDensity_overSubLoci() *** \n";
	  std::cout << "TreesWithMaxP should be either 1 or 0, but TreesWithMaxP = " 
		    << TreesWithMaxP <<"\n\n";
	}
    } 
  
  
 
  return;	
}

// YC 1/4/2016
void Chain::compute_partialJointPosteriorDensity_overSubLoci(popTree* poptree, IM im, unsigned int crrProcID, unsigned int nProcs)
{ 
  
  std::chrono::high_resolution_clock::time_point start_t, end_t;
  
  double migRateMax = im.get_migRateMax();

  // prepare_Lmode(poptree);
  if(migRateMax !=0) // not a single population
    {
      if(poptree->get_age()>0)
	compute_eigenValuesVectors_subMatOfEigenVectors_MPI(poptree);
    }
  else if(migRateMax==0)
    {
      compute_summaryOfCoalTrees(poptree);
    }

  // std::cout << "Done with compute_eigenValuesVectors_subMatOfEigenVectors_MPI\n";
  

  std::chrono::microseconds initial(0);
  eachComputingTime_condiProb = initial;
  eachComputingTime_condiProb_case1= initial;
  eachComputingTime_condiProb_case2_1 = initial;
  eachComputingTime_condiProb_case2_2 = initial;
  eachComputingTime_condiProb_case3_1 = initial;
  eachComputingTime_condiProb_case3_2 = initial;

  logExpectationOfCoalProb.resize(numSubLoci);

  // std::cout << "numSubLoci = " << numSubLoci <<"\n";
  // std::cout << "nSubSample = " << nSubSample <<"\n";

  for(unsigned int j=0; j<numSubLoci; j++)
    {
      logExpectationOfCoalProb.at(j) = 0.0;
      
      if(TreesWithMaxP==1)
	{
	  unsigned int sampleID = MaxSampleID.at(j);
	  long double eachLogProb = (long double) compute_logConditionalProb(sampleID,j,poptree,crrProcID);
	  long double ntrees = (long double) numTrees.at(treeIDs.at(sampleID).at(j)); 
	  logExpectationOfCoalProb.at(j) = eachLogProb-log(ntrees);
	}
      else if(TreesWithMaxP==0)
	{
	  long double maxLogCondPr = -1*numeric_limits<long double>::infinity();
	  std::vector<long double> logCondPr;
	  logCondPr.resize(nSubSample);
	  for(unsigned int i=0; i< nSubSample; i++)
	    {
	      long double logPriorTree = 0.0;
	      if(im.get_priorType()==1)
		logPriorTree = (long double) logPrior_each.at(i).at(j);
	      
	      //std::cout << "crrProcID = " << crrProcID 
	      //    << "logPriorTree = " << logPriorTree <<"\n";
	      
	      long double eachLogProb = 0;
	      eachLogProb = (long double) compute_logConditionalProb(i,j,poptree,crrProcID);	      
	      // std::cout <<"eachLogProb = "<<eachLogProb <<"\n";
	      // the number of coalescent trees that have the give tree topology with the special labels
	      long double ntrees = numTrees.at(treeIDs.at(i).at(j)); 
	      // eachterm_log = eachLogProb-log(ntrees);

	      logCondPr.at(i) = eachLogProb-log(ntrees)-logPriorTree;
	      if(maxLogCondPr < logCondPr.at(i))
		maxLogCondPr = logCondPr.at(i);

	      // REMOVE
	      // std::cout << "ntrees = " << ntrees <<"\n";
	    }
	  long double ratios =0.0;
	  for(unsigned int i=0; i< nSubSample; i++)
	    {
	      // expectationOfCoalProb.at(j) += exp(eachterm_log);
	      ratios += exp(logCondPr.at(i)-maxLogCondPr);
	    } 
	  logExpectationOfCoalProb.at(j) = maxLogCondPr + log(ratios)-log(n_MCMCgen);
	  // expectationOfCoalProb.at(j) /= n_MCMCgen;
	}
      else
	{
	  std::cout << "\n*** Error in Chain::compute_partialJointPosteriorDensity_overSubLoci() *** \n";
	  std::cout << "TreesWithMaxP should be either 1 or 0, but TreesWithMaxP = " 
		    << TreesWithMaxP <<"\n\n";
	}
    } 
  
  
 
  return;	
}


// YC 1/6/2016
void Chain::compute_partialJointPosteriorDensity_mutationScalars_overSubLoci(popTree* poptree, IM im, unsigned int crrProcID, unsigned int nProcs)
{ 
  
  std::chrono::high_resolution_clock::time_point start_t, end_t;
  
  double migRateMax = im.get_migRateMax();

  if(migRateMax !=0) // not a single population
    {
      if(poptree->get_age()>0)
	compute_eigenValuesVectors_subMatOfEigenVectors_MPI(poptree);
    }
  else if(migRateMax==0)
    {
      compute_summaryOfCoalTrees(poptree);
    }
  

  std::chrono::microseconds initial(0);
  eachComputingTime_condiProb = initial;
  eachComputingTime_condiProb_case1= initial;
  eachComputingTime_condiProb_case2_1 = initial;
  eachComputingTime_condiProb_case2_2 = initial;
  eachComputingTime_condiProb_case3_1 = initial;
  eachComputingTime_condiProb_case3_2 = initial;

  // expectationOfCoalProb.resize(numSubLoci);
  
  if(TreesWithMaxP==1)
    partialLogJointCoalProb.resize(1);
  else if(TreesWithMaxP==0)
    partialLogJointCoalProb.resize(nSubSample);
  else
    {
      cout << "\n\n *** Errors in Chain::compute_partialJointPosteriorDensity_mutationScalars_overSubLoci()\n";
      cout << "TreesWithMaxP should be either 0 or 1, but TreesWithMaxP = " << TreesWithMaxP <<".\n\n";
    }
  
  for(unsigned int i=0; i< nSubSample; i++)
    {
      partialLogJointCoalProb.at(i) = 0.0;
      
      for(unsigned int j=0; j<numSubLoci; j++)
	{
	  double logPriorTree = 0.0;
	  if(im.get_priorType()==1)
	    logPriorTree = logPrior_each.at(i).at(j);
	  if(TreesWithMaxP==1 && i==0)
	    {
	      unsigned int sampleID = MaxSampleID.at(j);
	      long double eachLogProb = compute_logConditionalProb(sampleID,j,poptree,crrProcID);
	      double ntrees = numTrees.at(treeIDs.at(sampleID).at(j)); 
	      partialLogJointCoalProb.at(i) += eachLogProb-log(ntrees)-logPriorTree;
	    }
	  else if(TreesWithMaxP==0)
	    {
	      long double eachLogProb = compute_logConditionalProb(i,j,poptree,crrProcID);
	      double ntrees = numTrees.at(treeIDs.at(i).at(j)); 
	      partialLogJointCoalProb.at(i) += eachLogProb-log(ntrees)-logPriorTree;
	    }
	}
    }

 
  return;	
}



// YC 3/25/2015
double Chain::compute_partialJointPosteriorDensity_overSubSample_expectationOfJoint(popTree* poptree, IM im, unsigned int crrProcID, unsigned int nProcs)
{
  // std::cout << "In Chain::compute_partialJointPosteriorDensity_overSubSample\n";

  std::chrono::high_resolution_clock::time_point start_t, end_t;
  
  count_condiProbFunctionCalls =0;
  count_case1 =0;
  count_case2_1=0;
  count_case2_2=0;
  count_case3_1=0;
  count_case3_2=0;

  double migRateMax = im.get_migRateMax();

  prepare_Lmode(poptree);
  if(migRateMax !=0) // not a single population
    {
      // prepare_Lmode(poptree);
      if(poptree->get_age()>0)
	compute_eigenValuesVectors_subMatOfEigenVectors_MPI(poptree);
    }
  else if(migRateMax==0)
    {
      compute_summaryOfCoalTrees(poptree);
    }
  
  double posteriorDensity_log = 0.0;
  double posteriorDensity_partial = 1;
  double posteriorDensity_partial_local = 1;
  double posteriorDensity_log_partial = 0.0;
  double posteriorDensity = 0.0;


  std::chrono::microseconds initial(0);
  eachComputingTime_condiProb = initial;
  eachComputingTime_condiProb_case1= initial;
  eachComputingTime_condiProb_case2_1 = initial;
  eachComputingTime_condiProb_case2_2 = initial;
  eachComputingTime_condiProb_case3_1 = initial;
  eachComputingTime_condiProb_case3_2 = initial;
  for(unsigned int i=0; i< nSubSample; i++)
    {
      double eachterm_log = 0.0;
      
      for(unsigned int j=0; j<n_loci; j++)
	{
	  // std::cout << "i= "<< i <<"j="<<j <<"\n";
	  double eachLogProb = 0;
	  
	  eachLogProb = compute_logConditionalProb(i,j,poptree,crrProcID);
	  
	  double ntrees = numTrees.at(treeIDs.at(i).at(j)); // the number of coalescent trees that have the give tree topology with the special labels
	  //std::cout << "ntrees = " << ntrees <<"\n";
	  
	  eachterm_log += eachLogProb-log(ntrees);
	      
	  if(eachterm_log == numeric_limits<double>::infinity())// || eachterm_log == -1*numeric_limits<double>::infinity())
	    {
	      std::cout << "\n*** Error *** in Chain::compute_partialJointPosteriorDensity_overSubSample()\n";
	      std::cout << "\t on " << i << "th sample and " << j<<"th loci: eachterm_log = " << eachterm_log <<"\n";
	      std::cout << "\t eachLogProb =  " << eachLogProb <<" and log(ntrees) = " << log(ntrees) <<"\n";
	      std::cout << "\t population tree is \n";
	      poptree->print_allPopTree();
	      std::cout <<"\n";
	    }
	  // YC 2/10/2015
	  // Here I need to add the following: if eachLogProb = -inf, then the probability 
	  // of the given coal. tree is (very close to) zero and then the joint prob 
	  // of coal. trees from all loci, which is the product of individual prob, is 
	  // very close to zero anyway. 	      
	  if(eachLogProb == -1*numeric_limits<double>::infinity())
	    {
	      j = n_loci;
	      eachterm_log =-1*numeric_limits<double>::infinity();
	    }
	  
	} 
      if(posteriorDensity_log_partial == 0.0)
	{
	  if(eachterm_log == -1*numeric_limits<double>::infinity())
	    posteriorDensity_log_partial = eachterm_log;
	  else
	    posteriorDensity_log_partial = eachterm_log-logPrior_trees.at(i);
	}
      else if(exp(eachterm_log) == 0 || exp(posteriorDensity_log_partial) == 0)
	{
	  double logLarge, logSmall;
	  if(posteriorDensity_log_partial < eachterm_log)
	    {
	      logLarge = eachterm_log-logPrior_trees.at(i);
	      logSmall = posteriorDensity_log_partial;
	    }
	  else
	    {
	      logLarge = posteriorDensity_log_partial;
	      logSmall = eachterm_log-logPrior_trees.at(i);
	    }
	  posteriorDensity_log_partial = logLarge+log(1+exp(logSmall-logLarge));
	  if(posteriorDensity_log_partial == numeric_limits<double>::infinity())
	    std::cout << "case1: posteriorDensity_log_partial = " << posteriorDensity_log_partial << " on procID = " << crrProcID <<"\n";
	}
      else
	{
	  double factor2sum = exp(eachterm_log -logPrior_trees.at(i));
	  posteriorDensity_log_partial = log(exp(posteriorDensity_log_partial)+factor2sum);
	  if(posteriorDensity_log_partial == numeric_limits<double>::infinity())
	    std::cout << "case2: posteriorDensity_log_partial = " << posteriorDensity_log_partial << "on procID = " << crrProcID <<"\n";
	}         
    } 
  
  posteriorDensity_log = posteriorDensity_log_partial;
  
  Eigen::Vector3d paraMax = im.get_paraMax();
  double priorPopTree = poptree->computeJointPrior(paraMax);
  if(posteriorDensity_log== -1*numeric_limits<double>::infinity())
    {
      posteriorDensity = 0.0;
    }
  else
    {
      posteriorDensity_log += log(priorPopTree/n_MCMCgen);
      
      posteriorDensity = exp(posteriorDensity_log);
    }
 
  return posteriorDensity;	
}



















double Chain::compute_jointPosteriorDensity_MPI(popTree* poptree, IM im, unsigned int crrProcID)
{
  //compute_eigenValuesVectors_rateMat_MPI(poptree,crrProcID);
  compute_eigenValuesVectors(poptree);

  double posteriorDensity_log = 0.0;
  double posteriorDensity_log_partial = 0.0;
  double posteriorDensity = 0.0;
  
  unsigned int nSample = im.get_nSampledTrees();
  unsigned int nloci = im.get_nLoci();
  
  for(unsigned int i=0; i< nSample ; i++)
    {
      double eachterm_log = 0.0;
		// double eachterm = 1.0;
      for(unsigned int j=0; j<nloci; j++)
	{
	  //std::cout <<"Compute the conditional prob for locus "<< j <<" from " << i<<"th generation on process "<< crrProcID <<"..\n";
	  double eachProb = compute_conditionalProb_MPI(i,j,poptree,crrProcID);
	  double ntrees = numTrees.at(treeIDs.at(i).at(j));
	  
	  //std::cout << "Done with computing the conditional prob for locus "<< j <<" from " << i<<"th generation on process "<< crrProcID <<".\n";
	  //eachterm *= eachProb/ntrees;
	  eachterm_log += log(eachProb)-log(ntrees);
	  
	  if(eachterm_log == numeric_limits<double>::infinity())
	    {
	      std::cout << "*** Error *** in Chain::compute_jointPosteriorDensity_MPI()\n";
	      std::cout << "\t on " << i << "th sample and " << j<<"th loci: eachterm_log = " << eachterm_log <<"\n";
		    std::cout << "\t log(eachProb) =  " << log(eachProb) <<" and log(ntrees) = " << log(ntrees) <<"\n";
		  }
		}
	
		posteriorDensity += exp(eachterm_log -logPrior_trees.at(i));
	
	}


	Eigen::Vector3d paraMax = im.get_paraMax();
	double priorPopTree = poptree->computeJointPrior(paraMax);
	posteriorDensity *= priorPopTree/nSample;
	
	return posteriorDensity;

}

