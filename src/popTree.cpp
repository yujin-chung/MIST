/* MIST copyright 2016 by Yujin Chung and Jody Hey */

/*
 * popTree.cpp
 */

#include <iostream>
#include "popTree.hpp"
#include "misc.hpp"

popTree::popTree()
{
	isRoot = 2;
	isTip = 2;
	popID = 0;
	par = 0;
	desc[0] = 0; desc[1] = 0;
	age = 0.0;
	populationSize = 0.0;
	no_assignedChildren = 3;
	// mig;
}



void popTree::initialization(IM im)
{
  isRoot = 2;
  isTip = 2;
  popID = 0;
  par = 0;
  desc[0] = 0; desc[1] = 0;
  age = 0.0;
  populationSize = 0.0;
  no_assignedChildren = 3;
  
  // Added by YC 5/9/2014
  pop2mig.resize(0);
  migRate.resize(0);
  
  ancPop = im.get_ancPop();
  samePopulationSizes = im.get_samePopulationSizes();
  sameMigrationRates = im.get_sameMigrationRates();
}


void popTree::deletePopTree()
{
	if(isTip == 0 )
	{
		desc[0]->deletePopTree();
		desc[1]->deletePopTree();
	}
	delete desc[0];
	delete desc[1];
	//delete par;
}




// void popTree::replacePara(unsigned int ancPop, Eigen::MatrixXd listPara)
void popTree::replacePara(Eigen::MatrixXd listPara)
{
  populationSize = listPara(0,2);
  age = listPara(0,5);
  
  if(age>0) // not a single population
    {      
      // 2018/07/17 YC
      durationOfSplitting = listPara(0,6);
      
      double popSize = listPara(0,0);
      desc[0]->assign_populationSize(popSize);
      popSize = listPara(0,1);
      desc[1]->assign_populationSize(popSize);
      std::vector<popTree*> pops;
      std::vector<double> mig;
      double mRate =listPara(0,3);
      mig.resize(0); pops.resize(0);
      mig.push_back(mRate);
      pops.push_back(desc[1]);
      desc[0]->replace_migRate(mig);
      desc[0]->replace_pop2mig(pops);
      mig.at(0) = listPara(0,4);
      pops.at(0) = desc[0];
      desc[1]->replace_migRate(mig);
      desc[1]->replace_pop2mig(pops);
    }
  /*
  if(ancPop == 0) // island model (no ancestral population)
    {
      if(listPara.cols() == 2) // same population sizes and same migration sizes
	{
	  double popSize = listPara(0,0);
	  desc[0]->assign_populationSize(popSize);
	  desc[1]->assign_populationSize(popSize);
	  std::vector<popTree*> pops;
	  std::vector<double> mig;
	  double mRate =listPara(0,1);
	  mig.push_back(mRate);
	  pops.push_back(desc[1]);
	  desc[0]->replace_migRate(mig);
	  desc[0]->replace_pop2mig(pops);
	  pops.at(0) = desc[0];
	  desc[1]->replace_migRate(mig);
	  desc[1]->replace_pop2mig(pops);	  
	}
    }
  else if(ancPop == 1) // isolation model
    {
      if(listPara.cols() == 6)
	{
	  populationSize = listPara(0,2);
	  age = listPara(0,5);
	  double popSize = listPara(0,0);
	  desc[0]->assign_populationSize(popSize);
	  popSize = listPara(0,1);
	  desc[1]->assign_populationSize(popSize);
	  std::vector<popTree*> pops;
	  std::vector<double> mig;
	  double mRate =listPara(0,3);
	  mig.push_back(mRate);
	  pops.push_back(desc[1]);
	  desc[0]->replace_migRate(mig);
	  desc[0]->replace_pop2mig(pops);
	  mig.at(0) = listPara(0,4);
	  pops.at(0) = desc[0];
	  desc[1]->replace_migRate(mig);
	  desc[1]->replace_pop2mig(pops);
	}
    }

  */

  return;
}



void popTree::initialize_popTree_recursion(IM im, std::string newickTree, Eigen::Vector3d paraMax)
{
  double para = 0.0;
  
  // 2018/07/16
  durationOfSplitting = im.get_durationOfSplitting();
  migband= im.get_migband();
  
  if(newickTree.size()>0 && newickTree.compare(0,1,";")!=0)
    {
      if(newickTree.compare(0,1,"(")==0)
	{
	  desc[no_assignedChildren] = new popTree;
	  desc[no_assignedChildren]->assign_isRoot_isTip(0,0);

	  //-- splitting time --//
	  if(im.get_ancPop() == 1) // isolation model
	    {
	      // 2018/07/16 - YC
	      if(migband == 1) // estimated
		{
		  para = runiform()*paraMax(1);
		  desc[no_assignedChildren]->assign_durationOfSplitting(runiform()*para);
		}
	      else if(migband == 2) //fixed
		{
		  para = durationOfSplitting + runiform()*(paraMax(1)-durationOfSplitting);
		  if(para < 0)
		    {
		      std::cout <<"\n\n Error: The upperbound of population splitting time should be larger than the duration of splitting.\n\n";
		    }			    
		}
	    }
	  else // island model (no ancestral population)
	    {
	      para = paraMax(1);
	    }
	  desc[no_assignedChildren]->assign_age(para);
	 

	  //-- population size --//	  
	  if(im.get_ancPop() == 1) // isolation model	
	    {
	      if(im.get_samePopulationSizes() == 1) // same population sizes
		para = populationSize;
	      else // allow different population sizes
		para = runiform() * paraMax(0);
	    }
	  else // island model (no ancestral population)
	    para =0;
	  desc[no_assignedChildren]->assign_populationSize(para);

	  desc[no_assignedChildren]->assign_popID(get_popID()+1);
	  desc[no_assignedChildren]->no_assignedChildren = 0;
	  desc[no_assignedChildren]->assign_par(this);
	  desc[no_assignedChildren]->initialize_popTree_recursion(im,newickTree.erase(0,1),paraMax);
	  no_assignedChildren++;
	}
      else if(newickTree.compare(0,1,")")==0)
	{	  
	  par->initialize_popTree_recursion(im,newickTree.erase(0,1),paraMax);
	}
      else if(newickTree.compare(0,1,",")==0)
	{
	  initialize_popTree_recursion(im,newickTree.erase(0,1),paraMax);
	}
      else
	{
	  desc[no_assignedChildren] = new popTree;
	  desc[no_assignedChildren]->assign_isRoot_isTip(0,1);
	  desc[no_assignedChildren]->assign_age(0.0);

	  //-- population size --//
	  if(im.get_samePopulationSizes() == 0) // allow different population sizes
	    para = runiform() * paraMax(0);
	  else
	    {
	      if(no_assignedChildren == 0)
		{
		  if(im.get_ancPop() == 1) // isolation model
		    para = populationSize;
		  else // island model (no ancestral population)
		    para = runiform() * paraMax(0);		    
		}
	      else
		para = desc[0]->get_popSize();		
	    }
	  desc[no_assignedChildren]->assign_populationSize(para);

	  char name = newickTree.at(0);
	  name = name - '0';
	  int n = name;
	  desc[no_assignedChildren]->assign_popID((unsigned)n);
	  desc[no_assignedChildren]->assign_par(this);
	  no_assignedChildren++;
	  initialize_popTree_recursion(im, newickTree.erase(0,1),paraMax);
	}
    }
  return;
}

/**
 * @param paraMax contains (popSizeMax,splittingTimeMax,migRateMax)
 */
void popTree::initialize_popTree(IM im, unsigned int processID)
{
  std::string newickTree = im.get_poptree_string();
  unsigned int ancPop = im.get_ancPop();
  Eigen::Vector3d paraMax = im.get_paraMax();

  // 2018/07/16
  durationOfSplitting = im.get_durationOfSplitting();
  migband= im.get_migband();
  
  if(processID == 0)
    {
	std::cout << "Initialization of population tree....\n";
    }

  double para = 0.0;
  if(newickTree.compare(0,1,"(")==0)
    {
      
      newickTree.erase(0,1);
      newickTree.erase(newickTree.size()-1,1);

      assign_isRoot_isTip(1,0);
      char name = newickTree.at(newickTree.size()-1);
      if(newickTree.compare(newickTree.size()-1,1,")")==0)
	{
	  name = newickTree.at(newickTree.size()-2);
	}
      name = name - '0';
      int n = name;
     
      popID = (unsigned) n;
      //-- splitting time --//
      if(im.get_ancPop() == 1) // isolation model
	{
	  if(name <= 2) // if ancestral popID =2, it is the case of a sigle population
	    {
	      age = 0;
	      isTip =1;
	      durationOfSplitting = 0;
	    }
	  else // single population
	    {
	      age = runiform()*paraMax(1);
	    }
	}
      else // island model (no ancestral population)
	{
	  age = paraMax(1);
	}

      //-- ancestral population size --//      
      if(im.get_ancPop() == 1) // isolation model
	para = runiform() * paraMax(0);
      else // island model (no ancestral population)
	para = 0;
      populationSize=para;
      no_assignedChildren = 0;
      
      
      newickTree.erase(newickTree.size()-2,2);
      
      if(age > 0)
	initialize_popTree_recursion(im, newickTree, paraMax);
    }
  else
    std::cout << "The input population tree is not in newick format.\n";
  
  // Initialize migration rates
  if(age>0)
    {
      para = paraMax(2);
      initialize_migrations(im,para, processID);
    }
  if(processID==0)
    {
      if(age>0)
	{
	  std::cout << "From\tTo\tMigrationRate\n";
	  std::cout << desc[0]->get_popID() << "\t" << desc[0]->get_pop2mig(0)->get_popID()
		    <<"\t" << desc[0]->get_migRate(0) <<"\n";
	  std::cout << desc[1]->get_popID() << "\t" << desc[1]->get_pop2mig(0)->get_popID()
		    <<"\t" << desc[1]->get_migRate(0) <<"\n";
	}
      print_poptree();
      
      std::cout << "End of population tree initialization.\n\n";
    }
  
  return;
}



void popTree::initialize_migrations_recursion(IM im, std::vector<popTree*> pops, double rateMax)
{
  // FIXME
  // now it works for 2 population IM model only.

  if(pops.size()==2)
    {
      double rate = runiform()*rateMax;
      pops.at(0)->add_pop2mig(pops.at(1));
      pops.at(0)->add_migRate(rate);
      if(im.get_sameMigrationRates() == 0) // allow different migration rates
	rate = runiform()*rateMax;
      pops.at(1)->add_pop2mig(pops.at(0));
      pops.at(1)->add_migRate(rate);
    }
}

void popTree::initialize_migrations(IM im, double rateMax, unsigned int processID)
{

  std::vector<popTree*> pops;
  pops.push_back(desc[0]);
  pops.push_back(desc[1]);
  initialize_migrations_recursion(im,pops,rateMax);

}





popTree* popTree::deepCopy_root()
{
	popTree* copyTr = new popTree;

	copyTr->isTip = isTip;
	copyTr->isRoot = isRoot;
	copyTr->popID = popID;
	copyTr->age = age;
	copyTr->populationSize = populationSize;
	copyTr->no_assignedChildren = no_assignedChildren;
	copyTr->mig = mig;
	if(isTip == 0)
	{
		copyTr->desc[0] = desc[0]->deepCopy_root();
		copyTr->desc[1] = desc[1]->deepCopy_root();
		copyTr->desc[0]->par = copyTr;
		copyTr->desc[1]->par = copyTr;
	}
	return copyTr;
}



void Migration::initialization(unsigned int ID, double rate)
{
	popID = ID;
	migRate = rate;
}

void Migration::print()
{
	std::cout << "The receiver popID: " << popID
			<< ", migration rate: " << migRate <<"\n";
}




// old version - YC 7/2/2014
/*
void popTree::print_migrationInfo()
{
	std::cout << "The current population ID is " << popID <<"\n";
	for(unsigned int i=0; i<mig.size(); i++)
		mig.at(i).print();
}
*/


unsigned int popTree::size()
{
	unsigned int s = 0;
	if(isTip)
		s++;
	else
	{
		s+= desc[0]->size();
		s+= desc[1]->size();
	}
	return s;
}

popTree* popTree::go2root()
{
	if(par->isRoot)
		return par;
	else
		par->go2root();
}

double popTree::get_minSplittingTime(double time)
{
	if(isTip == 0)
	{
		if(age > time)
			time = age;
		time = desc[0]->get_minSplittingTime(time);
		time = desc[1]->get_minSplittingTime(time);
	}

	return time;
}

int popTree::getID_beforeSplit(int id, double SplittingT, int &found)
{
	int newID = 0;
	if(age == SplittingT)
	{
		if(desc[0]->popID == id || desc[1]->popID == id)
		{
			newID = popID;
			found = 1;
		}
	}
	else if(isTip == 0)
	{
		int tmpID = 0;
		if(desc[0]->age < SplittingT)
		{
			newID = popID;
			found = 1;
		}
		else
		{
			tmpID = desc[0]->getID_beforeSplit(id, SplittingT, found);
			if(found ==1)
				newID = tmpID;
			else
			{
				tmpID = desc[1]->getID_beforeSplit(id, SplittingT, found);
				if(found == 1)
					newID = tmpID;
			}
		}
	}
	return newID;
}

void popTree::print_popSize()
{
	if(isRoot == 1)
		std::cout << "Population sizes:\n";

	std::cout << popID << ": " << populationSize <<"\n";

	if(isTip == 0)
	{
		desc[0]->print_popSize();
		desc[1]->print_popSize();
	}
}



void popTree::print_node()
{
	if(isTip == 1)
		std::cout << popID << ":" << par->age;
	else
	{
		std::cout << "(";
		desc[0]->print_node();
		std::cout << ",";
		desc[1]->print_node();
		std::cout << ")"<< popID <<":" << par->age - age ;
	}
}


void popTree::print_poptree()
{
	if(isTip == 1)
		std::cout << "(" << popID << ");\n";
	else
	{
		std::cout << "(";
		try{
			desc[0]->print_node();
		}catch (std::exception &e) {
			std::cout << "In popTree::print_poptree(). Can't access element desc[0] \n";
		}

		std::cout << ",";
		try{
			desc[1]->print_node();
		}catch (std::exception &e) {
			std::cout << "In popTree::print_poptree(). Can't access element desc[1] \n";
		}
		std::cout << ")"<< popID <<";\n";
	}

	//REMOVE
	//std::cout << "testing\n";

	print_popSize();
}

// print migration information
void popTree::print_migRate()
{
  if(isRoot == 1)
    {
      std::cout << "Migration rates\n";
    }
  if(migRate.size()>0)
    {
      for(unsigned int i=0; i<migRate.size(); i++)
	{
	  std::cout << "From pop " << popID << " to pop "<< pop2mig.at(i)->get_popID() << ": " << migRate.at(i) << "\n";
	}
    }
  if(isTip == 0)
    {
      desc[0]->print_migRate();
      desc[1]->print_migRate();
    }
}

void popTree::print_allPopTree()
{
  print_poptree();
  print_migRate();
}



/**
 * @param paraMax contains (popSizeMax,splittingTimeMax,migRateMax)
 This does not compute the exact prior of population tree, but this value is constant over population tree space. That is, the exact prior is proportional to this value.
 */
double popTree::computeJointPrior(Eigen::Vector3d paraMax)
{
  double prior = 1;
  //if(paraMax(0) != 0.0)
  //  prior = 1/paraMax(0);
  unsigned int nSamplingPops = size();
  unsigned int nTotalPops = (int) nSamplingPops*(nSamplingPops+1)/2;
  if(ancPop==1 && size()>1 && paraMax(1) != 0.0)
    prior *= 1/paraMax(1);

  if(paraMax(0) != 0.0)
    {
      if(samePopulationSizes ==0) // allow different population sizes
	{
	  if(ancPop == 0) // no ancestral populations
	    prior *= 1/pow(paraMax(0),nSamplingPops);
	  else		
	    prior *= 1/pow(paraMax(0),nTotalPops);		
	}
      else // same population sizes
	prior *= 1/paraMax(0); 
    }

  if(paraMax(2) !=0.0)
    {
      if(sameMigrationRates == 0) // allow different migration rates
	{
	  prior *= 1/pow(paraMax(2), nSamplingPops*(nSamplingPops-1) );
	}
      else // same migration rates
	prior *= 1/paraMax(2);
    }
 
  return prior;
}

double popTree::find_migrationRate(unsigned int popID_from, unsigned int popID_to)
{
  double mig = -1;
  if(popID == popID_from)
    {
      unsigned int count = 0;
      unsigned int found = 0;
      while(found ==0 && count < pop2mig.size())
	{
	  if(popID_to == pop2mig.at(count)->get_popID())
	    {
	      mig = migRate.at(count);
	      found = 1;
	    }
	  count++;
	}
    }
  else if(isTip == 0)
    {
      mig = desc[0]->find_migrationRate(popID_from, popID_to);
      if(mig < 0)
	{
	  mig = desc[1]->find_migrationRate(popID_from, popID_to);
	}
    }
  
  return mig;
}

double popTree::find_popSize(unsigned int pop_ID)
{
	double popSize = -1;
	if(popID == pop_ID)
	{
		popSize = populationSize;
	}
	else if(isTip==0)
	{
		popSize = desc[0]->find_popSize(pop_ID);
		if(popSize < 0)
		{
			popSize =  desc[1]->find_popSize(pop_ID);
		}
	}

	return popSize;

}


