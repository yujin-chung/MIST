/* MIST copyright 2016 by Yujin Chung and Jody Hey */

/*
 * popTree.cpp
 */

#include <iostream>
#include <algorithm>    // std::max
#include "popTree.hpp"
#include "misc.hpp"

popTree::popTree()
{
	isRoot = 2;
	isTip = 2;
	popID = 0;
	par = 0;
	// desc[0] = 0; desc[1] = 0;
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
  // desc[0] = 0; desc[1] = 0;
  desc.resize(0);
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
	  desc.at(0)->deletePopTree();
	  desc.at(1)->deletePopTree();
	}
	delete desc.at(0);
	delete desc.at(1);
	//delete par;
}




// void popTree::replacePara(unsigned int ancPop, Eigen::MatrixXd listPara)
void popTree::replacePara(Eigen::MatrixXd listPara)
{
  populationSize = listPara(0,2);
  age = listPara(0,5);
  if(age>0) // not a single population
    {
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
  #ifdef DEBUG
  std::cout <<"\nIn popTree::initialize_popTree_recursion()\n";
  std::cout <<"newickTree = "<<newickTree <<"\n";
#endif// DEBUG

  
  isRoot = 0;
  std::size_t loc_openP = newickTree.find_first_of("(");
  std::size_t loc_closeP = newickTree.find_first_of(")");
  std::size_t loc_comma = newickTree.find_first_of(",");
  if(loc_comma == loc_openP)
    {
      isTip = 1;
      popID = (unsigned) atoi(newickTree.c_str());      
    }
  else if(loc_comma < loc_openP)
    {
      std::string child = newickTree.substr(0, loc_comma);
      std::string newSubtr = newickTree.substr(loc_comma+1, newickTree.size()-loc_comma);
      popTree *ch = new popTree;
      desc.push_back(ch);     
      desc.at(no_assignedChildren)->initialize_popTree_recursion(im,child,paraMax);
      no_assignedChildren++;
      // the remaining string may have one or more children.
      initialize_popTree_recursion(im,newSubtr,paraMax);
    }
  else
    {
    }
  double para = 0.0;
  if(newickTree.size()>0 && newickTree.compare(0,1,";")!=0)
    {
      if(newickTree.compare(0,1,"(")==0)
	{
	  popTree *poptr = new popTree;
	  desc.push_back(poptr);
	  // desc[no_assignedChildren] = new popTree;
	  
	  desc.at(no_assignedChildren)->assign_isRoot_isTip(0,0);

	  //-- splitting time --//
	  if(im.get_ancPop() == 1) // isolation model			
	    para = runiform()*paraMax(1);
	  else // island model (no ancestral population)
	    para = paraMax(1);
	  desc.at(no_assignedChildren)->assign_age(para);

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
	  desc.at(no_assignedChildren)->assign_populationSize(para);

	  desc.at(no_assignedChildren)->assign_popID(get_popID()+1);
	  desc.at(no_assignedChildren)->no_assignedChildren = 0;
	  desc.at(no_assignedChildren)->assign_par(this);
	  desc.at(no_assignedChildren)->initialize_popTree_recursion(im,newickTree.erase(0,1),paraMax);
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
	  popTree *poptr = new popTree;
	  desc.push_back(poptr);
	  // desc[no_assignedChildren] = new popTree;
	  desc.at(no_assignedChildren)->assign_isRoot_isTip(0,1);
	  desc.at(no_assignedChildren)->assign_age(0.0);

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
		para = desc.at(0)->get_popSize();		
	    }
	  desc.at(no_assignedChildren)->assign_populationSize(para);

	  char name = newickTree.at(0);
	  name = name - '0';
	  int n = name;
	  desc.at(no_assignedChildren)->assign_popID((unsigned)n);
	  desc.at(no_assignedChildren)->assign_par(this);
	  no_assignedChildren++;
	  initialize_popTree_recursion(im, newickTree.erase(0,1),paraMax);
	}
    }
  return;
}

/**
 * @param paraMax contains (popSizeMax,splittingTimeMax,migRateMax)
 */
/**
 * ((1,2):4,3):5; Isolation models with and without migrations
 * (1,2,3):4; Isolation models but the ancestor has three descendants.
 * (1,2,3); Island models (with and without migrations)
 * ((1,2):4,3); A combination of isolation and island models.
 */
void popTree::initialize_popTree(IM im, unsigned int processID)
{
  #ifdef DEBUG
  std::cout <<"In popTree::initialize_popTree(IM im, unsigned int processID)\n";
#endif //DEBUG
  std::string newickTree = im.get_poptree_string();
  unsigned int ancPop = im.get_ancPop();
  Eigen::Vector3d paraMax = im.get_paraMax();

  /*
  if(processID == 0)
    {
	std::cout << "Initialization of population tree....\n";
    }
  */

  #ifdef DEBUG
  std::cout <<"newickTree = "<<newickTree <<"\n";
#endif// DEBUG

  double para = 0.0;
  if(newickTree.compare(0,1,"(")==0 & newickTree.compare(newickTree.size()-1,1,";")==0)
    {
      
      newickTree.erase(0,1);
      newickTree.erase(newickTree.size()-1,1);
	  
      isRoot = 1;
  #ifdef DEBUG
      // std::cout <<"newickTree = "<<newickTree <<"\n";
#endif// DEBUG
      
      std::size_t loc_comma = newickTree.find(",");
      if(loc_comma == string::npos)
	{
	  // single population e.g, (1);
	  isTip = 1;
	  age = 0;
	  populationSize= runiform() * paraMax(0);
	  
	  newickTree.erase(newickTree.size()-1,1);
	  popID = (unsigned) atoi(newickTree.c_str());
	  no_assignedChildren = 0;
	  
	  #ifdef DEBUG
	  std::cout <<"Single population\n";
	  std::cout << "popID = "<<popID <<"\n";
#endif// DEBUG
	}
      else
	{
	  if(newickTree.compare(newickTree.size()-1,1,")") == 0)
	    {
	      /** island model or a mix
		  e.g. (1,2,3); ((1,2):4,3); **/
	      isTip = 0;
	      age = numeric_limits<long double>::infinity();
	      populationSize= 0;
	      popID = 0;
	    }
	  else
	    {
	      /** islation models e.g., ((1,2):4,3):5; (1,2,3):4; **/
	      isTip = 0;
	      age = runiform()*paraMax(1);
	      populationSize= runiform() * paraMax(0);
	      
	      std::size_t loc = newickTree.find_last_of(":");
	      std::string str = newickTree.substr(loc+1,newickTree.size()-loc);
	      popID = (unsigned) atoi(str.c_str());  
	      
	      newickTree.erase(loc-1,str.size()+2);
	    }
#ifdef DEBUG
	      std::cout << "popID = "<< popID <<"\n";
	      std::cout <<"newickTree = "<<newickTree <<"\n";
#endif //DEBUG
	      no_assignedChildren = 0;

	      
	      std::size_t loc_comma = newickTree.find(",");
	      std::size_t loc_openP = newickTree.find("(");
	      std::size_t loc_closeP = newickTree.find(")");	      
	      /*
	      int loc_comma = newickTree.find(",");
	      int loc_openP = newickTree.find("(");
	      int loc_closeP = newickTree.find(")");
	      */
	      int count_open_minus_close = 0;	      
	      int loc_max =0;
	      while(newickTree.size()!=0)
		{
		  #ifdef DEBUG
		  /*
		  std::cout <<"loc_comma = "<< loc_comma <<" loc_openP =" << loc_openP <<" loc_closeP = "<< loc_closeP
			    <<" count_open_minus_close = "<< count_open_minus_close
			    << " newickTree = "<< newickTree
			    <<"\n";
		  */
#endif //DEBUG
		  // if(loc_comma < loc_openP & loc_openP < loc_closeP)
		  if(loc_openP < loc_closeP)
		    {
		      count_open_minus_close++;	
		      std::string stmp = newickTree.substr(loc_openP+1,newickTree.size()-loc_openP);
		      if(stmp.find("(")==std::string::npos)
			loc_openP = std::string::npos;
		      else
			loc_openP = 1+loc_openP + stmp.find("(");
		      #ifdef DEBUG
		      /*
		      std::cout << stmp << " loc of ( is "<< stmp.find("(")
			<<"\n";
		      */
#endif //DEBUG
		    }
		  else if(loc_openP > loc_closeP)
		    {
#ifdef DEBUG
		      // std::cout << newickTree.substr(loc_closeP+1,newickTree.size()-loc_closeP)<<"\n";
#endif //DEBUG
		      count_open_minus_close--;
		      std::string stmp =newickTree.substr(loc_closeP+1,newickTree.size()-loc_closeP);   
		      if(count_open_minus_close!=0)
			{
			  if(stmp.find(")")==std::string::npos)
			    loc_closeP = std::string::npos;
			  else
			    loc_closeP = 1+loc_closeP + stmp.find(")");
			}
		    }
		  if(count_open_minus_close ==0)
		    {
		      
		      std::string stmp =newickTree.substr(loc_closeP+1,newickTree.size()-loc_closeP);   
		      if(stmp.find(",")==std::string::npos)
			loc_comma = newickTree.size();
		      else
			loc_comma = 1+loc_closeP + stmp.find(","); 
		      std::string child = newickTree.substr(0, loc_comma);
		      #ifdef DEBUG
		      std::cout <<"child = "<< child <<"\n";
#endif// DEBUG
		      popTree *ch = new popTree;
		      desc.push_back(ch);     
		      //  desc.at(no_assignedChildren)->initialize_popTree_recursion(im,child,paraMax);
		      no_assignedChildren++;
		      newickTree.erase(0,loc_comma+1); // loc_comma+1, newickTree.size()-loc_comma);
		      loc_comma = newickTree.find_first_of(",");
		      loc_openP = newickTree.find_first_of("(");
		      loc_closeP = newickTree.find_first_of(")");
		      // the remaining string may have one or more children.
		      //initialize_popTree_recursion(im,newSubtr,paraMax);
		    }
		  /*
		  else 
		    {
		      int loc_max = std::max(std::max(loc_comma, loc_closeP), loc_openP);
		      #ifdef DEBUG
		      std::cout <<"loc_max = "<< loc_max;
#endif // DEBUG			
		      loc_comma = loc_max+newickTree.substr(loc_max+1,newickTree.size()-loc_max).find(",");
		    }
		  */
		}
	}      
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
	  std::cout << desc.at(0)->get_popID() << "\t" << desc.at(0)->get_pop2mig(0)->get_popID()
		    <<"\t" << desc.at(0)->get_migRate(0) <<"\n";
	  std::cout << desc.at(1)->get_popID() << "\t" << desc.at(1)->get_pop2mig(0)->get_popID()
		    <<"\t" << desc.at(1)->get_migRate(0) <<"\n";
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


