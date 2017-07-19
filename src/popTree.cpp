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



void popTree::initialize_popTree_recursion(IM im, std::string newickTree, Eigen::Vector3d paraMax, popTree* root)
{
  #ifdef DEBUG
  /*
  std::cout <<"\nIn popTree::initialize_popTree_recursion()\n";
  std::cout <<"newickTree = "<<newickTree <<"\n";
  */
#endif// DEBUG

  
  isRoot = 0;
  no_assignedChildren = 0;
  std::size_t loc_comma = newickTree.find(",");
  #ifdef DEBUG
  //  std::cout <<"loc_comma = "<< loc_comma <<"\n";
#endif //DEBUG
  if(loc_comma == std::string::npos)
    {
      isTip = 1;
      popID = (unsigned) atoi(newickTree.c_str());
      age = 0;
      populationSize= runiform() * paraMax(0);  
      if(root->get_size_of_splittingTimes()<popID)
	root->resize_splittingTimes(popID);
      root->add_splittingTimes(popID-1,age);
      if(root->get_size_of_pops()< popID)
	root->resize_pops(popID);
      root->add_pops(popID-1,this);
    }
  else
    {
      isTip = 0;
      #ifdef DEBUG
      /*
      std::cout <<"par->get_age() = \n";
      std::cout << par->get_age()<<"\n";
      */
#endif //DEBUG
      if(par->get_popID() ==0) // island model and `par' is the root
	age = runiform()*paraMax(1);
      else
	age = runiform()*par->get_age();
      populationSize= runiform() * paraMax(0);
      
      std::size_t loc = newickTree.find_last_of(":");
      std::string str = newickTree.substr(loc+1,newickTree.size()-loc);
      popID = (unsigned) atoi(str.c_str());
      
      if(root->get_size_of_splittingTimes()<popID)
	root->resize_splittingTimes(popID);
      root->add_splittingTimes(popID-1,age);
      if(root->get_size_of_pops()< popID)
	root->resize_pops(popID);
      root->add_pops(popID-1,this);
      
      newickTree.erase(loc-1,str.size()+2);
      newickTree.erase(0,1);
      #ifdef DEBUG
      /*
      std::cout <<"popID = "<< popID<<"\n";
      std::cout <<"newickTree = " << newickTree <<"\n";
      */
#endif //DEBUG
      
      loc_comma = newickTree.find_first_of(",");
      std::size_t loc_openP = newickTree.find_first_of("(");
      std::size_t loc_closeP = newickTree.find_first_of(")");
      int count_open_minus_close = 0;	      
      while(newickTree.size()!=0)
	{
	  if(count_open_minus_close ==0 && loc_comma < loc_openP && loc_comma < loc_closeP)
	    {
	      loc_openP = loc_closeP; 
	    }
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
	      
	      stmp =newickTree.substr(loc_closeP+1,newickTree.size()-loc_closeP);   
	      if(stmp.find(",")==std::string::npos)
		loc_comma = newickTree.size();
	      else
		loc_comma = 1+loc_closeP + stmp.find(","); 
	    }
	  if(count_open_minus_close ==0)
	    {
	      std::string child = newickTree.substr(0, loc_comma);
	      popTree *ch = new popTree;
	      desc.push_back(ch);   
	      desc.at(no_assignedChildren)->assign_par(this);  
	      desc.at(no_assignedChildren)->initialize_popTree_recursion(im,child,paraMax,root);
	      no_assignedChildren++;
	      if(loc_comma ==std::string::npos)
		newickTree.resize(0); // loc_comma+1, newickTree.size()-loc_comma);
	      else
		newickTree.erase(0,loc_comma+1); // loc_comma+1, newickTree.size()-loc_comma);
#ifdef DEBUG
	      /*
	         std::cout <<"child = "<< child <<"\n";
		 std::cout <<"loc_comma = "<< loc_comma <<"\n";
		 std::cout <<"newickTree = "<< newickTree <<"\n";
	      */
#endif// DEBUG
	      loc_comma = newickTree.find_first_of(",");
	      loc_openP = newickTree.find_first_of("(");
	      loc_closeP = newickTree.find_first_of(")");
	      // the remaining string may have one or more children.
	      //initialize_popTree_recursion(im,newSubtr,paraMax);
	    }
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
  //  std::cout <<"In popTree::initialize_popTree(IM im, unsigned int processID)\n";
#endif //DEBUG
  std::string newickTree = im.get_poptree_string();
  unsigned int ancPop = im.get_ancPop();
  Eigen::Vector3d paraMax = im.get_paraMax();


  #ifdef DEBUG
  // std::cout <<"newickTree = "<<newickTree <<"\n";
#endif// DEBUG

  double para = 0.0;
  if(newickTree.compare(0,1,"(")==0 && newickTree.compare(newickTree.size()-1,1,";")==0)
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
	  splittingTimes.push_back(age);
	  populationSize= runiform() * paraMax(0);
	  
	  newickTree.erase(newickTree.size()-1,1);
	  popID = (unsigned) atoi(newickTree.c_str());
	  no_assignedChildren = 0;
	  
	  #ifdef DEBUG
	  /*
	  std::cout <<"Single population\n";
	  std::cout << "popID = "<<popID <<"\n";
	  */
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
	      
	      newickTree.erase(newickTree.size()-1,1);
#ifdef DEBUG
	      //  std::cout <<"newickTree = "<<newickTree <<"\n";
#endif //DEBUG
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

	      splittingTimes.resize(popID);
	      splittingTimes.at(popID-1) = age;
	      
	      pops.resize(popID);
	      pops.at(popID-1) = this;
	      
	      newickTree.erase(loc-1,str.size()+2);
	    }
	      
#ifdef DEBUG
	  //     std::cout <<"newickTree = "<<newickTree <<"\n";
#endif //DEBUG
	  no_assignedChildren = 0;  
	  
	  std::size_t loc_comma = newickTree.find(",");
	  std::size_t loc_openP = newickTree.find("(");
	  std::size_t loc_closeP = newickTree.find(")");
	  int count_open_minus_close = 0;	      
	  int loc_max =0;
	  while(newickTree.size()!=0)
	    {
	      if(count_open_minus_close ==0 && loc_comma < loc_openP && loc_comma < loc_closeP)
		{
		  loc_openP = loc_closeP; 
		}
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
		      stmp =newickTree.substr(loc_closeP+1,newickTree.size()-loc_closeP);   
		      if(stmp.find(",")==std::string::npos)
			loc_comma = newickTree.size();
		      else
			loc_comma = 1+loc_closeP + stmp.find(",");
		}
	      if(count_open_minus_close ==0)
		{ 
		  std::string child = newickTree.substr(0, loc_comma);
		  popTree *ch = new popTree;
		  desc.push_back(ch);    
		  desc.at(no_assignedChildren)->assign_par(this); 
		  desc.at(no_assignedChildren)->initialize_popTree_recursion(im,child,paraMax, this);
		  no_assignedChildren++;
		  if(loc_comma ==std::string::npos)
		    newickTree.resize(0); // loc_comma+1, newickTree.size()-loc_comma);
		  else
		    newickTree.erase(0,loc_comma+1);
		      #ifdef DEBUG
		  /*
		  std::cout <<"child = "<< child <<"\n";
		  std::cout <<"loc_comma = "<< loc_comma <<"\n";
		  std::cout <<"newickTree = "<< newickTree <<"\n";
		  */
#endif// DEBUG
		  loc_comma = newickTree.find_first_of(",");
		  loc_openP = newickTree.find_first_of("(");
		  loc_closeP = newickTree.find_first_of(")");
		  // the remaining string may have one or more children.
		  //initialize_popTree_recursion(im,newSubtr,paraMax);
		}
	    }
	}
    }   
  else
    std::cout << "The input population tree is not in newick format.\n";


  // clean up the splitting times.
  #ifdef DEBUG
  /*
  for(unsigned int i=0; i<splittingTimes.size(); i++)
    std::cout <<"i = "<< i+1 << " splttingTimes.at(i) = "<< splittingTimes.at(i) <<"\n";
  print_poptree();
  std::cout <<"pops.size() = "<< pops.size()<<"\n";
  */
#endif //DEBUG
  for(unsigned int i=0; i<splittingTimes.size()-1; i++)
    {
      if(splittingTimes.at(i) > splittingTimes.at(i+1))
	{
	  double lowerB = 0;
	  if(i != 0)
	    lowerB =  splittingTimes.at(i-1);
	  splittingTimes.at(i) = lowerB + runiform()*(splittingTimes.at(i+1)-lowerB);
	  pops.at(i)-> assign_age(splittingTimes.at(i));
	}
    }  
  #ifdef DEBUG
  /*
  for(unsigned int i=0; i<splittingTimes.size(); i++)
    std::cout <<"i = "<< i+1 << " splttingTimes.at(i) = "<< splittingTimes.at(i) <<"\n";
  print_poptree();
  */
#endif //DEBUG
  
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
	  print_migRate();
	  /*
	  std::cout << "From\tTo\tMigrationRate\n";
	  std::cout << desc.at(0)->get_popID() << "\t" << desc.at(0)->get_pop2mig(0)->get_popID()
		    <<"\t" << desc.at(0)->get_migRate(0) <<"\n";
	  std::cout << desc.at(1)->get_popID() << "\t" << desc.at(1)->get_pop2mig(0)->get_popID()
		    <<"\t" << desc.at(1)->get_migRate(0) <<"\n";
	  */
	}
      std::cout <<"\nPopulation tree\n";
      print_poptree();
      
      std::cout << "End of population tree initialization.\n\n";
    }
  
  return;
}


/*
void popTree::initialize_migrations_recursion(IM im, std::vector<popTree*> list_pops, double rateMax)
// void popTree::initialize_migrations_recursion(IM im, popTree* pop, double rateMax)
{
  // migration from `this' to pop
  pop2mig.push_back(pop);
  migRate.push_back(runiform()*rateMax);
  std::vector<double> epoch;
  epoch.push_back(max(age,pop->get_age()));
  epoch.push_back(min(par->get_age(),pop->par->get_age()));
  epoch2mig.push_back(epoch);

  // migration from pop to `this'
  pop->add_pop2mig(this);
  pop->add_migRate(runiform()*rateMax);
  pop->add_epoch2mig(epoch);
  
  
  #ifdef DEBUG
  std::cout <<"age =" << age <<" pop->get_age() = "<< pop->get_age()<<"\n";
#endif //DEBUG
  if(age < pop->get_age())
    {
      popTree* tempPop = pop;
      for(unsigned int i=0; i<tempPop->get_no_assignedChildren(); i++)
	{
	  
	  pop2mig.push_back(tempPop->get_child(i));
	  migRate.push_back(runiform()*rateMax);
	  std::vector<double> epoch;
	  epoch.push_back(max(age,tempPop->get_child(i)->get_age()));
	  epoch.push_back(min(par->get_age(),tempPop->get_age()));
	  epoch2mig.push_back(epoch);
	  epoch.resize(0);
	}
    }
  
  if(no_assignedChildren >= 2)
    {
      for(unsigned int i=0; i<no_assignedChildren; i++)
	{
	  for(unsigned int j=i+1; j<no_assignedChildren; j++)
	    {
	      // if(i != j)
		desc.at(i)->initialize_migrations_recursion(im,desc.at(j), rateMax);
	    }
	  
	}
    }
  
}
*/

void popTree::initialize_migrations(IM im, double rateMax, unsigned int processID)
{
#ifdef DEBUG
  /*
  std::cout <<"\nIn popTree::initialize_migrations()\n";
  std::cout <<"no_assignedChildren = "<< no_assignedChildren <<"\n";
  std::cout <<"pops.size() = "<< pops.size()<<"\n";
  */
#endif //DEBUG

  for(unsigned int i=0; i<pops.size(); i++)
    {
      popTree* source = pops.at(i);
      if(source->get_isRoot()==0)
	{
	  for(unsigned int j=i+1; j<pops.size(); j++)
	    {
	      popTree* receiver = pops.at(j);
	      unsigned int overlap = 0;
	      double LB, UB;
#ifdef DEBUG
	      //  std::cout <<"i = "<< i <<"j=" << j<<"\n";
#endif// DEBUG
	      if(source->get_age() == receiver->get_age())
		{
		  LB = receiver->get_age();
		  UB = min(source->par->get_age(), receiver->par->get_age());
		  overlap = 1;
		}
	      else if(source->get_age() < receiver->get_age())
		{
		  LB = receiver->get_age();
		  if(source->par->get_age() > receiver->get_age())
		    {
		      UB = min(source->par->get_age(), receiver->par->get_age());
		      overlap = 1;
		    }
		}
	      else
		{
		  LB = source->get_age();
		  if(receiver->par->get_age() > source->get_age())
		    {
		      UB = min(source->par->get_age(), receiver->par->get_age());
		      // UB = source->get_age();
		      overlap = 1;
		    }
		}
	      if(overlap)
		{
		  #ifdef DEBUG
		  // std::cout <<"i = "<< i <<"j=" << j<<"\n";
#endif// DEBUG
		  std::vector<double> epoch;
		  epoch.push_back(LB); epoch.push_back(UB);
		  source->add_pop2mig(receiver);
		  source->add_migRate(runiform()*rateMax);
		  source->add_epoch2mig(epoch);
		  receiver->add_pop2mig(source);
		  receiver->add_migRate(runiform()*rateMax);
		  receiver->add_epoch2mig(epoch);		  
		}
	    }	  
	}
    }

  /*
  if(no_assignedChildren >= 2)
    {
      for(unsigned int i=0; i<no_assignedChildren; i++)
	{
	  std::vector<popTree*> list_pops;
	  for(unsigned int j=0; j<no_assignedChildren; j++)
	    {
	      if(j != i)
		{
		  list_pops.push_back(desc.at(j));
		}
	    }
	  desc.at(i)->initialize_migrations_recursion(im,list_pops, rateMax);
	}
    }
  */
 
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
    {
      std::cout << "Population sizes:\n";
      for(unsigned int i=0; i<pops.size(); i++)
	std::cout << pops.at(i)->get_popID() << ": " << pops.at(i)->get_popSize() <<"\n";
	
    }
  else if(isTip == 0)
    {
      std::cout << popID << ": " << populationSize <<"\n";
      for(unsigned int i=0; i<no_assignedChildren; i++)
	desc.at(i)->print_popSize();
      //	  desc[1]->print_popSize();
    }
}



void popTree::print_node()
{
  double intv =0;
  if(isTip == 1)
    {
      if(par->get_popID()==0) // island model and `par' is the root
	intv = par->get_age();
      else
	intv = par->get_age()-age;
      
      std::cout << popID << ":" << intv;
    }
  else
    {
      std::cout << "(";
      for(unsigned int i=0; i<no_assignedChildren; i++)
	{
	  desc.at(i)->print_node();
	  if(i != no_assignedChildren-1)
	    std::cout << ",";
	  // desc[1]->print_node();
	}
      if(par->get_popID()==0) // island model and `par' is the root
	intv = par->get_age();
      else
	intv = par->get_age()-age;
      std::cout << ")"<< popID <<":" << intv ;
    }
}


void popTree::print_poptree()
{
	if(isTip == 1)
		std::cout << "(" << popID << ");\n";
	else
	{
	  std::cout << "(";
	  for(unsigned int i=0; i<no_assignedChildren; i++)
	    {
	      try{
		desc.at(i)->print_node();
	      }catch (std::exception &e) {
		std::cout << "In popTree::print_poptree(). Can't access element desc[0] \n";
	      }
	      if(i != no_assignedChildren-1)
		std::cout << ",";
	    }
	  /*		    
			    try{
			    desc[1]->print_node();
			    }catch (std::exception &e) {
			    std::cout << "In popTree::print_poptree(). Can't access element desc[1] \n";
			    }
	  */
	  std::cout << ")"<< popID <<";\n";
	}	
	print_popSize();
}

// print migration information
void popTree::print_migRate()
{
  if(isRoot == 1)
    {
      std::cout << "\nMigration rates\n";

      for(unsigned int i=0; i<pops.size(); i++)
	{
	  if(pops.at(i)->get_isRoot()==0)
	    {
	      for(unsigned int j=0; j<pops.at(i)->get_size_of_migRate(); j++)
		{
		  std::cout << "From pop " << pops.at(i)->get_popID()
			    << " to pop "<< pops.at(i)->get_pop2mig(j)->get_popID()
			    << ": " << pops.at(i)->get_migRate(j)
			    << " during ("<<pops.at(i)->get_epoch2mig(j,0)
			    << ","<< pops.at(i)->get_epoch2mig(j,1)<<")"
			    << "\n";
		}	      
	    }
	}
    }
  /*
  if(migRate.size()>0)
    {
      for(unsigned int i=0; i<migRate.size(); i++)
	{
	  std::cout << "From pop " << popID << " to pop "<< pop2mig.at(i)->get_popID() << ": " << migRate.at(i)
		    << " during ("<<epoch2mig.at(i).at(0)<< ","<< epoch2mig.at(i).at(1)<<")"
		    << "\n";
	}
    }
  if(no_assignedChildren >= 2)
    {
      for(unsigned int i=0; i<no_assignedChildren; i++)
	desc.at(i)->print_migRate();
    }
*/
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


