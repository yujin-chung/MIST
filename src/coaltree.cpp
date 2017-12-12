/* MIST copyright 2016 by Yujin Chung and Jody Hey */

#include "coaltree.hpp"
#include <iostream>
#include <fstream>
#include <list>





// YC 2/23/2015
// See Eq (13) from Hey and Nielsen (2007)
// When 'totalCoalRate' is defined as Eq (12) and Eq (16), then
// the correct expression of Eq (13) is
// 2^(n-1) * f^{-(n-2)} * Gamma(n-2, fn/theta_max) / theta_max, where 
// n is the number of sequences.
// double node::compute_prior_indepLoci(locus lc, double popSizeMax, double splittingTimeMax,unsigned int priorType, double changeOfRatePoint)
double node::compute_prior_indepLoci(locus lc, double popSizeMax, unsigned int priorType, double changeOfRatePoint)
{
  //std::cout << "In double node::compute_prior_indepLoci\n";

  unsigned int numCoalEvents = size_tree()-1;  
  double log_prior = 0.0;
  compute_totalCoalescentRate(); // \sum_{i=2}^n t_i*i*(i-1)

  double log_prior1,log_prior2,log_prior3,log_prior4;

  // Marginal coalescent of a single population after integrating out population size
  log_prior1 = (double) numCoalEvents *log(2.0) 
    - (double) (numCoalEvents-1) * log(totalCoalRate)
    + log_incompleteUpperGamma (numCoalEvents-1, totalCoalRate/popSizeMax)
    - log(popSizeMax);
  
  if(priorType == 0)
    {
      log_prior = 0;
    }
  else if(priorType ==1)
    {
      log_prior = log_prior1;
    }
  else if(priorType ==2)
    {
      unsigned int numCoalEvents_before =0; // the number of coalescent events before changeOfRatePoint
      unsigned int numCoalEvents_after =0; // the number of coalescent events after changeOfRatePoint
      double totalCoalRate_before =0;
      double totalCoalRate_after =0;
      unsigned int n_lineages = size_tree( );
      list<double> coalTimes;
      coalTimes = get_coalescentTimes(coalTimes);
      coalTimes.sort(); // sort the elements in ascending order
      list<double>::iterator iter;
      for(iter=coalTimes.begin(); iter!=coalTimes.end(); iter++)
	{
	  list<double>::iterator iter_prev = iter;
	  double waitingTime = 0.0;
	  if(*iter <= changeOfRatePoint)
	    {
	      if(n_lineages == size_tree())
		waitingTime = *iter;
	      else
		{	 
		  iter_prev--;
		  waitingTime = *iter - *iter_prev;
		}
	      totalCoalRate_before += waitingTime* n_lineages *(n_lineages-1);
	      numCoalEvents_before++;	     
	    }
	  else
	    {
	      if(numCoalEvents_after ==0)
		{
		  if(n_lineages == size_tree())
		    {
		      waitingTime = changeOfRatePoint;
		    }
		  else
		    { 
		      iter_prev--;		      
		      waitingTime = changeOfRatePoint-*iter_prev;
		    }
		  totalCoalRate_before += waitingTime* n_lineages *(n_lineages-1);
		  
		  waitingTime = *iter-changeOfRatePoint;
		  totalCoalRate_after += waitingTime* n_lineages *(n_lineages-1);
		  
		}
	      else
		{ 
		  iter_prev--;
		  waitingTime = *iter - *iter_prev;
		  totalCoalRate_after += waitingTime* n_lineages *(n_lineages-1);		  
		}
	      numCoalEvents_after++;
	    }
	  n_lineages--;

	  //std::cout << "*iter = " << *iter << " waitingTime = " << waitingTime <<"\n";
	}
      /*
      std::cout << "numCoalEvents_after = " << numCoalEvents_after 
		<< " numCoalEvents_before = " << numCoalEvents_before <<"\n";
      std::cout << "totalCoalRate_before = " << totalCoalRate_before
		<< "totalCoalRate_after = " << totalCoalRate_after <<"\n";
      */
	
      double log_prior_before = 0;
      double log_prior_after =0;
      if(numCoalEvents_before ==0)
	{
	  // max: the upper bound of population size
	  // f: the total coalescent rate
	  // \int_0^max exp(-f/theta) dtheta 
	  //      = max*exp(-f/max)-f*\int_{1/max}^\infty (1/x)*exp(-f*x)dx
	  log_prior_before 
	    =  log( popSizeMax * exp(-1*totalCoalRate_before/popSizeMax)
		    - totalCoalRate_before * exp(log_incompleteUpperGamma (0, totalCoalRate_before/popSizeMax)) ) - log(popSizeMax);
	}
      else
	log_prior_before = (double) numCoalEvents_before *log(2.0) 
	  - (double) (numCoalEvents_before-1) * log(totalCoalRate_before)
	  + log_incompleteUpperGamma (numCoalEvents_before-1, totalCoalRate_before/popSizeMax)
	  -log(popSizeMax);
      if(numCoalEvents_after>0)
	log_prior_after = (double) numCoalEvents_after *log(2.0) 
	  - (double) (numCoalEvents_after-1) * log(totalCoalRate_after)
	  + log_incompleteUpperGamma (numCoalEvents_after-1, totalCoalRate_after/popSizeMax)
	  -log(popSizeMax);
      log_prior = log_prior_before + log_prior_after;

   
  }
  
  return log_prior;
  
}


/***
 * Convert a string in newick form to a gene tree.
   From root
 */
void node::convertFromNewick(string tree, unsigned int root)
{
  // std::cout << "\nreading tree: "<< tree <<"\n";

  isRoot = root;

  unsigned int stringSize = tree.size();
  
  if( stringSize > 1)
    {   
      // Remove the parentheses
      std::string subTree = tree.substr(1,stringSize-2);
      if(isRoot == 1)
	subTree.resize(subTree.size()-1);
	//subTree.erase(subTree.end()-1,subTree.end());

      // std::cout << "tree = " << subTree <<"\n";

      isTip = 0;
        

      // Find the position to divide the string into two strings for firstChild and secondChild
      unsigned int found =0;      
      unsigned int found_comma = subTree.find_first_of(",");
      unsigned int found_openParenthesis = subTree.find_first_of("(");
      unsigned int count_openParenthesis = 0;
      unsigned int loc_split =0;
      string tmpTopo;
      while(found == 0)
	{
	  std::string child1;
	  std::string child2;
	  if(found_comma < found_openParenthesis)
	    {
	      child1 = subTree.substr(0,found_comma);
	      std::size_t cutoutBrLen = child1.find_last_of(":");
	      std::string age_string = child1.substr(cutoutBrLen+1);
	      age = atof(age_string.c_str());
	      child1 = child1.substr(0,cutoutBrLen);
	      child2 = subTree.substr(found_comma+1);
	      cutoutBrLen = child2.find_last_of(":");
	      child2 = child2.substr(0,cutoutBrLen);
	      desc[0] = new node;
	      desc[1] = new node;
	      
	      // std::cout << "child1 = " << child1 << " and child2 = " << child2 <<"\n";

	      desc[0]->convertFromNewick(child1, 0);
	      desc[1]->convertFromNewick(child2, 0);
	      desc[0]->par = this;
	      desc[1]->par = this;  

	      age += desc[0]->get_age();

	      found = 1;
	      //std::cout << "child1 = "<<child1 <<" and child2 = "<< child2 <<"\n";
	    }
	  else
	    {
	      count_openParenthesis++;	      
	      loc_split = found_openParenthesis+1;
	      tmpTopo = subTree.substr(found_openParenthesis+1);
	      while(count_openParenthesis >0)
		{
		  // std::cout << "tmpTopo =" << tmpTopo <<"and count_openParenthesis = "<< count_openParenthesis << "\n";
		  unsigned int found_closeParenthesis = tmpTopo.find_first_of(")");
		  found_openParenthesis = tmpTopo.find_first_of("(");
		  if(found_closeParenthesis < found_openParenthesis)
		    {
		      tmpTopo = tmpTopo.substr(found_closeParenthesis+1);
		      loc_split += found_closeParenthesis+1;
		      count_openParenthesis--;
		    }
		  else
		    {
		      count_openParenthesis++;
		      tmpTopo = tmpTopo.substr(found_openParenthesis+1);
		      loc_split += found_openParenthesis+1;
		    }
		}	      
	      found_comma = loc_split + tmpTopo.find_first_of(",");
	      found_openParenthesis =  found_comma +1;
	    }
	}
    }
  else
    { 
      age = 0;    
      isTip = 1;
      // size = 1;   
      tipID = (unsigned int)  atoi(tree.c_str());
    }

  // REMOVE
  // std::cout << "reading tree: " << tree << " and saving as \n";
  // print_coaltree();
  
  return;
}

/*
void nodeSimple::assignNodeLabel()
{
  if(isTip)
    nodeLabel = popID;
  else
    {
      firstChild->assignNodeLabel();
      secondChild->assignNodeLabel();
      int label1 = firstChild->get_nodeLabel();
      int label2 = secondChild->get_nodeLabel();
      if(label1==label2)
	nodeLabel = label1;
      else
	nodeLabel = -1;
    }
  std::cout <<"nodeLabel = "<< nodeLabel <<"\n";
  return;
}
*/


/***
 * Convert a string in newick form to a ranked topology.
   From root
 */
void nodeSimple::convertFromNewick(string topo, unsigned int root)
{
  isRoot = root;

  unsigned int stringSize = topo.size();
 
  if(stringSize > 1)
    {
      // Extract the rank
      unsigned int found = topo.find_last_of(":");
      std::string subTopo = topo.substr(0,found);
      std::string rank_string = topo.substr(found+1);

      isTip = 0;
      rank = (unsigned) atoi(rank_string.c_str());
      popID =0;
      
      // Remove the parentheses
      subTopo.erase(subTopo.end()-1,subTopo.end());
      subTopo.erase(subTopo.begin(),subTopo.begin()+1);

      // Find the position to divide the string into two strings for firstChild and secondChild
      found =0;      
      unsigned int found_comma = subTopo.find_first_of(",");
      unsigned int found_openParenthesis = subTopo.find_first_of("(");
      unsigned int count_openParenthesis = 0;
      unsigned int loc_split =0;
      string tmpTopo;
      while(found == 0)
	{
	  std::string child1;
	  std::string child2;
	  if(found_comma < found_openParenthesis)
	    {
	      child1 = subTopo.substr(0,found_comma);
	      child2 = subTopo.substr(found_comma+1);
	      firstChild = new nodeSimple;
	      secondChild = new nodeSimple;
	      firstChild->convertFromNewick(child1, 0);
	      secondChild->convertFromNewick(child2, 0);
	      firstChild->par = this;
	      secondChild->par = this;
	      
	      found = 1;
	      //std::cout << "child1 = "<<child1 <<" and child2 = "<< child2 <<"\n";
	    }
	  else
	    {
	      count_openParenthesis++;	      
	      loc_split = found_openParenthesis+1;
	      tmpTopo = subTopo.substr(found_openParenthesis+1);
	      while(count_openParenthesis >0)
		{
		  // std::cout << "tmpTopo =" << tmpTopo <<"and count_openParenthesis = "<< count_openParenthesis << "\n";
		  unsigned int found_closeParenthesis = tmpTopo.find_first_of(")");
		  found_openParenthesis = tmpTopo.find_first_of("(");
		  if(found_closeParenthesis < found_openParenthesis)
		    {
		      tmpTopo = tmpTopo.substr(found_closeParenthesis+1);
		      loc_split += found_closeParenthesis+1;
		      count_openParenthesis--;
		    }
		  else
		    {
		      count_openParenthesis++;
		      tmpTopo = tmpTopo.substr(found_openParenthesis+1);
		      loc_split += found_openParenthesis+1;
		    }
		}	      
	      found_comma = loc_split + tmpTopo.find_first_of(",");
	      found_openParenthesis =  found_comma +1;
	    }
	}
    }
  else
    {
      isTip = 1;
      size = 1;      
      rank =0;
      popID = (unsigned int)  atoi(topo.c_str());
    }
  return;
}



// For class node
/***
 * Convert a coalescent tree to a string in newick form.
 */
std::string node::convert2string_newick_root()
{
	std::string tree;
	if(isTip == 1)
	{
		std::stringstream tmp;
		tmp << popID;
		tree = tmp.str();
	}
	else
	{
		tree = "(";
		tree += desc[0]->convert2string_newick();
		tree += ",";
		tree += desc[1]->convert2string_newick();
		tree += ")";
	}
	return tree;
}
// For class node
/***
 * Convert a ranked topology to a string in newick form.
 */
std::string node::convert2string_newick()
{
	std::string tree;
	if(isTip == 1)
	{
		std::stringstream tmp;
		tmp << popID;
		tree = tmp.str();
	}
	else
	{
		tree = "(";
		tree += desc[0]->convert2string_newick();
		tree += ",";
		tree += desc[1]->convert2string_newick();
		tree += "):";

		std::stringstream tmp;
		tmp << par->age-age;
		tree += tmp.str();
	}
	return tree;
}



// For class nodeSimple
/***
 * Convert a ranked topology to a string in newick form.
 */
std::string nodeSimple::convert2Newick_topo()
{
	std::string tree;
	if(isTip == 1)
	{
		std::stringstream tmp;
		tmp << popID;
		tree = tmp.str();
	}
	else
	{
		tree = "(";
		tree += firstChild->convert2Newick_topo();
		tree += ",";
		tree += secondChild->convert2Newick_topo();
		tree += "):";

		std::stringstream tmp;
		tmp << rank;
		tree += tmp.str();
	}
	return tree;
}


/***
 * Convert a ranked topology to a string in newick form.
 */
std::string nodeSimple::convert2Newick_topo_root()
{
	std::string tree;
	if(isTip == 1)
	{
		std::stringstream tmp;
		tmp << popID;
		tree = tmp.str();
	}
	else
	{
		tree = "(";
		tree += firstChild->convert2Newick_topo();
		tree += ",";
		tree += secondChild->convert2Newick_topo();
		tree += ")";
	}
	return tree;
}

std::string nodeSimple::convert2Newick_recursion(std::list<double> coalT, std::vector<unsigned int> tipIDs, double parentAge)
{
	std::string tree;
	if(isTip == 1)
	{
		std::stringstream tmp;
		tmp << tipIDs.at(0) << ":" << parentAge ;
		tree = tmp.str();
		if(tipIDs.size() > 1)
			tipIDs.erase(tipIDs.begin());
	}
	else
	{
		double age = coalT.front();

		tree = "(";

		std::list<double> new_coalT = coalT;
		new_coalT.pop_front();
		tree += firstChild->convert2Newick_recursion(new_coalT,tipIDs,age);

		tree += ",";

		if(firstChild->size >= 2)
		{
			std::list<double>::iterator it = new_coalT.begin();
			std::advance(it,firstChild->size - 1);
			//it +=(firstChild->size - 1);
			new_coalT.erase(new_coalT.begin(),it);
		}

		std::vector<unsigned int> new_tipIDs;
		std::vector<unsigned int>::iterator it_tips = tipIDs.begin();
		std::advance(it_tips,firstChild->size);
		new_tipIDs.assign(it_tips,tipIDs.end());

		tree += secondChild->convert2Newick_recursion(new_coalT,new_tipIDs,age);

		tree += "):";

		std::stringstream tmp;
		tmp << parentAge - age;
		tree += tmp.str();
	}

	return tree;
}

std::string nodeSimple::convert2Newick(std::list<double> coalT, std::vector<unsigned int> tipIDs)
{
	std::string tree;
	if(isTip == 1)
	{
		std::stringstream tmp;
		tmp << "(" << tipIDs.at(0) << ");\n" ;
		tree = tmp.str();
		// tipIDs.erase(tipIDs.begin());
	}
	else
	{
		tree = "(";

		std::list<double> new_coalT = coalT;
		new_coalT.pop_front();
		tree += firstChild->convert2Newick_recursion(new_coalT,tipIDs,coalT.front());

		tree += ",";

		if(firstChild->size >= 2)
		{
			std::list<double>::iterator it = new_coalT.begin();
			std::advance(it,firstChild->size - 1);
			new_coalT.erase(new_coalT.begin(),it);
		}

		std::vector<unsigned int> new_tipIDs;
		std::vector<unsigned int>::iterator it_tips = tipIDs.begin();
		std::advance(it_tips,firstChild->size);
		new_tipIDs.assign(it_tips,tipIDs.end());

		tree += secondChild->convert2Newick_recursion(new_coalT,new_tipIDs,coalT.front());

		tree += ");\n";
	}

	return tree;
}



std::string nodeSimple::convert2Newick_originalTopo(std::vector<unsigned int> tipIDs)
{
  std::string tree;
  if(isTip == 1)
    {
      std::stringstream tmp;
      tmp << "(" << tipIDs.at(0) << ");\n" ;
      tree = tmp.str();
    }
  else
    {
      tree = "(";
      
      //std::list<double> new_coalT = coalT;
      //new_coalT.pop_front();
      tree += firstChild->convert2Newick_originalTopo_recursion(tipIDs);
      
      tree += ",";
      
      /*
      if(firstChild->size >= 2)
	{
	  std::list<double>::iterator it = new_coalT.begin();
	  std::advance(it,firstChild->size - 1);
	  new_coalT.erase(new_coalT.begin(),it);
	}
      */

      std::vector<unsigned int> new_tipIDs;
      std::vector<unsigned int>::iterator it_tips = tipIDs.begin();
      std::advance(it_tips,firstChild->size);
      new_tipIDs.assign(it_tips,tipIDs.end());
      
      tree += secondChild->convert2Newick_originalTopo_recursion(new_tipIDs);
      // tree += secondChild->convert2Newick_originalTopo_recursion(new_coalT,new_tipIDs,coalT.front());
      
      //tree += ");\n";
      tree += "):";
      std::stringstream tmp;
      tmp << rank;
      tree += tmp.str();
      tree += ";\n";
    }
  
  return tree;
}


// std::string nodeSimple::convert2Newick_originalTopo_recursion(std::list<double> coalT, std::vector<unsigned int> tipIDs, double parentAge)
std::string nodeSimple::convert2Newick_originalTopo_recursion(std::vector<unsigned int> tipIDs)
{
	std::string tree;
	if(isTip == 1)
	{
		std::stringstream tmp;
		tmp << tipIDs.at(0); // << ":" << parentAge ;
		tree = tmp.str();
		if(tipIDs.size() > 1)
			tipIDs.erase(tipIDs.begin());
	}
	else
	{
	  // double age = coalT.front();

	  tree = "(";
	  
	  // std::list<double> new_coalT = coalT;
	  // new_coalT.pop_front();
	  // tree += firstChild->convert2Newick_originalTopo_recursion(new_coalT,tipIDs,age);
	  tree += firstChild->convert2Newick_originalTopo_recursion(tipIDs);
	  
	  tree += ",";
	  
	  /*
	    if(firstChild->size >= 2)
	    {
	    std::list<double>::iterator it = new_coalT.begin();
	    std::advance(it,firstChild->size - 1);
	    //it +=(firstChild->size - 1);
	    new_coalT.erase(new_coalT.begin(),it);
	    }
	  */

	  std::vector<unsigned int> new_tipIDs;
	  std::vector<unsigned int>::iterator it_tips = tipIDs.begin();
	  std::advance(it_tips,firstChild->size);
	  new_tipIDs.assign(it_tips,tipIDs.end());

	  // tree += secondChild->convert2Newick_recursion(new_coalT,new_tipIDs,age);
	  tree += secondChild->convert2Newick_originalTopo_recursion(new_tipIDs);
	  
	  //tree += ")";
	  tree += "):";
	  
	  std::stringstream tmp;
	  tmp <<rank;
	  // tmp << parentAge - age;
	  tree += tmp.str();
	}
	
	return tree;
}



// Updated by YC 5/8/2014
// Argument 'nPops' is added.
/***
 * @param coaltimes list of coalescent times from a tree, in ascending order.
 * @param nLineages the ranks for internal nodes are larger than nLineages.
 */
void nodeSimple::convert(node* tree, std::list<double> coaltimes, unsigned int nLineages)
{
  isRoot = tree->isRoot;
  isTip  = tree->isTip;
  
  unsigned int foundSameCoalTime = 0;
  std::list<double>::iterator iter = coaltimes.begin();
  unsigned int count =nLineages;
  while(foundSameCoalTime == 0 && count < coaltimes.size()+nLineages)
    {
      if(tree->age == *iter)
	{
	  foundSameCoalTime = 1;
	}
      count++;
      iter++;
    }
  if(foundSameCoalTime == 1)
    rank = count;
  
  if(isRoot == 1)
    {
      par =0;
    }
  if(isTip == 1)
    {
      popID = tree->popID;
      rank = 0;
    }
  else
    {
      firstChild = new nodeSimple;
      secondChild = new nodeSimple;
      if(tree->desc[0]->age < tree->desc[1]->age)
	{
	  firstChild->convert(tree->desc[0],coaltimes,nLineages);
	  secondChild->convert(tree->desc[1],coaltimes,nLineages);
	}
      else if(tree->desc[0]->age > tree->desc[1]->age)
	{
	  firstChild->convert(tree->desc[1],coaltimes,nLineages);
	  secondChild->convert(tree->desc[0],coaltimes,nLineages);
	}
      else if(tree->desc[0]->get_isTip() == 1&&  tree->desc[1]->get_isTip()==1)
	{
	  if(tree->desc[0]->get_tipID() < tree->desc[1]->get_tipID())
	    {
	      firstChild->convert(tree->desc[0],coaltimes,nLineages);
	      secondChild->convert(tree->desc[1],coaltimes,nLineages);
	    }
	  else
	    {
	      firstChild->convert(tree->desc[1],coaltimes,nLineages);
	      secondChild->convert(tree->desc[0],coaltimes,nLineages);
	    }
	}
      firstChild->par = this;
      secondChild->par = this;
    }
  return;
}



void node::MPIsend_coaltimes_tipIDs()
{
  if (size_tree() > 1) 
    {
      MPI::COMM_WORLD.Send(&age, 1, MPI::DOUBLE, 0, 298);
    }
  if (isTip == 0) 
    {
      if (desc[0]->age < desc[1]->age) 
	{
	  desc[0]->MPIsend_coaltimes_tipIDs();
	  desc[1]->MPIsend_coaltimes_tipIDs();
	}
      else if (desc[0]->age > desc[1]->age)
	{
	  desc[1]->MPIsend_coaltimes_tipIDs();
	  desc[0]->MPIsend_coaltimes_tipIDs();
	}
      else if (desc[0]->get_isTip() ==1 && desc[1]->get_isTip() ==1)
	{
	  if(desc[0]->get_tipID() < desc[1]->get_tipID())
	    {
	      desc[0]->MPIsend_coaltimes_tipIDs();
	      desc[1]->MPIsend_coaltimes_tipIDs();
	    }
	  else if(desc[0]->get_tipID() > desc[1]->get_tipID())
	    {
	      desc[1]->MPIsend_coaltimes_tipIDs();
	      desc[0]->MPIsend_coaltimes_tipIDs();
	    }
	  else
	    {
	      std::cout << "\n*** Error in node::MPIsend_coaltimes_tipIDs() ***\n";
	      std::cout << "(sub)tree is ";
	      print_coaltree();
	      std::cout << "desc[0]->age = " << desc[0]->age
			<< "desc[1]->age = " << desc[1]->age
			<< "desc[0]->get_isTip() = " << desc[0]->get_isTip()
			<< "desc[1]->get_isTip() = " << desc[1]->get_isTip()
			<< "desc[0]->get_tipID() = " << desc[0]->get_tipID()
			<< "desc[1]->get_tipID() = " << desc[1]->get_tipID()
			<<"\n";	      
	    }
	}
      else
	{
	  std::cout << "\n*** Error in node::MPIsend_coaltimes_tipIDs() ***\n";
	  std::cout << "(sub)tree is ";
	  print_coaltree();
	  std::cout << "desc[0]->age = " << desc[0]->age
		    << "desc[1]->age = " << desc[1]->age
		    << "desc[0]->get_isTip() = " << desc[0]->get_isTip()
		    << "desc[1]->get_isTip() = " << desc[1]->get_isTip()
		    <<"\n";
	}
    } 
  else 
    {
      MPI::COMM_WORLD.Send(&tipID, 1, MPI::UNSIGNED, 0, 211);
    }
  return;
}

/*
void node::send_tipids()
{
	if (isTip == 1) {
		MPI::COMM_WORLD.Send(&tipID, 1, MPI::UNSIGNED, 0, 370);
		//	MPI::COMM_WORLD.Send(&popID, 1, MPI::UNSIGNED, 0, 3);
	} else {
		desc[0]->send_tipids();
		desc[1]->send_tipids();
	}
	return;
}
*/

void nodeSimple::MPIreceive_tipids(int coldprocess)
{
	MPI::Status status;
	if (isTip == 1) {
		MPI::COMM_WORLD.Recv(&popID, 1, MPI::UNSIGNED, coldprocess, 370, status);
	} else {
		firstChild->MPIreceive_tipids(coldprocess);
		secondChild->MPIreceive_tipids(coldprocess);
	}
	return;
}

void node::send_size()
{
	unsigned int treesize = size_tree();
	std::cout << "Sending treesize\n";
	MPI::COMM_WORLD.Send(&treesize, 1, MPI::UNSIGNED, 0, 2000);
	//MPI::COMM_WORLD.Send(&treesize, 1, MPI::UNSIGNED, 0, 2);
	return;
}

unsigned int node::receive_size()
{
	unsigned int toposize = 0;
	MPI::Status status;
	std::cout << "Receiving toposize\n";
	MPI::COMM_WORLD.Recv(&toposize, 1, MPI::UNSIGNED, 0, 5000, status);
	//MPI::COMM_WORLD.Recv(&toposize, 1, MPI::UNSIGNED, 0, 5, status);
	return toposize;
}


void nodeSimple::send_size(int coldprocess)
{
	unsigned int toposize = size;
	std::cout << "Sending toposize\n";
	MPI::COMM_WORLD.Send(&toposize, 1, MPI::UNSIGNED, coldprocess, 5000);
	//	MPI::COMM_WORLD.Send(&toposize, 1, MPI::UNSIGNED, coldprocess, 5);
	return;
}
	

unsigned int nodeSimple::receive_size(int coldprocess)
{
	unsigned int treesize = 0;
	MPI::Status status;
	std::cout << "Receiving treesize\n";
	MPI::COMM_WORLD.Recv(&treesize, 1, MPI::UNSIGNED, coldprocess, 2000, status);
	//MPI::COMM_WORLD.Recv(&treesize, 1, MPI::UNSIGNED, coldprocess, 2, status);
	return treesize;
}


//AS: Adding MPI version
//This has to happen simultaneously between the head node, which is saving the nodesimple stuff
//and the cold chain node. 
void node::MPIsend_coaltree()
{
	//std::cout << "Sending isRoot\n";
  MPI::COMM_WORLD.Send(&isRoot, 1, MPI::UNSIGNED, 0, 190);
  //std::cout << "Sending isTip\n";
  MPI::COMM_WORLD.Send(&isTip, 1, MPI::UNSIGNED, 0, 191);
  //std::cout << "Sending age\n";
  MPI::COMM_WORLD.Send(&age, 1, MPI::DOUBLE, 0, 192);
  //	MPI::COMM_WORLD.Send(&isRoot, 1, MPI::UNSIGNED, 0, 0);
  //std::cout << "Sending isTip\n";
  //	MPI::COMM_WORLD.Send(&isTip, 1, MPI::UNSIGNED, 0, 1);
  //std::cout << "Sending age\n";
  //	MPI::COMM_WORLD.Send(&age, 1, MPI::DOUBLE, 0, 2);
  if (isTip == 0) 
    {
      if(desc[0]->isTip ==1 && desc[1]->isTip == 1)
	{	
	  if(desc[0]->popID <= desc[1]->popID)
	    {
	      desc[0]->MPIsend_coaltree();
	      desc[1]->MPIsend_coaltree();
	    } 
	  else 
	    {
	      desc[1]->MPIsend_coaltree();
	      desc[0]->MPIsend_coaltree();
	    }
	}
      else
	{
	  if (desc[0]->age <= desc[1]->age) 
	    {
	      desc[0]->MPIsend_coaltree();
	      desc[1]->MPIsend_coaltree();
	    } 
	  else 
	    {
	      desc[1]->MPIsend_coaltree();
	      desc[0]->MPIsend_coaltree();
	    }
	}
    } 
  else 
    {
      //std::cout << "Sending popID\n";
      MPI::COMM_WORLD.Send(&popID, 1, MPI::UNSIGNED, 0, 193);
      //	MPI::COMM_WORLD.Send(&popID, 1, MPI::UNSIGNED, 0, 3);
    }
  return;
}

void node::send_topoinfo()
{
	send_size();
	unsigned int toposize = receive_size();
	unsigned int treesize = size_tree();
	unsigned int foundSameTopo = 0;
	MPI::Status status;
	if (toposize == treesize) {
		if (isTip == 1) {
			std::cout << "Sending popID inside send_topoinfo\n";
			MPI::COMM_WORLD.Send(&popID, 1, MPI::UNSIGNED, 0, 400);
			//MPI::COMM_WORLD.Send(&popID, 1, MPI::UNSIGNED, 0, 4);
		} else if (treesize == 2) {
			std::cout << "Sending desc0popid\n";
			MPI::COMM_WORLD.Send(&desc[0]->popID, 1, MPI::UNSIGNED, 0, 600);
			//	MPI::COMM_WORLD.Send(&desc[0]->popID, 1, MPI::UNSIGNED, 0, 6);
			std::cout << "Sending desc1popid\n";
			MPI::COMM_WORLD.Send(&desc[1]->popID, 1, MPI::UNSIGNED, 0, 700);
			//MPI::COMM_WORLD.Send(&desc[1]->popID, 1, MPI::UNSIGNED, 0, 7);
		} else {
			std::cout << "Sending desc0age\n";
			MPI::COMM_WORLD.Send(&desc[0]->age, 1, MPI::DOUBLE, 0, 1010);
			//MPI::COMM_WORLD.Send(&desc[0]->age, 1, MPI::DOUBLE, 0, 10);
			std::cout << "Sending desc1age\n";
			MPI::COMM_WORLD.Send(&desc[1]->age, 1, MPI::DOUBLE, 0, 1111);
			//MPI::COMM_WORLD.Send(&desc[1]->age, 1, MPI::DOUBLE, 0, 11);
			if (desc[0]->age <= desc[1]->age) {
				desc[0]->send_topoinfo();
				std::cout << "Receiving foundsametopo inside send_topoinfo\n";
				MPI::COMM_WORLD.Recv(&foundSameTopo, 1, MPI::UNSIGNED, 0, 800, status);
				//	MPI::COMM_WORLD.Recv(&foundSameTopo, 1, MPI::UNSIGNED, 0, 8, status);
				if(foundSameTopo == 1) {
					desc[1]->send_topoinfo();
				}
			} else {
				desc[1]->send_topoinfo();
				std::cout << "Receinving foudnsametopon inside send_tpopinfo\n";
				MPI::COMM_WORLD.Recv(&foundSameTopo, 1, MPI::UNSIGNED, 0, 900, status);
				//	MPI::COMM_WORLD.Recv(&foundSameTopo, 1, MPI::UNSIGNED, 0, 9, status);
				if (foundSameTopo == 1) {
					desc[0]->send_topoinfo();
				}
			}

		}
	} else {
	std::cout << "returning nothing inside send_topoinfo...\n";
	}
	return;
}

unsigned int nodeSimple::receive_topoinfo(int coldprocess)
{
	MPI::Status status;
	unsigned int treesize = receive_size(coldprocess);
	send_size(coldprocess);
	unsigned int foundSameTopo = 0;
	unsigned int desc0pid = 0;
	unsigned int desc1pid = 0;
	double desc0age = 0.0;
	double desc1age = 0.0;
	unsigned int treepid = 0;
	if (size != treesize) {
		std::cout << "FoundSameTopo inside receive topoinfo is 0\n";
		foundSameTopo = 0;
	} else {
		if (isTip == 1) {
			std::cout << "Receiving treepid\n";
			MPI::COMM_WORLD.Recv(&treepid, 1, MPI::UNSIGNED, coldprocess, 400, status);
			//MPI::COMM_WORLD.Recv(&treepid, 1, MPI::UNSIGNED, coldprocess, 4, status);
			if (treepid != popID) {
				foundSameTopo = 0;
			}
			else if (size == 2) {
				std::cout << "Receiving desc0pid\n";
				MPI::COMM_WORLD.Recv(&desc0pid, 1, MPI::UNSIGNED, coldprocess, 600, status);
				//	MPI::COMM_WORLD.Recv(&desc0pid, 1, MPI::UNSIGNED, coldprocess, 6, status);
				std::cout << "Receiving desc1pid\n";
				MPI::COMM_WORLD.Recv(&desc1pid, 1, MPI::UNSIGNED, coldprocess, 700, status);
				//MPI::COMM_WORLD.Recv(&desc1pid, 1, MPI::UNSIGNED, coldprocess, 7, status);
				if ((firstChild->getPopID() == desc0pid && secondChild->getPopID() == desc1pid)
					|| (firstChild->getPopID() == desc1pid && secondChild->getPopID() == desc0pid)) {
					foundSameTopo = 1;
				} else {
					foundSameTopo = 0;
				}
			} else {
				std::cout << "Receiving desc0age\n";
				MPI::COMM_WORLD.Recv(&desc0age, 1, MPI::DOUBLE, coldprocess, 1010, status);
				//	MPI::COMM_WORLD.Recv(&desc0age, 1, MPI::DOUBLE, coldprocess, 10, status);
				std::cout << "Receving desc1age\n";
				MPI::COMM_WORLD.Recv(&desc1age, 1, MPI::DOUBLE, coldprocess, 1111, status);
				//MPI::COMM_WORLD.Recv(&desc1age, 1, MPI::DOUBLE, coldprocess, 11, status);
				if (desc0age <= desc1age) {
					foundSameTopo = firstChild->receive_topoinfo(coldprocess);
					std::cout << "Sending foundSameTopo inside receive_topoinfo\n";
					MPI::COMM_WORLD.Send(&foundSameTopo, 1, MPI::UNSIGNED, coldprocess, 800);
					//MPI::COMM_WORLD.Send(&foundSameTopo, 1, MPI::UNSIGNED, coldprocess, 8);
					if (foundSameTopo == 1) {
						foundSameTopo = secondChild->receive_topoinfo(coldprocess);
					}
				} else {
					foundSameTopo = firstChild->receive_topoinfo(coldprocess);
					std::cout << "Sending foundSameTopo inside receive_topoinfo\n";
					MPI::COMM_WORLD.Send(&foundSameTopo, 1, MPI::UNSIGNED, coldprocess, 900);
					//MPI::COMM_WORLD.Send(&foundSameTopo, 1, MPI::UNSIGNED, coldprocess, 9);
					if (foundSameTopo == 1) {
						foundSameTopo = secondChild->receive_topoinfo(coldprocess);
					}
				}
			}				
		}
	}
	std::cout << "Returning foundsampetopo = " << foundSameTopo << "\n";
	return foundSameTopo;
}



void nodeSimple::MPIreceive_coaltree(std::list<double> ct, int coldprocess)
{	
  MPI::Status status;
  //std::cout << "Receiving isroot\n";
  MPI::COMM_WORLD.Recv(&isRoot, 1, MPI::UNSIGNED, coldprocess, 190, status);
  //std::cout << "Receiving istip\n";
  MPI::COMM_WORLD.Recv(&isTip, 1, MPI::UNSIGNED, coldprocess, 191, status);
  double age = 0.0;
  unsigned int foundSameCoalTime = 0;
  unsigned int nLineages = ct.size()+1;
  unsigned int count = nLineages;
  //std::cout << "Receiving age\n";
  MPI::COMM_WORLD.Recv(&age, 1, MPI::DOUBLE, coldprocess, 192, status);
  std::list<double>::iterator iter = ct.begin();
  while (foundSameCoalTime == 0 && count < ct.size()+nLineages) {
    if (age == *iter) 
      {
	foundSameCoalTime = 1;
      }
    count++;
    iter++;
  }
  if (foundSameCoalTime == 1) 
    {
      rank = count;
    }
  if (isRoot == 1) 
    {
      par = 0;
    }
  if (isTip == 1) 
    {
      //std::cout << "Receiving popID\n";
      MPI::COMM_WORLD.Recv(&popID, 1, MPI::UNSIGNED, coldprocess, 193, status);
    }
  if (isTip == 0) 
    {
      firstChild = new nodeSimple;
      secondChild = new nodeSimple;
      firstChild->MPIreceive_coaltree(ct, coldprocess);
      secondChild->MPIreceive_coaltree(ct, coldprocess);
      firstChild->par = this;
      secondChild->par = this;
    }
  return;
}


unsigned int nodeSimple::sameTopo(nodeSimple* topo)
{
	unsigned int same = 1;

	if(size != topo->size)
	{
		same = 0;
	}
	else if(isTip==1)
	{
		if(popID != topo->popID)
		{	
			same = 0;
		}
	}
	else if(rank != topo->rank)
	{
		same =0;
	}
	else 
	{
			same = firstChild->sameTopo(topo->firstChild);
			if(same ==1)
				same = secondChild->sameTopo(topo->secondChild);	
	}
	return same;
}





void nodeSimple::computeSizes()
{
	unsigned int count = 0;
	//if(size > 0)
	//	count = size;
	//else
	//{
		if(isTip == 0)
		{
			firstChild->computeSizes();
			count += firstChild->getSize();
			secondChild->computeSizes();
			count +=secondChild->getSize();
		}
		else
		{
			count++;
		}
		size = count;
	//}
	return;
}

unsigned int nodeSimple::sameTopo(node* tree)
{
  unsigned int same = 1;
  
  //std::cout << "nodeSimple::sameTopo().\t size = " << size <<"\n";
  //std::cout << "tree->size_tree() = "<<tree->size_tree() << "\n";
  //std::cout << "Comparing..";
  //tree->print_coaltree();
  //print_topo();
  //std::cout << "\n";
  
  //if(size <= 0)
  //		computeSizes();
  
  
  if(size != tree->size_tree())
    {
      same = 0;
    }
  else
    {
      if(isTip==1)
	{
	  if(popID != tree->popID)
	    same = 0;
	}
      else if(size == 2)
	{
	  if((firstChild->getPopID() == tree->desc[0]->popID && secondChild->getPopID() == tree->desc[1]->popID)
	     || (firstChild->getPopID() == tree->desc[1]->popID && secondChild->getPopID() == tree->desc[0]->popID))
	    same = 1;
	  else
	    same = 0;
	}
      else
	{
	  
	  if(tree->desc[0]->age < tree->desc[1]->age)
	    {
	      same = firstChild->sameTopo(tree->desc[0]);
	      if(same == 1)
		same = secondChild->sameTopo(tree->desc[1]);
	    }
	  else if(tree->desc[0]->age > tree->desc[1]->age)
	    {
	      same = firstChild->sameTopo(tree->desc[1]);
	      if(same == 1)
		same = secondChild->sameTopo(tree->desc[0]);
	    }
	  //}
	}
    }
  return same;
}

void saveCoalTimes2File(vector<vector<node* > > trees)
{
  ofstream treefile;
  treefile.open ("coalescentTimes_treeSample.txt");
  treefile << "iter\tlocus\tcoalescent times\n";
  for(unsigned int i=0; i<trees.size(); i++)
    for(unsigned int j=0; j<trees.at(i).size(); j++)
      {
	treefile << i <<"\t" <<j <<"\t";
	list<double> coaltime;
	coaltime.resize(0);
	coaltime = trees.at(i).at(j)->get_coalescentTimes(coaltime);
	coaltime.sort();
	for(list<double>::iterator it= coaltime.begin() ;it != coaltime.end(); ++it)
	  {
	    treefile << *it <<" ";
	  }
	treefile <<"\n";
      }
  treefile.close();
}


void saveTrees2File(vector<vector<node* > > trees)
{
	ofstream treefile;
	treefile.open ("treeSample.txt");
	treefile << "iter\tlocus\tTree\n";
	for(unsigned int i=0; i<trees.size(); i++)
		for(unsigned int j=0; j<trees.at(i).size(); j++)
		{
			treefile << i <<"\t" <<j <<"\t";
			trees.at(i).at(j)->saveTree(treefile);
		}
	treefile.close();
}


void node::saveNode(ofstream& filename)
{
	if(isTip == 1)
		filename << tipID << ":" << par->age;
	else if(isTip == 0)
	{
		filename << "(";
		desc[0]->saveNode(filename);
		filename << ",";
		desc[1]->saveNode(filename);
		filename << "):" << par->age - age ;
	}
}

void node::saveTree(ofstream& filename)
{
	if(isTip == 1)
		filename << "(" << tipID << ");\n";
	else if(isTip == 0)
	{
		filename << "(";
		desc[0]->saveNode(filename);
		filename << ",";
		desc[1]->saveNode(filename);
		filename << ");\n";
	}
}

node::node()
{
	label = 0;
	totalCoalRate = 0.0;
	isRoot = -1;
	isTip = -1;
	isLikelihoodNULL = 0;
	tipID = 0;
	age = 0;
	popID = 0;
	// valarray<double> val(1);
	lik.resize(1,1);
	lik.setZero();
	//lik.initialize(1,1);
	n_mut = 0;
	par = 0; // null constant pointer
	desc[0] =0; desc[1] =0;
	siblingOrder = 2;
}

node_old::node_old()
{
	label = 0;
	totalCoalRate = 0.0;
	isRoot = -1;
	isTip = -1;
	isLikelihoodNULL = 0;
	tipID = 0;
	age = 0;
	popID = 0;
	// valarray<double> val(1);
	lik.initialize(1,1);
	n_mut = 0;
	par = 0; // null constant pointer
	desc[0] =0; desc[1] =0;
	siblingOrder = 2;
}



void node::initialization()
{
  label = 0;
  totalCoalRate = 0.0;
  isRoot = -1;
  isTip = -1;
  isLikelihoodNULL = 0;
  tipID = 0;
  age = 0;
  popID = 0;
  lik.resize(1,1);
  lik.setZero();
  //lik.initialize(0.0,1,1);
  n_mut = 0;
  par = 0; // null constant pointer
  siblingOrder = 2;
  desc[0] =0; desc[1] =0;
  
  totalLengthsLineages = 0;
}

void node_old::initialization()
{

	label = 0;
	totalCoalRate = 0.0;
	isRoot = -1;
	isTip = -1;
	isLikelihoodNULL = 0;
	tipID = 0;
	age = 0;
	popID = 0;
	lik.initialize(0.0,1,1);
	n_mut = 0;
	par = 0; // null constant pointer
	siblingOrder = 2;
	desc[0] =0; desc[1] =0;
}



void node::deleteCoalTree()
{
	if(isTip == 0 )
	{
		desc[0]->deleteCoalTree();
		desc[1]->deleteCoalTree();
	}
	delete this;
	//delete desc[0];
	//delete desc[1];
}


void node::print_labelsPopIDs()
{
	if(age == 0.0)
		std::cout << "tipID: " << tipID
		<< ", label: " << label
		<<", popID: " << popID <<"\n";
	else if(isTip ==0)
	{
		desc[0]->print_labelsPopIDs();
		desc[1]->print_labelsPopIDs();
	}
}

int node::validTree()
{
	// REMOVE
	//std::cout << "in node::validTree()\n The given tree is";
	//print_coaltree();

	int isValid = 1;
	if(isTip == 1) // tip
	{
		if(age !=0)
		{
			isValid = 0;
			std::cout << "In node::validTree()\n"
					"The tip's age is not 0, but" << age <<"\n";
		}
	}
	else
	{
		if(desc[0]->get_siblingOrder() != 0)
		{
			isValid = 0;
			std::cout << "In node::validTree()\n"
					"desc[0]'s siblingOrder is not 0, but" << desc[0]->get_siblingOrder() <<"\n";
		}
		if(desc[1]->get_siblingOrder() != 1)
		{
			isValid = 0;
			std::cout << "In node::validTree()\n"
					"desc[1]'s siblingOrder is not 0, but" << desc[1]->get_siblingOrder() <<"\n";
		}
		isValid = isSameTree(desc[0]->par);
		if(isValid == 0)
		{
			std::cout << "In node::validTree()\n"
					"The given tree ";
			print_coaltree();
			std::cout << " is invalid, because desc[0]->par is ";
			desc[0]->par->print_coaltree();
			std::cout <<"\n";
		}
		else if(desc[0]->validTree()== 0)
		{
			std::cout << "In node::validTree()\n"
					"The given tree ";
			desc[0]->print_coaltree();
			std::cout << " is invalid\n";
		}
		else
		{
			isValid = isSameTree(desc[1]->par);
			if(isValid == 0)
			{
				std::cout << "In node::validTree()\n"
							"The given tree ";
				print_coaltree();
				std::cout << " is invalid, because desc[1]->par is ";
				desc[1]->par->print_coaltree();
				std::cout <<"\n";
			}
			else if(desc[1]->validTree()== 0)
			{
				std::cout << "In node::validTree()\n"
						"The given tree ";
				desc[1]->print_coaltree();
				std::cout << " is invalid\n";
			}
		}
	}
	return isValid;
}

void node::print_nodeMembers()
{
	std::cout << "\nPrint all the node members\n";
	print_coaltree();
	std::cout << "isRoot: "<<isRoot<<"\n";
	std::cout << "isTip: "<<isTip<<"\n";
	std::cout << "tipID: "<<tipID<<"\n";
	std::cout << "age: "<<age<<"\n";
	std::cout << "isLikelihoodNULL: "<< isLikelihoodNULL <<"\n";
	if(!isTip & !isLikelihoodNULL)
	{
		std::cout << "Likelihoods:\n";
		std::cout << lik <<"\n";
		/*
		for(int i=0;i<lik.get_nrow();i++)
		{
			for(int j=0;j<lik.get_ncol();j++)
				std::cout << lik.val(i,j) <<" ";
			std::cout <<"\n";
		}
		*/
	}
	std::cout <<"\n";
}

void node_old::print_nodeMembers()
{
	std::cout << "\nPrint all the node members\n";
	print_coaltree();
	std::cout << "isRoot: "<<isRoot<<"\n";
	std::cout << "isTip: "<<isTip<<"\n";
	std::cout << "tipID: "<<tipID<<"\n";
	std::cout << "age: "<<age<<"\n";
	std::cout << "isLikelihoodNULL: "<< isLikelihoodNULL <<"\n";
	if(!isTip & !isLikelihoodNULL)
	{
		std::cout << "Likelihoods:\n";
		for(int i=0;i<lik.get_nrow();i++)
		{
			for(int j=0;j<lik.get_ncol();j++)
				std::cout << lik.val(i,j) <<" ";
			std::cout <<"\n";
		}
	}
	std::cout <<"\n";
}

void node::maxTipID(int &maxID)
{
	if(isTip==0)
	{
		desc[0]->maxTipID(maxID);
		desc[1]->maxTipID(maxID);
	}
	else if(isTip == 1)
	{
		if(tipID > maxID)
			maxID = tipID;
	}
	return;
}

/**
 *  Return 0 if the argument 'tree' is same as the current tree; 0 otherwise.
 *	@param node
 */
int node::isSameTree(node* tree)
{
	// REMOVE
	//std::cout << "In node::isSameTree().\n";
	//std::cout << "Compare the following two trees:\n";
	//print_coaltree();
	//std::cout << "and ";
	//tree->print_coaltree();

	int sameTree = 1;

	if((isTip ==1 && tree->isTip ==1) && tipID != tree->tipID)
	{
		sameTree = 0;
		return sameTree;
	}
	else if(size_tree() != tree->size_tree() || isTip != tree->isTip || isRoot != tree->isRoot
			|| age != tree->age)// || (lik.get_valarray() != (tree->lik).get_valarray()).sum() > 0 )
	{
		// class members are different
		sameTree = 0;
		return sameTree;
	}
	else
	{
		// FIXME YC 4/16/2014
		// So far, I haven't seen any two different trees belong to this case.
		//std::cout << "Warning in node::isSameTree(). Double check whether these trees are the same or not:\n";
		//print_coaltree();
		//std::cout << "and ";
		//tree->print_coaltree();

		// Warning!
		// Actually we need to compare the class members of descendant nodes of two trees we want to compare.
		// However, the likelihood values and ages of two trees are highly likely different
		// if two trees are different.
		return sameTree;
	}
}

int node::find_sisterID()
{
	// REMOVE
	//std::cout << "In node::find_sisterID()\n";
	//print_coaltree();

	if(!isSameTree(par->desc[0]))
		return 0;
	else
		return 1;
}



node* node::findNode(int tipid, int &found)
{
	if(tipID == tipid)
	{
		found = 1;
		return get_node();
	}
	else if(isTip == 0)
	{
		node* nd = desc[0]->findNode(tipid, found);
		if(found == 1)
			return nd;
		else // if(found == 0)
		{
			nd = desc[1]->findNode(tipid, found);
			return nd;
		}
	}
}

node* node::findNode_fromRoot(int tipid)
{
	if(tipID == tipid)
		return par->desc[siblingOrder];
	else
	{
		int found = 0;
		node* nd = desc[0]->findNode(tipid, found);
		if(found == 1)
			return nd;
		else // if(found == 0)
		{
			nd = desc[1]->findNode(tipid, found);
			if(found == 0)
				std::cout << "Error in node::findNode_fromRoot().\n";
			return nd;
		}

	}
}



node* node::deepCopy_root()
{
	node* copy_node = new node;
	copy_node->initialization();
	//node tree;
	//node* copy_node = &tree;
	//node* copy_node;

	copy_node->isTip = isTip;
	copy_node->isRoot = isRoot;
	copy_node->tipID = tipID;
	copy_node->age = age;
	copy_node->isLikelihoodNULL = isLikelihoodNULL;
	copy_node->totalCoalRate = totalCoalRate;
	copy_node->label = label;
	copy_node->popID = popID;
	copy_node->assign_siblingOrder(siblingOrder);
	if(!isLikelihoodNULL)
	{
		copy_node->lik.resize(4,lik.cols());
		copy_node->lik.setZero();
		copy_node->lik = lik;
		// copy_node->lik.replace(lik.get_valarray(),4,lik.get_ncol());
		// copy_node->lik=lik;
	}
	if(!isTip)
	{
		copy_node->desc[0] = desc[0]->deepCopy_root();
		copy_node->desc[1] = desc[1]->deepCopy_root();
		copy_node->desc[0]->par = copy_node;
		copy_node->desc[1]->par = copy_node;
	}
	return copy_node;
}


node_old* node_old::deepCopy_root()
{
	node_old* copy_node = new node_old;
	copy_node->initialization();
	//node tree;
	//node* copy_node = &tree;
	//node* copy_node;

	copy_node->isTip = isTip;
	copy_node->isRoot = isRoot;
	copy_node->tipID = tipID;
	copy_node->age = age;
	copy_node->isLikelihoodNULL = isLikelihoodNULL;
	copy_node->totalCoalRate = totalCoalRate;
	copy_node->label = label;
	copy_node->popID = popID;
	copy_node->assign_siblingOrder(siblingOrder);
	if(!isLikelihoodNULL)
	{
		copy_node->lik.replace(lik.get_valarray(),4,lik.get_ncol());
		// copy_node->lik=lik;
	}
	if(!isTip)
	{
		copy_node->desc[0] = desc[0]->deepCopy_root();
		copy_node->desc[1] = desc[1]->deepCopy_root();
		copy_node->desc[0]->par = copy_node;
		copy_node->desc[1]->par = copy_node;
	}
	return copy_node;
}

/**
 *  Deep copy of a tree (the current node and its descendants)
 */
node* node::deepCopy()
{
	node* copy_node = new node;
	copy_node->initialization();
	copy_node->assign_siblingOrder(siblingOrder);
	//node tree;
	//node* copy_node = &tree;
	//node* copy_node;
	if(isRoot == 1)
		copy_node = deepCopy_root();
	else if(isRoot == 0 )
	{
		int crrTipID = tipID;
		tipID = -1;

		node* root = go2root();
		copy_node = root->deepCopy_root();
		copy_node = copy_node->findNode_fromRoot(tipID);

		tipID = crrTipID;
		copy_node->tipID = crrTipID;
	}
	return copy_node;
}


/**
 *  Deep copy of a tree (the current node and its descendants)
 */
void node::deepCopy(node* tr)
{
	// REMOVE
	//std::cout << "In node::deepCopy()\n";
	//tr->print_coaltree();
	//tr->desc[0]->print_coaltree();
	//tr->desc[1]->print_coaltree();

	//initialization();

  assign_siblingOrder(tr->get_siblingOrder());
	//node tree;
	//node* copy_node = &tree;
	//node* copy_node;
  if(tr->isRoot == 1)
    deepCopy_root(tr);
  else if(tr->isRoot == 0 )
    {
		//int crrTipID = tr->tipID;
		//tr->tipID = -1;

      par->deepCopy(tr->par);
		//node* root;
		//root->deepCopy_root(tr->go2root());
		//copy_node = copy_node->findNode_fromRoot(tr->tipID);

		//tr->tipID = crrTipID;
		//tipID = crrTipID;
    }
}



void node::deepCopy_obj(node* tr)
{
	// REMOVE
	//std::cout << "In node::deepCopy_obj()\n";
	//tr->print_coaltree();
	//node* copy_node = new node;
	//initialization();
	//node tree;
	//node* copy_node = &tree;
	//node* copy_node;

	isTip = tr->isTip;
	isRoot = tr->isRoot;
	tipID = tr->tipID;
	age = tr->age;
	isLikelihoodNULL = tr->isLikelihoodNULL;
	totalCoalRate = tr->totalCoalRate;
	label = tr->label;
	popID = tr->popID;
	assign_siblingOrder(tr->get_siblingOrder());

	if(!tr->isLikelihoodNULL)
	{
		lik.resize(4,tr->lik.cols());
		lik = tr->lik;
		//lik.replace(tr->lik.get_valarray(),4,tr->lik.get_ncol());
		// copy_node->lik=lik;
	}
	if(tr->isTip == 0)
	{

		desc[0] = new node;
		desc[1] = new node;
		desc[0]->deepCopy_root(tr->desc[0]);
		desc[1]->deepCopy_root(tr->desc[1]);

		if(tr->isRoot == 0)
		{
			desc[1]->par = this;
			desc[0]->par = this;
		}
	}
}

void node_old::deepCopy_obj(node_old copiedTree, node_old* tr)
{
	// REMOVE
	//std::cout << "In node::deepCopy_obj()\n";
	//tr->print_coaltree();
	//node* copy_node = new node;
	//initialization();
	//node tree;
	//node* copy_node = &tree;
	//node* copy_node;

	isTip = tr->isTip;
	isRoot = tr->isRoot;
	tipID = tr->tipID;
	age = tr->age;
	isLikelihoodNULL = tr->isLikelihoodNULL;
	totalCoalRate = tr->totalCoalRate;
	label = tr->label;
	popID = tr->popID;
	assign_siblingOrder(tr->get_siblingOrder());

	if(!tr->isLikelihoodNULL)
	{
		lik.replace(tr->lik.get_valarray(),4,tr->lik.get_ncol());
		// copy_node->lik=lik;
	}
	if(tr->isTip == 0)
	{

		desc[0] = new node_old;
		desc[1] = new node_old;
		desc[0]->deepCopy_root(desc[0],tr->desc[0]);
		desc[1]->deepCopy_root(desc[1],tr->desc[1]);

		if(tr->isRoot == 0)
		{
			desc[1]->par = &copiedTree;
			desc[0]->par = &copiedTree;
		}
	}
}



void node::deepCopy_root(node* tr)
{
  isTip = tr->isTip;
  isRoot = tr->isRoot;
  tipID = tr->tipID;
  age = tr->age;
  isLikelihoodNULL = tr->isLikelihoodNULL;
  totalCoalRate = tr->totalCoalRate;
  label = tr->label;
  popID = tr->popID;
  assign_siblingOrder(tr->get_siblingOrder());

  numMutations_overLocus = tr->get_numMutations_overLocus();
  totalLengthsLineages = tr->get_totalLengthsLineages();
  std::vector<double> wt = tr->get_waitingTime();
  waitingTime.resize(wt.size());
  for(unsigned int i; i< wt.size(); i++)
    waitingTime.at(i) = wt.at(i);

  if(!tr->isLikelihoodNULL)
    {
      lik.resize(4,tr->lik.cols());
      lik.setZero();
      lik = tr->lik;
      //lik.replace(tr->lik.get_valarray(),4,tr->lik.get_ncol());
      // copy_node->lik=lik;
    }
  if(tr->isTip == 0)
    {
      delete desc[0];
      delete desc[1];
      desc[0] = new node;
      desc[0]->deepCopy_root(tr->desc[0]);
      desc[1] = new node;
      desc[1]->deepCopy_root(tr->desc[1]);
      desc[0]->par = this;
      desc[1]->par = this;
    }
}


void node_old::deepCopy_root(node_old* copiedTree, node_old* tr)
{
	// REMOVE
	//std::cout << "In node::deepCopy_root()\n";
	//tr->print_coaltree();
	//node* copy_node = new node;
	//initialization();
	//node tree;
	//node* copy_node = &tree;
	//node* copy_node;

	isTip = tr->isTip;
	isRoot = tr->isRoot;
	tipID = tr->tipID;
	age = tr->age;
	isLikelihoodNULL = tr->isLikelihoodNULL;
	totalCoalRate = tr->totalCoalRate;
	label = tr->label;
	popID = tr->popID;
	assign_siblingOrder(tr->get_siblingOrder());

	if(!tr->isLikelihoodNULL)
	{
		lik.replace(tr->lik.get_valarray(),4,tr->lik.get_ncol());
		// copy_node->lik=lik;
	}
	if(tr->isTip == 0)
	{
		delete desc[0];
		delete desc[1];
		desc[0] = new node_old;
		desc[0]->deepCopy_root(desc[0],tr->desc[0]);
		desc[1] = new node_old;
		desc[1]->deepCopy_root(desc[1],tr->desc[1]);
		desc[0]->par = copiedTree;
		desc[1]->par = copiedTree;
	}
}



/**
 *  shallow copy of a tree (the current node and its descendants)
 */
node* node::get_node()
{
	// FIXME
	node* copy_node = new node;
	copy_node->initialization();
	//node tmp;
	//node* copy_node = &tmp;

	copy_node->isTip = isTip;
	copy_node->isRoot = isRoot;
	if(!isRoot)
		copy_node->par = par;

	if(!isTip)
	{
		copy_node->desc[0] = desc[0];
		copy_node->desc[1] = desc[1];
		copy_node->isLikelihoodNULL = isLikelihoodNULL;
		if(!isLikelihoodNULL)
			copy_node->lik = lik;
	}
	else
		copy_node->tipID = tipID;
	copy_node->age = age;

	return copy_node;
}

Eigen::MatrixXd node::get_likMatrix()
{
	return lik;
}

Matrix node_old::get_likMatrix()
{
	return lik;
}

/**
 * Count the number of tips below the node
 */
unsigned int node::size_tree()
{
	unsigned int ntips = 0;
	if(isTip == 1)
		ntips++;
	else if(isTip == 0)
		ntips += desc[0]->size_tree()+desc[1]->size_tree();

	return ntips;
}


/**
 * Called by print_coaltree() for printing a tree in the Newick format.
 */
void node::print_node()
{
	if(isTip == 1)
	{
		std::cout << tipID << ":" << par->age;
	}
	else if(isTip == 0)
	{
		std::cout << "(";
		desc[0]->print_node();
		std::cout << ",";
		desc[1]->print_node();
		std::cout << "):" << par->age - age ;
	}
}

/**
 * Called by print_coaltree() for printing a tree in the Newick format.
 */
void node_old::print_node()
{
	if(isTip == 1)
	{
		std::cout << tipID << ":" << par->age;
	}
	else if(isTip == 0)
	{
		std::cout << "(";
		desc[0]->print_node();
		std::cout << ",";
		desc[1]->print_node();
		std::cout << "):" << par->age - age ;
	}
}


/**
 * Print out a coalescent tree in the Newick format
 */
void nodeSimple::print_topo()
{
	if(isTip == 1)
		std::cout << popID ;
	else if(isTip == 0)
	{
		std::cout << "(";
		firstChild->print_topo();
		std::cout << ",";
		secondChild->print_topo();
		std::cout << "):";
		std::cout << rank;
	}
	if(isRoot ==1)
		std::cout << ";\n";
}

/**
 * Prints out coalescent tree in the newick format to a file pointed by fp
 */
void nodeSimple::fileprint_topo(std::ofstream& fp)
{
	if(isTip == 1)
		fp << popID ;
	else if(isTip == 0)
	{
		fp << "(";
		firstChild->fileprint_topo(fp);
		fp << ",";
		secondChild->fileprint_topo(fp);
		fp << "):";
		fp << rank;
	}
	if(isRoot ==1)
		fp << ";\n";
}

/**
 * Print out a coalescent tree in the Newick format
 */
void node::print_coaltree()
{
	if(isTip == 1)
		std::cout << "(" << tipID << ");\n";
	else if(isTip == 0)
	{
		std::cout << "(";
		desc[0]->print_node();
		std::cout << ",";
		desc[1]->print_node();
		std::cout << ");\n";
	}
}
/**
 * Print out a coalescent tree in the Newick format
 */
void node_old::print_coaltree()
{
	if(isTip == 1)
		std::cout << "(" << tipID << ");\n";
	else if(isTip == 0)
	{
		std::cout << "(";
		desc[0]->print_node();
		std::cout << ",";
		desc[1]->print_node();
		std::cout << ");\n";
	}
}


/*
void node::initialize_lik(int n_uniqSites)
{
	lik.initialize(-1.0, 4 , n_uniqSites);
}
*/

unsigned int locus::get_uniqSeq(int idx_site, int idx_gene)
{
	return seq_uniq.at(idx_site).at(idx_gene);
}


valarray<double> locus::get_pi()
{
	return pi;
}

vector<unsigned int> locus::get_freq_uniqueSeq()
{
	return freq_uniqueSeq;
}



/**
 * Extract the smallest coalescent time from a tree in descending order
 */
double node::get_minCoalTime(double coalT)
{
	if(isTip)
		return coalT;
	else
	{
		if(age <= coalT || coalT < 0.0)
			coalT = age;

		coalT = desc[0]->get_minCoalTime(coalT);
		coalT = desc[1]->get_minCoalTime(coalT);
		return coalT;
	}
}




/**
 * Extract the coalescent times from the root of the tree to the tips and from left to right
 */
list<double> node::get_coalescentTimes(list<double> coalTimes)
{
  //list<double> coalTimes;
  coalTimes.push_back(age);
  if(isTip == 0)
    {
      unsigned int s1 = desc[0]->size_tree();
      unsigned int s2 = desc[1]->size_tree();
      if(s1 >1 && s2 <=1)
	coalTimes = desc[0]->get_coalescentTimes(coalTimes);
      else if(s1 <=1 && s2 >1)
	coalTimes=desc[1]->get_coalescentTimes(coalTimes);
      else if(s1>1 && s2 >1)
	{
	  if(s1 < s2)
	    {
	      coalTimes  = desc[0]->get_coalescentTimes(coalTimes);
	      coalTimes=desc[1]->get_coalescentTimes(coalTimes);
	    }
	  else if(s1 > s2)
	    {
	      coalTimes=desc[1]->get_coalescentTimes(coalTimes);
	      coalTimes = desc[0]->get_coalescentTimes(coalTimes);
	    }
	  else // if s1 == s2
	    {
	      if(desc[0]->age >= desc[1]->age) // older first
		{
		  coalTimes=desc[0]->get_coalescentTimes(coalTimes);
		  coalTimes =desc[1]->get_coalescentTimes(coalTimes);
		}
	      else
		{
		  coalTimes =desc[1]->get_coalescentTimes(coalTimes);
		  coalTimes=desc[0]->get_coalescentTimes(coalTimes);
		}
	    }
	}
    }
  return coalTimes;
}

void nodeSimple::deleteTopo()
{
	if(isTip == 0 )
	{
		firstChild->deleteTopo();
		secondChild->deleteTopo();
	}
	delete this;
}

std::vector<unsigned int> nodeSimple::get_nodePopID_withRank(unsigned int R)
{  
  std::vector<unsigned int> res;
  unsigned int id = 0;
  unsigned int found = 0;
  res.push_back(id); res.push_back(found);
  if(rank == R)
    {
      res.at(0) = nodePopID;
      res.at(1) = 1;
    }
  else
    {    
      if(isTip == 0)
	{
	  res = firstChild->get_nodePopID_withRank(R);
	  if(res.at(1) == 0)
	    res = secondChild->get_nodePopID_withRank(R);
	}
    }
  
  return res;
}




void nodeSimple::compute_areTipsFromSamePop()
{
  if(isTip==0)
    {
      firstChild->compute_areTipsFromSamePop();
      secondChild->compute_areTipsFromSamePop();
      unsigned int popID_1stChild = firstChild->get_nodePopID();
      unsigned int popID_2ndChild = secondChild->get_nodePopID();
      if(popID_1stChild !=0 && popID_1stChild ==popID_2ndChild)
	{
	  AreTipsFromSamePop = 1;
	  nodePopID = popID_1stChild;
	}
      else
	{
	  AreTipsFromSamePop = 0;
	  nodePopID = 0;
	}
    }
  else
    nodePopID = popID;
  return;
}

void nodeSimple::compute_nSamples_nSymNode(Eigen::MatrixXi &nSamples,unsigned int &nSymNode)
{
	if(isTip == 1)
	{
		nSamples(0,popID-1)++;
	}
	else
	{
		if(firstChild->get_isTip() == 1 && secondChild->get_isTip() == 1)
			if(firstChild->getPopID() == secondChild->getPopID())
				nSymNode++;
		firstChild->compute_nSamples_nSymNode(nSamples,nSymNode);
		secondChild->compute_nSamples_nSymNode(nSamples,nSymNode);
	}
	return;
}


unsigned int nodeSimple::compute_nTrees_sameTopo(unsigned int nPops)
{
  unsigned int ntrees = 1;
  Eigen::MatrixXi nSamples(1,nPops);
  nSamples.setZero();
  unsigned int nSymNode = 0;
  
  compute_nSamples_nSymNode(nSamples,nSymNode);
  for(unsigned int i = 0; i<nPops; i++)
    {
      if(nSamples(0,i)!=0)
	ntrees *= factorial((unsigned) nSamples(0,i));
    }
  ntrees /= pow(2,nSymNode);
  
  return ntrees;
}

