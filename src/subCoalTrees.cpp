/* MIST copyright 2017 by Yujin Chung and Jody Hey */

#include <iostream>
#include "Chain.hpp"
#include "misc.hpp"

/* This file contains functions to have subtrees of coalescent trees */

// subtree of "subSize" of "tree"
void Chain::getSubtree( unsigned int subSize, node* tree, unsigned int sampleID, unsigned int subLocusID)
{
  unsigned int size = tree->size_tree(); // the number of tips
  #ifdef DEBUG
  //std::cout <<"In Chain::getSubtree()\n";
  /*
  std::cout <<"subSize = " << subSize <<" size of the give tree = "<< size <<".\n";
  tree->print_coaltree();
  */
#endif //DEBUG
   
  std::vector<node*> subtrees;  
  if(subSize >= size)
    {
      std::cout <<"\n\nError in Chain::getSubtree. 'subSize' should be smaller than the size of 'tree'.\n";
      std::cout <<"'subSize' = " << subSize <<" and the size of 'tree' = "<< size <<"\n\n";
    }
  else // subsize < size
    {
      unsigned int leftsize = min(subSize,tree->desc[0]->size_tree());
      unsigned int rightsize = min(subSize-leftsize,tree->desc[1]->size_tree());
      #ifdef DEBUG
      //    std::cout <<"leftsize = "<< leftsize <<" rightsize = "<< rightsize <<"\n";
#endif //DEBUG
      for(int ll=leftsize; ll>=0; ll--)
	{	  
	  node* subtr = new node;
	  int rr=subSize-ll;
      #ifdef DEBUG
	  //	  	   std::cout <<"ll = "<< ll <<" rr = "<< rr <<"\n";
#endif //DEBUG
	  if(rr<0 || rr> tree->desc[1]->size_tree())
	    ll=-1;
	  else
	    {
	      subtr = tree->deepCopy_root();
	      if(ll==0)
		{
		  if(rr == subtr->desc[1]->size_tree())
		    {
		      subtr->desc[0]->deleteCoalTree();
		      subtr = subtr->desc[1];
		      delete subtr->par;
		      subtr->set_isRoot(1);
		      subtr->assignPopulations2Tips(SeqPop);
		      subtrees.push_back(subtr);		  
		    }
		  else
		    {	      
		      std::vector<node*> subsubtrees;
		      subsubtrees = subtr->desc[1]->getSubtree(rr);
		      for(unsigned int i=0; i<subsubtrees.size(); i++)
			{
			  subsubtrees.at(i)->assignPopulations2Tips(SeqPop);
			  subtrees.push_back(subsubtrees.at(i));
			  #ifdef DEBUG
			  /*
			  subsubtrees.at(i)->print_coaltree();
			  std::cout <<"subtrees.size() = "<<subtrees.size() << "\n";
			  for(unsigned int i=0; i<subtrees.size(); i++)
			    {
			      std::cout <<"i="<<i<<" ";
			      subtrees.at(i)->print_coaltree();
			    }
			  */
#endif //DEBUG
			}
		      subtr->deleteCoalTree();
		    }
		}
	      else if(rr==0)
		{
		  if(ll == subtr->desc[0]->size_tree())
		    {
		      subtr->desc[1]->deleteCoalTree();
		      subtr = subtr->desc[0];		  
		      delete subtr->par;
		      subtr->set_isRoot(1);
		      subtr->assignPopulations2Tips(SeqPop);
		      subtrees.push_back(subtr);
		    }
		  else
		    {		      
		      std::vector<node*> subsubtrees;
		      subsubtrees = tree->desc[0]->getSubtree(ll);
		      for(unsigned int i=0; i<subsubtrees.size(); i++)
			{
			  subsubtrees.at(i)->assignPopulations2Tips(SeqPop);
			  subtrees.push_back(subsubtrees.at(i));
			}
		      subtr->deleteCoalTree();
		    }
		}
	      else if(ll>0 & rr>0)
		{
		  #ifdef DEBUG
		  // std::cout <<"subtr->desc[0] is ";
		  // subtr->desc[0]->print_coaltree();
#endif //DEBUG
		  std::vector<node*> Lsubtrees;
		  if(ll ==  subtr->desc[0]->size_tree())
		    {
		      node* Ltree = new node;
		      Ltree = subtr->desc[0]->deepCopy_root();
		      Lsubtrees.push_back(Ltree);		      
		    }
		  else if(ll < subtr->desc[0]->size_tree())
		    {
		      Lsubtrees = subtr->desc[0]->getSubtree(ll);
		    }
		  
		  #ifdef DEBUG
		  // std::cout <<"subtr->desc[1] is ";
		  // subtr->desc[1]->print_coaltree();
#endif //DEBUG
		  std::vector<node*> Rsubtrees;
		  if(rr ==  subtr->desc[1]->size_tree())
		    {
		      node* Rtree = new node;
		      Rtree = subtr->desc[1]->deepCopy_root();
		      Rsubtrees.push_back(Rtree);		      
		    }
		  else if(rr < subtr->desc[1]->size_tree())
		    {		    		      
		      Rsubtrees = subtr->desc[1]->getSubtree(rr);
		    }
#ifdef DEBUG
		  //   std::cout <<"Lsubtrees.size() = "<<Lsubtrees.size()
		  //	    <<" Rsubtrees.size() = "<< Rsubtrees.size() <<"\n";
#endif //DEBUG
		  for(unsigned int i=0; i<Lsubtrees.size(); i++)
		    {
		      for(unsigned int j=0; j<Rsubtrees.size(); j++)
			{
			  node* tempTr = new node;
			  tempTr = subtr->deepCopy_root();
			  tempTr->desc[0]->deleteCoalTree();
			  tempTr->desc[1]->deleteCoalTree();
			  tempTr->desc[0] = Lsubtrees.at(i);
			  tempTr->desc[1] = Rsubtrees.at(j);
			  tempTr->set_isRoot(1);
			  tempTr->desc[0]->par = tempTr;
			  tempTr->desc[1]->par = tempTr;
			  tempTr->assignPopulations2Tips(SeqPop);
			  subtrees.push_back(tempTr);

#ifdef DEBUG
			  /*
			  std::cout <<"subtrees.size() = "<<subtrees.size() << "\n";
			  Rsubtrees.at(j)->print_coaltree();
			      subtrees.at(subtrees.size()-1)->print_coaltree();
			  */
			  /*
			  for(unsigned int k=0; k<subtrees.size(); k++)
			    {
			      std::cout <<"k="<<k<<" ";
			      subtrees.at(k)->print_coaltree();
			    }
			  */
#endif //DEBUG		  
			}
		    }
		  subtr->deleteCoalTree();		  
		}   		  
	    }
	}
    }
  
  
   #ifdef DEBUG
  /*
  std::cout << "\nThere are " <<subtrees.size()<< " subtrees:\n";
  for(unsigned int i=0; i<subtrees.size(); i++)
    {
      std::cout <<"i="<<i<<" ";
      subtrees.at(i)->print_coaltree();
    }
  std::cout <<"\n\n";
  */
   #endif //DEBUG

  Coalsubtrees.at(sampleID).at(subLocusID) = subtrees;
  sizeForest = subtrees.size();

  #ifdef DEBUG
  // std::cout <<"sizeForest = "<< sizeForest <<"\n";
#endif //DEBUG

  return;
}


std::vector<node*> node::getSubtree(unsigned int subSize)
{  
  unsigned int size = size_tree(); // the number of tips
  
#ifdef DEBUG
  //  std::cout << "\tin node::getSubtree()\n";
  //  std::cout <<"\tsubSize = " << subSize <<" size of the give tree = "<< size <<".\n\t";
  //  print_coaltree();
#endif //DEBUG
  
  std::vector<node*> subtrees;
 
  if(subSize > size_tree())
    {
      std::cout <<"\n\nError in node::getSubtree. 'subSize' should be smaller than the size of 'tree'.\n";
      std::cout <<"'subSize' = " << subSize <<" and the size of 'tree' = "<< size_tree() <<"\n\n";
    }
  else
    {
      unsigned int leftsize = min(subSize,desc[0]->size_tree());
      unsigned int rightsize = min(subSize-leftsize,desc[1]->size_tree());
#ifdef DEBUG
      //      std::cout <<"\tleftsize = "<< leftsize <<" rightsize = "<< rightsize <<"\n";
#endif //DEBUG
      for(int ll=leftsize; ll>=0; ll--)
	{	  
	  int rr=subSize-ll;
#ifdef DEBUG
	  // 	    std::cout <<"\tll = "<< ll <<" rr = "<< rr <<"\n";
#endif //DEBUG
	  if(rr<0 || rr> desc[1]->size_tree())
	    ll=-1;
	  else
	    {
	      node* tempTr = new node;
	      tempTr = deepCopy_root();
	      if(ll==0) // need to delete 'this'
		{	      
		  if(rr == tempTr->desc[1]->size_tree())
		    {
		      node* subtr = new node;
		      subtr = tempTr->desc[1]->deepCopy_root();
		      subtrees.push_back(subtr);
		    }
		  else
		    {      
		      std::vector<node*> subsubtrees;
		      subsubtrees = tempTr->desc[1]->getSubtree(rr);
		      for(unsigned int i=0; i<subsubtrees.size(); i++)
			{
			  subtrees.push_back(subsubtrees.at(i));
			}

		    }
		}
	      else if(rr==0) // need to delete 'this'
		{
		  if(ll == tempTr->desc[0]->size_tree())
		    {
		      node* subtr = new node;
		      subtr = tempTr->desc[0]->deepCopy_root();
		      subtrees.push_back(subtr);      
		    }
		  else
		    { 
		      std::vector<node*> subsubtrees;
		      subsubtrees = tempTr->desc[0]->getSubtree(ll);
		      for(unsigned int i=0; i<subsubtrees.size(); i++)
			{
			  subtrees.push_back(subsubtrees.at(i));
			}
		    }
		}
	      else if(ll>0 & rr>0) // need to keep 'this'
		{
		  std::vector<node*> Lsubtrees;
		  if(ll == desc[0]->size_tree())
		    {
		      node* Ltree = new node;
		      Ltree = desc[0]->deepCopy_root();
		      Lsubtrees.push_back(Ltree);		      
		    }
		  else if(ll < desc[0]->size_tree())
		    {
		      Lsubtrees = tempTr->desc[0]->getSubtree(ll);
		    }
		  std::vector<node*> Rsubtrees;
		  if(rr == desc[1]->size_tree())
		    {
		      node* Rtree = new node;
		      Rtree = desc[1]->deepCopy_root();
		      Rsubtrees.push_back(Rtree);		      
		    }
		  else if(rr < desc[1]->size_tree())
		    {		    		      
		      Rsubtrees = tempTr->desc[1]->getSubtree(rr);
		    }
		  for(unsigned int i=0; i<Lsubtrees.size(); i++)
		    {
		      for(unsigned int j=0; j<Rsubtrees.size(); j++)
			{
			  node* subsubtr = new node;
			  subsubtr = tempTr->deepCopy_root();
			  subsubtr->desc[0]->deleteCoalTree();
			  subsubtr->desc[1]->deleteCoalTree();
			  subsubtr->desc[0] = Lsubtrees.at(i);
			  subsubtr->desc[1] = Rsubtrees.at(j);
			  subsubtr->desc[0]->par = subsubtr;
			  subsubtr->desc[1]->par = subsubtr;
			  subtrees.push_back(subsubtr);			  
			}
		    }	  
		}
	      tempTr->deleteCoalTree();
	    }
	}
    }
  #ifdef DEBUG
  /*
  std::cout << "\n\tThere are " <<subtrees.size()<< " subtrees:\n";
  for(unsigned int i=0; i<subtrees.size(); i++)
    {
      std::cout <<"\ti="<<i<<" ";
      subtrees.at(i)->print_coaltree();
    }
  std::cout <<"\n\n";
  */
#endif // DEBUG
  
  // this->deleteCoalTree(); // doublecheck the memory leak - fixit YC 3/10/2017
  return subtrees;
}





void Chain::collectAllUpdates_Lmode_subtrees( unsigned int savingID)
{
 
  subtreeIDs.at(savingID).resize(numSubLoci);
  subcoalTimes.at(savingID).resize(numSubLoci);
  subtipIDs.at(savingID).resize(numSubLoci);

  if(lociInParallel ==0) 
    logPrior_trees.at(savingID) = logPriorTrees_atPrev; // which is 0.
 
  for(unsigned int i=0; i < numSubLoci; i++)
    {
      subtreeIDs.at(savingID).at(i).resize(sizeForest);
      subcoalTimes.at(savingID).at(i).resize(sizeForest);
      subtipIDs.at(savingID).at(i).resize(sizeForest);
      
      for(unsigned int j=0; j< sizeForest; j++)
	{
	  node *subtr = Coalsubtrees.at(savingID).at(i).at(j);
	  compute_coalTimes_tipIDs_subtrees(savingID,i,j,subtr);
	  unsigned int list_size = list_trees.size();
#ifdef DEBUG
	  /*
	  std::cout <<"i(locus id) = "<< i <<" j(subtree id)="<<j;
	  std::cout <<" list_size = "<< list_size <<"\n";
	  std::cout <<" subtree is ";
	  subtr->print_coaltree();
	  */
#endif //DEBUG
	  if(list_size == 0)
	    {
	      nodeSimple *topo = new nodeSimple;
	      std::list<double> ct = subcoalTimes.at(savingID).at(i).at(j);
	      ct.sort();
	      topo->convert(subtr, ct, numTotalSeq); // Updated by YC 5/8/2014
	      topo->computeSizes();
	      list_trees.push_back(topo);
	      subtreeIDs.at(savingID).at(i).at(j) = 0;
	    }
	  else
	    {
	      unsigned int found_sampleTopo = 0;
	      unsigned int count = 0;
	      while(found_sampleTopo == 0 && count < list_size)
		{
		  found_sampleTopo = list_trees.at(count)->sameTopo(subtr);
		  count++;
		}

	      #ifdef DEBUG
	      // std::cout <<"found_sampleTopo = "<< found_sampleTopo <<"\n";
#endif //DEBUG

	      if(found_sampleTopo == 0)
		{
		  nodeSimple *topo = new nodeSimple;
		  std::list<double> ct = subcoalTimes.at(savingID).at(i).at(j);
		  ct.sort();
		  topo->convert(subtr,ct,numTotalSeq);// updated by YC 5/8/2014
		  topo->computeSizes();
		  list_trees.push_back(topo);
		  subtreeIDs.at(savingID).at(i).at(j) = list_trees.size()-1;
		}
	      else if(found_sampleTopo == 1)
		{
		  subtreeIDs.at(savingID).at(i).at(j) = count -1;
		}
	      else
		{
		  std::cout << "\n*** Error in Chain::collectAllUpdates ***\n";
		  std::cout << "found_sampleTopo = " << found_sampleTopo << "\n";
		  std::cout << "subtree = ";
		  subtr->print_coaltree();
		  nodeSimple *topo = new nodeSimple;
		  std::list<double> ct = subcoalTimes.at(savingID).at(i).at(j);
		  ct.sort();
		  topo->convert(subtr,ct, numTotalSeq);
		  topo->computeSizes();
		  std::cout << "topo = " ;
		  topo->print_topo();
		  std::cout <<  "list_size = " << list_size <<"\n";
		  for(unsigned int ll=0; ll<list_size; ll++)
		    {
		      std::cout << "list topo ll = " << ll << " ";
		      list_trees.at(ll)->print_topo();
		    }
		}
	    }
	  #ifdef DEBUG
	  // std::cout <<"The size of listSubtrees is " << listSubtrees.size() <<"\n";
	  // listSubtrees.at(listSubtrees.size()-1)->print_topo();
#endif //DEBUG
	}
    }

  #ifdef DEBUG
  /*
  for(unsigned int i=0; i< listSubtrees.size(); i++)
    {
      std::cout <<"i = "<< i <<" ";
      listSubtrees.at(i)->print_topo();
    }
  */
#endif // DEBUG
  
  return;
}

  
void Chain::compute_coalTimes_tipIDs_subtrees( unsigned int iter, unsigned int id_locus, unsigned int id_subtr, node* tree)
{
  if(tree->size_tree() > 1)
    {
      if(subcoalTimes.at(iter).size() == 0)
	subcoalTimes.at(iter).resize(n_loci);
      if(subcoalTimes.at(iter).at(id_locus).size() == 0)
	subcoalTimes.at(iter).at(id_locus).resize(Coalsubtrees.at(iter).at(id_locus).size());
      subcoalTimes.at(iter).at(id_locus).at(id_subtr).push_back(tree->age);
    }

  if(tree->isTip == 0)
    {
      if(tree->desc[0]->age < tree->desc[1]->age) // younger child is the first
	{
	  compute_coalTimes_tipIDs_subtrees(iter,id_locus,id_subtr,tree->desc[0]);
	  compute_coalTimes_tipIDs_subtrees(iter,id_locus,id_subtr,tree->desc[1]);
	}
      else if(tree->desc[0]->age > tree->desc[1]->age) 
	{
	  compute_coalTimes_tipIDs_subtrees(iter,id_locus,id_subtr,tree->desc[1]);
	  compute_coalTimes_tipIDs_subtrees(iter,id_locus,id_subtr,tree->desc[0]);
	}
      else if(tree->desc[0]->isTip == 1 && tree->desc[1]->isTip ==1) 
	{
	  if(tree->desc[0]->tipID < tree->desc[1]->tipID) // smaller tipID is the first
	    {
	      compute_coalTimes_tipIDs_subtrees(iter,id_locus,id_subtr,tree->desc[0]);
	      compute_coalTimes_tipIDs_subtrees(iter,id_locus,id_subtr,tree->desc[1]);	      
	    }
	  else if(tree->desc[0]->tipID > tree->desc[1]->tipID) 
	    {
	      compute_coalTimes_tipIDs_subtrees(iter,id_locus,id_subtr,tree->desc[1]);
	      compute_coalTimes_tipIDs_subtrees(iter,id_locus,id_subtr,tree->desc[0]);
	    }
	  else
	    {
	      std::cout << "\n*** Error in Chain::compute_coalTimes_tipIDs() ***\n";
	      std::cout << "(sub)tree is ";
	      tree->print_coaltree();
	      std::cout << "tree->desc[0]->age = " << tree->desc[0]->age
			<< "tree->desc[1]->age = " << tree->desc[1]->age
			<< "tree->desc[0]->isTip = "<< tree->desc[0]->isTip
			<< "tree->desc[1]->isTip = "<< tree->desc[1]->isTip
			<< "\n";
	    }
	}
      else if(tree->desc[0]->age == tree->desc[1]->age && tree->desc[0]->isTip == 0 && tree->desc[1]->isTip ==0)
	{
	  // same ages of two coalescent events
	  // It can happen when the mles of trees are analyzed in step 2. - YC 8/17/2016
	  tree->desc[1]->age += pow(10,-10);
	  if(tree->age < tree->desc[1]->age)
	    {
	      tree->desc[1]->age =(tree->age+ tree->desc[1]->age)/2;
	    }

	  // REMOVE
	  // std::cout <<" tree->desc[0]->age = " << tree->desc[0]->age 
	  //	    <<" tree->desc[1]->age =" << tree->desc[1]->age <<"\n";
	  compute_coalTimes_tipIDs_subtrees(iter,id_locus,id_subtr,tree->desc[0]);
	  compute_coalTimes_tipIDs_subtrees(iter,id_locus,id_subtr,tree->desc[1]);	 
	}
      else
	{
	  std::cout << "/n*** Error in Chain::compute_coalTimes_tipIDs() ***\n";
	  std::cout << "(sub)tree is ";
	  tree->print_coaltree();
	  std::cout << "tree->isTip = "<< tree->isTip
		    << "tree->desc[0]->age = " << tree->desc[0]->age
		    << "tree->desc[1]->age = " << tree->desc[1]->age
		    << "tree->desc[0]->isTip = "<< tree->desc[0]->isTip
		    << "tree->desc[1]->isTip = "<< tree->desc[1]->isTip
		    << "\n";  	  
	}
    }
  else if (tree->isTip == 1)
    {
      if(subtipIDs.at(iter).size() == 0)
	subtipIDs.at(iter).resize(n_loci);
      if(subtipIDs.at(iter).at(id_locus).size() == 0)
	subtipIDs.at(iter).at(id_locus).resize(Coalsubtrees.at(iter).at(id_locus).size());
      subtipIDs.at(iter).at(id_locus).at(id_subtr).push_back(tree->tipID);
    }
  else
    {
      std::cout << "/n*** Error in Chain::compute_coalTimes_tipIDs() ***\n";
      std::cout << "(sub)tree is ";
      tree->print_coaltree();
      std::cout << "tree->isTip = "<< tree->isTip
		<< "\n";      
    }

  return;
}

void Chain::collectAllUpdates_bwProcs_Lmode_subtrees()
{
  #ifdef DEBUG
  // std::cout <<"\n In Chain::collectAllUpdates_bwProcs_Lmode_subtrees()\n";
#endif //DEBUG
  for(unsigned int p=1; p< nProcesses; p++)
    {
      MPI::COMM_WORLD.Barrier();
      unsigned int list_size = list_trees.size();
      unsigned int listSize_new = 0;
      if(cpuID==0) // send list_trees
	{
	  MPI::COMM_WORLD.Send(&list_size, 1, MPI::UNSIGNED, p, 3168);
	  
	  for(unsigned int ll=0; ll<list_size; ll++)
	    list_trees.at(ll)->MPIsend_nodeSimple(p);	  
	}
      else if(cpuID==p) // receive list_trees
	{
	  list_size = 0;
	  MPI::COMM_WORLD.Recv(&list_size, 1, MPI::UNSIGNED, 0, 3168);
	  
	  for(unsigned int ll=0; ll<list_size; ll++)
	    {
	      nodeSimple *topo = new nodeSimple;
	      topo->MPIreceive_coaltree(0);
	      list_trees.push_back(topo);
	    }

	  collectAllUpdates_Lmode_subtrees(0);
	  listSize_new = list_trees.size();
	}
      
      
      MPI::COMM_WORLD.Bcast(&listSize_new, 1, MPI::UNSIGNED, p);
      MPI::COMM_WORLD.Barrier();

      if(cpuID <= p)
	{
	  if(listSize_new > list_size) // new unique topologies
	    {
	      if(cpuID==p)
		{
		  for(unsigned int qq=0; qq<p; qq++)
		    for(unsigned int tt=list_size; tt<listSize_new; tt++)
		      list_trees.at(tt)->MPIsend_nodeSimple(qq);		    
		}
	      for(unsigned int qq=0; qq<p; qq++)	      
		{
		  if(cpuID==qq)
		    {
		      for(unsigned int tt=list_size; tt<listSize_new; tt++)
			{
			  nodeSimple *topo = new nodeSimple;
			  topo->MPIreceive_coaltree(p);
			  list_trees.push_back(topo);
			}
		    }
		}
	    }	  
	}
    }
}






// YC 5/17/2017
// 'crrProcID' is required for debugging only. This function is called on each process.
double Chain::compute_logConditionalProb_subtrees(unsigned int id_sample, unsigned int id_locus, popTree* poptree, unsigned int crrProcID)
{

  #ifdef DEBUG
  // std::cout <<"Chain::compute_logConditionalProb_subtrees()\n";
#endif //DEBUG
  /*
  std::chrono::high_resolution_clock::time_point start_t, end_t;
  std::chrono::high_resolution_clock::time_point start_t_case, end_t_case;
  start_t= std::chrono::high_resolution_clock::now();
  */


  unsigned int nPops = poptree->size();
  // std::cout << "nPops = " << nPops <<"\n";
  
  // FIXME YC 5/9/2014
  // It works for up to 2 populations
  // Note: It might be better to have a migration rate matrix and a vector of population sizes
  // as members of class 'popTree'.
  Eigen::MatrixXd popSize(1,nPops);
  popSize.setZero();
  for(unsigned int p1=0; p1< nPops; p1++)
    {
      popSize(0,p1) = poptree->find_popSize(p1+1);
    }
  double splittingTime = poptree->get_age(); // Sorts the elements in ascending order

  #ifdef DEBUG
  // std::cout << "popSize = " << popSize <<"\n";
  // std::cout << "splittingTime = " << splittingTime <<"\n";
#endif //DEBUG

  double logPseudoProb = 0;

  for(unsigned int id_st=0; id_st<sizeForest; id_st++)
    {
      double logProb = 0;
      unsigned int trID = subtreeIDs.at(id_sample).at(id_locus).at(id_st);
      std::list<double> eventT = subcoalTimes.at(id_sample).at(id_locus).at(id_st);
      
      eventT.sort();
      
      #ifdef DEBUG
      //std::cout <<"id_st = "<< id_st << " trID = " << trID <<"\n";
#endif //DEBUG
      // std::cout <<"nGeneCopies = " << coalTimes.at(id_sample).at(id_locus).size()+1 <<"\n";
      
      std::list<double>::iterator iter = eventT.begin();
      
      unsigned int nGeneCopies = subcoalTimes.at(id_sample).at(id_locus).at(id_st).size()+1; 
      unsigned int ancestralPop = 0;
      Eigen::MatrixXcd probMat;
      unsigned int count_events=0;
      iter = eventT.begin();
      while(count_events <nGeneCopies-1 && iter!=eventT.end())
	{
	  #ifdef DEBUG
	  // std::cout <<"count_events = "<< count_events << " *iter = " << *iter << "\n";
#endif //DEBUG
	  if(*iter <=splittingTime)
	    {
	      #ifdef DEBUG
	      //	      std::cout <<"--- Case 3: coalescent events in sampling populations ---\n";	      
#endif //DEBUG
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
		  // count_case3_1++;	      
		  // start_t_case= std::chrono::high_resolution_clock::now();
		  
		  //--- Computing V*exp(D*t) ---//
		  V = V.cwiseProduct((eigenval*waitingTime).exp().matrix().transpose());
		  probMat = V*V_inv;
		  
		  // end_t_case= std::chrono::high_resolution_clock::now();	      
		  // eachComputingTime_condiProb_case3_1 += std::chrono::duration_cast<std::chrono::microseconds>(end_t_case - start_t_case);
		}
	      else // non-first coalescent event
		{
		  //std::cout << "Case3-2 : count_events = "<< count_events<< "\n";
		  // count_case3_2++;
		  // start_t_case= std::chrono::high_resolution_clock::now();
		  
		  //Eigen::MatrixXcd expD = (eigenval*waitingTime).exp().matrix().asDiagonal();
		  if(probMat.rows() != 1 || probMat.cols() != V.rows())
		    {
		  std::cout << "\n ***** Error in Chain::compute_conditionalProb_MPI *****\n";	
		  std::cout << "\n ***** Case 3-2: different dimensions *****\n";
		  std::cout << "\n ***** id_sample = "<< id_sample
			    << ", id_locus = "<< id_locus <<", crrProcID = "<< crrProcID << " and poptree is ";
		  poptree->print_allPopTree();
		  #ifdef DEBUG
		  // std::cout <<"V = "<< V<<"\n";
#endif //DEBUG
		    }
		  probMat = probMat*V; // row-vector
		  probMat = probMat.cwiseProduct((eigenval*waitingTime).exp().matrix().transpose());
		  probMat = probMat*V_inv;	      
		  
		  // end_t_case= std::chrono::high_resolution_clock::now();	      
		  // eachComputingTime_condiProb_case3_2 += std::chrono::duration_cast<std::chrono::microseconds>(end_t_case - start_t_case);
		}
	      
		  #ifdef DEBUG
	      //  std::cout <<" waitingTime= " << waitingTime << " and probMat=" << probMat <<"\n";
#endif //DEBUG
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

	      #ifdef DEBUG
	      /*
	      std::cout << "Case 2\n";
	  std::cout <<"count_events = "<< count_events << " *iter = " << *iter << "\n";
	      */
#endif //DEBUG
	      
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
		      #ifdef DEBUG
	      /*
	      std::cout << "V = "<< V <<"\n";
	      std::cout << "V_inv = "<< V_inv <<"\n";
		      std::cout << "eigenval = "<< eigenval <<"\n";
		      std::cout << "waitingTime = "<< waitingTime <<"\n";
	      */
#endif //DEBUG
	      if(count_events==0)
		{
		  // count_case2_1++;
		  // start_t_case= std::chrono::high_resolution_clock::now();
		  
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
		  
		  // end_t_case= std::chrono::high_resolution_clock::now();	      
		  // eachComputingTime_condiProb_case2_1 += std::chrono::duration_cast<std::chrono::microseconds>(end_t_case - start_t_case);
		}
	      else // if(count_events != 0 && waitingTime > 0.0)
		{
		  #ifdef DEBUG		  
		  //std::cout << "Case2-2 : count_events = "<< count_events<< "\n";
#endif //DEBUG
		  // count_case2_2++;
		  // start_t_case= std::chrono::high_resolution_clock::now();
		  
		  
		  if(probMat.rows()!=1 || probMat.cols() != V.rows())
		    {
		      std::cout << "\n ***** Error in Chain::compute_conditionalProb_MPI *****\n";	
		      std::cout << "\n ***** Case 2-2: different dimensions *****\n";
		      #ifdef DEBUG
		      /*
		      std::cout << "probMat = "<< probMat <<"\n";
		      std::cout << "V = "<< V <<"\n";
		      std::cout << "V_inv = "<< V_inv <<"\n";
		      std::cout << "eigenval = "<< eigenval <<"\n";
		      std::cout << "waitingTime = "<< waitingTime <<"\n";
		      */
#endif //DEBUG
		    }
		  probMat *= V; // row-vector
		  probMat = probMat.cwiseProduct((eigenval*waitingTime).exp().matrix().transpose());
		  probMat *= V_inv;	   
		  
		  // end_t_case= std::chrono::high_resolution_clock::now();	      
		  // eachComputingTime_condiProb_case2_2 += std::chrono::duration_cast<std::chrono::microseconds>(end_t_case - start_t_case);
		}
	      
	      double prob_noAbsorbing = probMat.sum().real();
	      /*
		for(unsigned int id_No_AbsorbingState = 0; id_No_AbsorbingState < probMat.cols()-1; id_No_AbsorbingState++)
		prob_noAbsorbing += probMat(0,id_No_AbsorbingState).real();
	      */
	      
	      logProb += log(prob_noAbsorbing);

	      #ifdef DEBUG
	      // std::cout <<"logProb = "<<logProb <<"\n";
#endif //DEBUG
	      
	  // The flag that 'Case 1' should be called.
	      ancestralPop =1;
	    }
	  else if(splittingTime==0 || *iter > splittingTime && ancestralPop == 1)
	    {
	      #ifdef DEBUG
	      // std::cout <<"--- Case 1: coalescents in the ancestral population. ---\n";
#endif //DEBUG
	      //--- Case 1: coalescents in the ancestral population. ---//
	      //  Some or all coalescent events happened in the ancestral population.
	      	      
	      // count_case1++;
	      // start_t_case= std::chrono::high_resolution_clock::now();
	      
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
	      // unsigned int trID = subtreeIDs.at(id_sample).at(id_locus).at(id_st);
	      double log_numSameCoalEvents =0;
	      for(unsigned int c=count_events; c < nGeneCopies-2; c++)
		{
		  log_numSameCoalEvents += log(static_cast<double> (numSameCoalEvents_inAncPop.at(trID).at(c)) );
		}
	      logProb += log_numSameCoalEvents;
	      
	      count_events = nGeneCopies-1;
	      ancestralPop = 1;

	      /*
	      end_t_case= std::chrono::high_resolution_clock::now();	      
	      eachComputingTime_condiProb_case1 += std::chrono::duration_cast<std::chrono::microseconds>(end_t_case - start_t_case);
	      */
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
      logPseudoProb += logProb;
    }
  /*
    end_t= std::chrono::high_resolution_clock::now();
    eachComputingTime_condiProb += std::chrono::duration_cast<std::chrono::microseconds>(end_t - start_t);
    count_condiProbFunctionCalls++;
  */
  
  // std::cout << "logProb = " << logProb <<"\n";
  
  #ifdef DEBUG
  //  std::cout <<"Exiting Chain::compute_logConditionalProb_subtrees()\n";
#endif //DEBUG
  return logPseudoProb;
}


long double Chain::compute_logConditionalProb_zeroMig_subtrees(unsigned int id_sample, unsigned int id_locus, popTree* poptree, unsigned int crrProcID)
{
  #ifdef DEBUG
  // std::cout << "Chain::compute_logConditionalProb_zeroMig_subtrees()\n";
  // poptree->print_poptree();
#endif //DEBUG
  
  unsigned int nPops = poptree->size();
  unsigned int nTotalPops = nPops*(nPops+1)/2;
  std::vector<double> popSize;
  for(unsigned int p1=0; p1< nTotalPops; p1++)
    {
      popSize.push_back(poptree->find_popSize(p1+1));
    }
  double splittingTime = poptree->get_age(); // Sorts the elements in ascending order

      #ifdef DEBUG
  // std::cout <<"sizeForest = " <<sizeForest <<"\n";
#endif //DEBUG
  long double logPseudoProb =0;
  for(unsigned int id_st=0; id_st<sizeForest; id_st++)
    {
      #ifdef DEBUG
      // std::cout <<"id_st = " <<id_st <<"\n";
#endif //DEBUG
      std::list<double> eventT = subcoalTimes.at(id_sample).at(id_locus).at(id_st);
      eventT.sort();
      std::vector<double> eventT_inc;
      std::list<double>::iterator iter = eventT.end();
      for(iter=eventT.begin(); iter!=eventT.end(); iter++)
	{
	  eventT_inc.push_back(*iter);
      #ifdef DEBUG
	  //  std::cout << "coaltime = "<< *iter <<"\n";
	  #endif
	}
      
      unsigned int trID = subtreeIDs.at(id_sample).at(id_locus).at(id_st);
      #ifdef DEBUG
      // std::cout <<"trID = "<<trID <<"\n";
#endif //DEBUG
      nodeSimple* topo=list_trees.at(trID);
      topo->assign_age_nodeLabel_popSize(eventT_inc,popSize, splittingTime);
      long double logProb = topo->compute_logProb_zeroMig(eventT_inc,popSize, splittingTime);

      logPseudoProb += logProb;
    }

  // std::cout << "logProb = " << logProb <<"\n";
  
  #ifdef DEBUG
  // std::cout << "Exiting Chain::compute_logConditionalProb_zeroMig_subtrees()\n";
#endif //DEBUG
  return logPseudoProb;
}
