/* MIST copyright 2016 by Yujin Chung and Jody Hey */

#include <cmath>
#include <iostream>
#include <fstream>
#include "Chain.hpp"



void Chain::print_listTrees()
{
	unsigned int list_size = list_trees.size();
	if(list_size==0)
	{
		std::cout << "No saved tree topologies\n";
	}
	else
	{
		std::cout << "List of sampled tree topologies:\n";
		std::cout << "ID\tTopology\n";
		for(unsigned int i=0; i<list_size; i++)
		{
			std::cout <<i<< "\t";
			list_trees.at(i)->print_topo();
		}

	}
}


void Chain::fileprint_listTrees(std::ofstream& fp)
{
	unsigned int list_size = list_trees.size();
	if(list_size==0)
	{
		fp << "No saved tree topologies\n";
	}
	else
	{
		fp << "List of sampled tree topologies:\n";
		fp << "ID\tTopology\n";
		for(unsigned int i=0; i<list_size; i++)
		{
			fp <<i<< "\t";
			list_trees.at(i)->fileprint_topo(fp);
		}

	}
}


void Chain::read_LmodeInputFile( IM im)
{
  if(im.get_newickTreeFormat() ==1)
    {     
      // read_newickTrees(im);
      read_newickTrees_partial(im);

      // save the list of ranked tree topologies
      if(cpuID ==0)
	{
	  // print_listTrees();
	  saveListTopo2File();
	}
    }
  else
    {
      // read_logJointPriorTrees_partial();
      read_logEachLikelihood_partial(); // logLikelihood should be read first.
      if(im.get_priorType() == 1)
	read_eachLogPrior_partial(); 
      // read_logJointLikelihood_partial(); 
      read_treeIDs_partial();
      read_coalTimes_partial(im);
      // read_mutationScaler_partial();
      read_tipIDs_partial(im);
      read_listTrees();
    }

  find_print_theMaxHeightsTrees();
  
  return;
}

void Chain::find_print_theMaxHeightsTrees()
{
  #ifdef DEBUG
  // std::cout <<"nSubSample = "<< nSubSample << " numSubLoci = "<< numSubLoci <<" sizeForest = "<<sizeForest <<"\n";
#endif //DEBUG
  
  double max_age=0;
  double max_test =0;
  for(unsigned int i=0; i<nSubSample; i++)
    for(unsigned int j=0; j<numSubLoci; j++)
      {
	if(Forest ==0||Forest ==2)
	  {
	    double age = coalTimes.at(i).at(j).front();
	    // std::cout <<"cpuID =" << cpuID  << " age = " << age <<"\n";
	    if(max_age < age)
	      {
		max_age = age;
	      }
	  }
	else if(Forest==1||Forest==3)
	  {
	    for(unsigned int k=0; k< sizeForest; k++)
	      {
		double age = subcoalTimes.at(i).at(j).at(k).front();
		if(max_age < age)	      
		  max_age = age;
		/*
		double age_test = *max_element(subcoalTimes.at(i).at(j).at(k).begin(), subcoalTimes.at(i).at(j).at(k).end());
		if(max_test < age_test)
		  max_test = age_test;
		*/
		#ifdef DEBUG
		/*
		std::cout <<"cpuID = "<<cpuID;
		std::cout << " i="<<i <<" j="<<j <<" k="<<k<< " height(age) = "<<Coalsubtrees.at(i).at(j).at(k)->get_age()<<"\n";
		std::cout <<"age = " <<age <<" age_test = "<< age_test <<"\n";
		*/
#endif //DEBUG
	      }
	  }
	    /*
	    std::cout <<"cpuID = " << cpuID <<" max_age=" 
		      << max_age <<" age = " << age 
		      << "locusID = " << locusID_start+j 
		      <<" " ;
	    trees_atPrev.at(j)->print_coaltree();
	    */	    
	  
      }

  double max_age_global = 0;
  MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Reduce(&max_age, &max_age_global, 1, MPI_DOUBLE, MPI_MAX, 0);

  if(cpuID ==0)
    {
      std::cout <<"\nThe maximum of tree heigths is " << max_age_global <<". If the upper bound of the splitting time is higher than this maximum, please use an upper bound smaller than the maximum height.\n\n";
      //  std::cout <<"\nThe maximum of tree heigths is " << max_test <<". If the upper bound of the splitting time is higher than this maximum, please use an upper bound smaller than the maximum height.\n\n";
    }

  return;
}



void Chain::read_newickTrees_partial( IM im)
{
  char* inputFile = im.get_newickTreeFileName();
  char* SeqPopFile = im.get_SeqPopFileName();
  if(cpuID ==0)
    std::cout << "\n Reading "<< inputFile <<"\n";

  #ifdef DEBUG
  /*
  std::cout <<"numSubLoci = "<< numSubLoci 
	    <<" nSubSample = "<< nSubSample <<"\n";
  */
#endif //DEBUG

  trees_atPrev.resize(numSubLoci); 
  logPrior_trees.resize(nSubSample);
  logPriorTrees_atPrev = 0.0;
  if(Forest == 0||Forest ==2)
    {
      treeIDs.resize(nSubSample);
      coalTimes.resize(nSubSample);
      tipIDs.resize(nSubSample);
    }
  else if(Forest == 1||Forest ==3)
    {
      Coalsubtrees.resize(nSubSample);
      Coalsubtrees.at(0).resize(numSubLoci);
      subtreeIDs.resize(nSubSample);
      subcoalTimes.resize(nSubSample);
      subtipIDs.resize(nSubSample);
    }

  ifstream inFile; 

  // seq-population names
  unsigned int word_counter = 0;
  unsigned int label = 0;
  inFile.open(SeqPopFile);
  while(inFile >> label){
    SeqPop.push_back(label);
    word_counter++;
  }
  inFile.close();

  #ifdef DEBUG
  /*
  std::cout <<"word_counter = " << word_counter <<"\n";
  for(unsigned int i=0; i<word_counter; i++)
    std::cout << SeqPop.at(i) <<" ";
  std::cout <<"\n";
  */
#endif //DEBUG
  
  string word;
  unsigned int iter = 0;
  unsigned int locusID = 0;
  unsigned int line_counter = 1;
  std::string tree_string;

  inFile.open(inputFile);
  
  unsigned int count_nSample_perProcess = 0;
  while(locusID < n_loci)
    {
      inFile >> locusID;
      inFile >> tree_string;

      if(locusID >=locusID_start && locusID <= locusID_end)
	{
	  #ifdef DEBUG
	  /*
	  std::cout << "n_loci = " << n_loci 
		    << " locusID = " << locusID <<"\n";
	    std::cout << " tree_string = " << tree_string <<"\n";
	    */
#endif //DEBUG
	  
	  node* tree = new node;
	  tree->convertFromNewick(tree_string,1);
	  tree->assignPopulations2Tips(SeqPop);
	  trees_atPrev.at(locusID-locusID_start) = tree;
	  /*
	  if(locusID == locusID_end)
	    {
	      collectAllUpdates_Lmode(0);
	    }	        
	  */

	  if(Forest == 1||Forest ==3)
	    {
	      if(sizeCoalsubtree == 0)
		{
		  std::cout << "\n**** Error in Chain::read_newickTrees_partial( IM im) ***\n";
		  std::cout <<" sizeCoalsubtree is zero.\n";
		}
	      else
		{
		  if(locusID == locusID_start)
		    numTotalSeq = trees_atPrev.at(0)->size_tree();
		  getSubtree(sizeCoalsubtree,trees_atPrev.at(locusID-locusID_start), 0, locusID-locusID_start);
		}
	    }	  
	  
	}
      if(locusID==n_loci-1)
	locusID++;
    }

  inFile.close();
  #ifdef DEBUG
  // std::cout <<"Done with reading the file.\n";
#endif //DEBUG
  
  if(Forest ==0 || Forest == 2)
    {
      if(cpuID == 0)
	collectAllUpdates_Lmode(0);
      MPI::COMM_WORLD.Barrier();
      collectAllUpdates_bwProcs_Lmode();
    }
  else if(Forest == 1||Forest ==3) // forest of subtrees
    {
      if(cpuID == 0)
	collectAllUpdates_Lmode_subtrees(0);
      MPI::COMM_WORLD.Barrier();
      collectAllUpdates_bwProcs_Lmode_subtrees();
    }
  else
    std::cout <<"\n*** Error *** Forest is "<<Forest <<".\n";

  


  
  // REMOVE
  /*
  if(cpuID==0)
    {
      unsigned int list_size = list_trees.size();
      for(unsigned int i=0; i<list_size; i++)
	{
	  std::cout <<i <<" ";
	  list_trees.at(i)->print_topo();
	}
    }
  */
  
  return;
}


void Chain::collectAllUpdates_bwProcs_Lmode()
{  
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

	  collectAllUpdates_Lmode(0);
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

void Chain::read_newickTrees(IM im)
{
  //std::cout << "\nIn Chain::read_newickTrees()\n";
  char* inputFile = im.get_newickTreeFileName();
  std::cout << "\n Reading "<< inputFile <<"\n";

  try
    {
      trees_atPrev.resize(n_loci);
      treeIDs.resize(n_MCMCgen);
      coalTimes.resize(n_MCMCgen);
      tipIDs.resize(n_MCMCgen);
      logPrior_trees.resize(n_MCMCgen);
      logPriorTrees_atPrev = 0.0;
    }
  catch(std::exception &e) 
    {
      std::cout << "Error In Chain::read_newickTrees()\n";
    }

  string word;
  unsigned int iter = 0;
  unsigned int locusID = 0;
  unsigned int line_counter = 1;
  std::string tree_string;

  ifstream inFile; 
  inFile.open(inputFile);
  /*
  for(unsigned int i=0; i<3; i++)
    {
      inFile >> word;
    }
  */
  
  unsigned int count_nSample_perProcess = 0;
  while(locusID < n_loci)
    {
      inFile >> locusID;
      inFile >> tree_string;
      std::cout << "locusID = " << locusID  << " tree_string = " << tree_string <<"\n";
     
      node* tree = new node;
      tree->convertFromNewick(tree_string,1);
      tree->assignPopulations2Tips(im.getLoci()[0]);
      trees_atPrev.at(locusID) = tree;
      	 
      if(locusID == n_loci-1)
	{
	  collectAllUpdates(count_nSample_perProcess, 1);	   
	  locusID++;
	  // count_nSample_perProcess++;   
	}	        
    }
     
  inFile.close();
  
  return;
}



void Chain::read_logJointPriorTrees_partial()
{
  ifstream inFile;
  inFile.open("MCMCsample_logJointPriorTrees.txt");
  
  string word;
  double logP = 0.0;
  unsigned int iter = 0;
  // unsigned int line_counter = 1;
	
  inFile >> word;
  inFile >> word;

  while(inFile >> iter &&  iter <= sampleID_end)
    {
      inFile >> logP;
      // if(line_counter >= 1)
      if(iter >= sampleID_start && iter <= sampleID_end)
	{
	  logPrior_trees.push_back(logP);
	}
      // line_counter++;
    }  
  inFile.close();


  return;
}

// YC 1/6/2016
// This function should be called 
// if the prior for trees in MCMC was NOT uniform (improper).
void Chain::read_eachLogPrior_partial()
{
  ifstream inFile;
  inFile.open("MCMCsample_eachLogPrior.txt");
  
  string word;
  double logPrior = 0.0;
  unsigned int locusID = 0;
  unsigned int iter = 0;
  std::vector<double> logPrior_iter;

  inFile >> word;
  inFile >> word;
  inFile >> word;

  // REMOVE
  // std::cout << "sampleID_start = " << sampleID_start <<" and sampleID_end = " << sampleID_end <<"\n";

  while(inFile >> iter &&  iter <= sampleID_end)
    {
      inFile >> locusID;
      inFile >> logPrior;
      if(iter >= sampleID_start && iter <= sampleID_end)
	{
	  if(locusID >= locusID_start && locusID <=locusID_end)
	    {
	      logPrior_iter.push_back(logPrior);
	    }
	  if(locusID ==locusID_end)
	    {
	      logPrior_each.push_back(logPrior_iter);
	      logPrior_iter.resize(0);
	    }
	}
    }  
  inFile.close();

  return;
}




// YC 7/17/2015
// This function is always called, but it
// needs to be called only if "TreeWithMaxP=1".
void Chain::read_logEachLikelihood_partial()
{
  ifstream inFile;
  inFile.open("MCMCsample_logEachLikelihoods.txt");
  
  string word;
  double logLik = 0.0;
  unsigned int locusID = 0;
  unsigned int iter = 0;
  std::vector<double> loglik_iter;

  inFile >> word;
  inFile >> word;
  inFile >> word;

  // REMOVE
  // std::cout << "sampleID_start = " << sampleID_start <<" and sampleID_end = " << sampleID_end <<"\n";

  while(inFile >> iter &&  iter <= sampleID_end)
    {
      inFile >> locusID;
      inFile >> logLik;
      if(iter >= sampleID_start && iter <= sampleID_end)
	{
	  if(locusID >= locusID_start && locusID <=locusID_end)
	    {
	      if(TreesWithMaxP==1)
		{
		  // The lengths of MaxLogLik and MaxSampleID are the same as numSubLoci
		  if(iter == sampleID_start)
		    {
		      MaxLogLik.push_back(logLik);
		      MaxSampleID.push_back(iter);
		    }
		  else
		    {
		      if(logLik > MaxLogLik.at(locusID-locusID_start))
			{
			  MaxLogLik.at(locusID-locusID_start) = logLik;
			  MaxSampleID.at(locusID-locusID_start) = iter;
			}
		    }
		}
	      loglik_iter.push_back(logLik);
	    }
	  if(locusID ==locusID_end)
	    {
	      logLikelihood.push_back(loglik_iter);
	      loglik_iter.resize(0);
	    }
	}
    }  
  inFile.close();

  if(TreesWithMaxP==1)
    {
      locusID = locusID_start;
      std::cout << "SampleIDs with largest likelihoods\n";
      for(locusID = locusID_start; locusID <=locusID_end; locusID++)
	{
	  std::cout << "locusID = " << locusID << ": " << MaxSampleID.at(locusID-locusID_start) 
		    << "( " << MaxLogLik.at(locusID-locusID_start) <<")  ";
	}
      std::cout << "\n";
    }

  return;
}




void Chain::read_logJointLikelihood_partial()
{
  ifstream inFile;
  inFile.open("MCMCsample_logJointLikelihoods.txt");
  
  string word;
  double logJointLik = 0.0;
  unsigned int iter = 0;
	
  inFile >> word;
  inFile >> word;

  while(inFile >> iter &&  iter <= sampleID_end)
    {
      inFile >> logJointLik;
      if(iter >= sampleID_start && iter <= sampleID_end)
	{
	  logJointLikelihood.push_back(logJointLik);
	}
    }  
  inFile.close();


  return;
}





void Chain::read_mutationScaler_partial()
{
  ifstream inFile;
  inFile.open("mutationScaler_MCMCsample.txt");
  
  string word;
  unsigned int iter = 0;
  unsigned int locusID = 0;
  double mt = 0;
  unsigned int line_counter = 1;
  std::vector<double> mt_iter;

  for(unsigned int i=0; i<3; i++)
    {
      inFile >> word;
    }

  
  while(inFile >> iter &&  iter <= sampleID_end)
    {
      inFile >> locusID;
      inFile >> mt;
      if(iter >= sampleID_start && iter <= sampleID_end)
	{
	  mt_iter.push_back(mt);		 
	  if(locusID == n_loci-1)
	    {
	      mutationScaler.push_back(mt_iter);
	      mt_iter.resize(0);
	      
	    }	  
	}
    } 
  inFile.close();

   return;
}


void Chain::read_kappa()
{
  ifstream inFile;
  inFile.open("kappa_MCMCsample.txt");
  
  string word;
  unsigned int iter = 0;
  unsigned int locusID = 0;
  double kp = 0;
  unsigned int line_counter = 1;
  std::vector<double> kappa_iter;

  for(unsigned int i=0; i<3; i++)
    {
      inFile >> word;
    }

  while(inFile >> iter)
    {
      if(line_counter == 1)
	{
	  //inFile >> iter;
	  inFile >> locusID;
	  inFile >> kp;
	  kappa_iter.push_back(kp);
	}
      else
	{     
	  inFile >> locusID;
	  inFile >> kp;	  
	  if(iter == kappa.size())
	    {
	      kappa_iter.push_back(kp);
	    }
	  else
	    {
	      kappa.push_back(kappa_iter);
	      kappa_iter.resize(0);
	      kappa_iter.push_back(kp);
	    }
	}
      line_counter++;
    } 
  kappa.push_back(kappa_iter); 
  inFile.close();

 

  return;
}




void Chain::read_treeIDs_partial()
{ 

  ifstream inFile;
  inFile.open("MCMCsample_treeIDs.txt");
  
  string word;
  unsigned int iter = 0;
  unsigned int locusID = 0;
  unsigned int trID = 0;
  // unsigned int line_counter = 1;
  std::vector<unsigned int> treeIDs_iter;

  for(unsigned int i=0; i<3; i++)
    {
      inFile >> word;
    }

  // std::cout << "sampleID_start = " << sampleID_start <<"\n";
  // std::cout << "sampleID_end = " << sampleID_end <<"\n";
  // std::cout << "locusID_end = " << locusID_end <<"\n";

  while(inFile >> iter && iter <= sampleID_end)
    {
      inFile >> locusID;
      inFile >> trID;
      if(iter >= sampleID_start && iter <= sampleID_end)
	{
	  if(locusID >= locusID_start && locusID <=locusID_end)
	    treeIDs_iter.push_back(trID);		 
	  if(locusID ==locusID_end)
	    {
	      treeIDs.push_back(treeIDs_iter);
	      // std::cout << "treeIDs_iter.size() = " <<treeIDs_iter.size() <<"\n";
	      treeIDs_iter.resize(0);
	    }
	}
    } 
  // treeIDs.push_back(treeIDs_iter); 
  inFile.close();

  return;
}



void Chain::read_treeIDs()
{
  ifstream inFile;
  inFile.open("MCMCsample_treeIDs.txt");
  
  string word;
  unsigned int iter = 0;
  unsigned int locusID = 0;
  unsigned int trID = 0;
  unsigned int line_counter = 1;
  std::vector<unsigned int> treeIDs_iter;

  for(unsigned int i=0; i<3; i++)
    {
      inFile >> word;
    }

  while(inFile >> iter)
    {
      if(line_counter == 1)
	{
	  //inFile >> iter;
	  inFile >> locusID;
	  inFile >> trID;
	  treeIDs_iter.push_back(trID);
	}
      else
	{
	  //inFile >> iter;
	  inFile >> locusID;
	  inFile >> trID;
	  if(iter == treeIDs.size())
	    {
	      treeIDs_iter.push_back(trID);
	    }
	  else
	    {
	      treeIDs.push_back(treeIDs_iter);
	      treeIDs_iter.resize(0);
	      treeIDs_iter.push_back(trID);
	    }
	}
      line_counter++;
    } 
  treeIDs.push_back(treeIDs_iter); 
  inFile.close();


  return;
}




void Chain::read_coalTimes(IM im)
{
  unsigned int nloci = im.get_nGeneCopies();
 
  ifstream inFile;
  inFile.open("coalescentTimes.txt");
  
  string word;
  double coalT =0.0;
  unsigned int iter = 0;
  unsigned int locusID = 0;
  unsigned int line_counter = 1;
  std::list<double> coalT_locus;
  std::vector<std::list<double> >  coalT_iter;

  for(unsigned int i=0; i< nloci+1; i++)
    {
      inFile >> word;
    }

  while(inFile >> iter)
    {
      if(line_counter == 1)
	{
	  //inFile >> iter;
	  inFile >> locusID;
	  for(unsigned int i=0; i< nloci-1; i++)
	    {
	      inFile >> coalT;
	      coalT_locus.push_back(coalT);
	    }
	  coalT_iter.push_back(coalT_locus);
	  coalT_locus.resize(0);
	}
      else
	{
	  //inFile >> iter;
	  inFile >> locusID;
	  for(unsigned int i=0; i< nloci-1; i++)
	    {
	      inFile >> coalT;
	      coalT_locus.push_back(coalT);
	    }
	  if(iter == coalTimes.size())
	    {
	      coalT_iter.push_back(coalT_locus);
	      coalT_locus.resize(0);
	    }
	  else
	    {  
	      coalTimes.push_back(coalT_iter);
	      coalT_iter.resize(0);
	      coalT_iter.push_back(coalT_locus);
	      coalT_locus.resize(0);
	    }
	}
      line_counter++;
    } 
  coalTimes.push_back(coalT_iter); 
  inFile.close();

  return;
}





void Chain::read_coalTimes_partial(IM im)
{

  ifstream inFile;
  inFile.open("coalescentTimes.txt");
  
  unsigned nCoalEvents = im.get_nGeneCopies() -1;
 
  string word;
  double coalT =0.0;
  unsigned int iter = 0;
  unsigned int locusID = 0;
  // unsigned int line_counter = 1;
  std::list<double> coalT_locus;
  std::vector<std::list<double> >  coalT_iter;

  for(unsigned int i=0; i< nCoalEvents+2; i++)
    {
      inFile >> word;
    }

  while(inFile >> iter && iter <= sampleID_end)
    {
      inFile >> locusID;
      for(unsigned int i=0; i< nCoalEvents; i++)
	{
	  inFile >> coalT;
	  if(iter >= sampleID_start && iter <= sampleID_end
	     && locusID >= locusID_start && locusID <= locusID_end)
	    coalT_locus.push_back(coalT);
	}
      if(iter >= sampleID_start && iter <= sampleID_end)
	{
	  if(locusID >= locusID_start && locusID <= locusID_end)
	    coalT_iter.push_back(coalT_locus);
	  if(locusID == locusID_end)
	    {
	      coalTimes.push_back(coalT_iter);	      
	      coalT_iter.resize(0);
	    }    
	}
      coalT_locus.resize(0);  
    } 
  inFile.close();

  return;
}




void Chain::read_tipIDs_partial(IM im)
{
  unsigned int nGeneCopies = im.get_nGeneCopies();
 
  ifstream inFile;
  inFile.open("MCMCsample_tipIDs.txt");
  
  string word;
  unsigned int id =0;
  unsigned int iter = 0;
  unsigned int locusID = 0;
  // unsigned int line_counter = 1;
  std::vector<unsigned int> id_locus;
  std::vector<std::vector<unsigned int> >  id_iter;

  for(unsigned int i=0; i< nGeneCopies+2; i++)
    {
      inFile >> word;
    }

  while(inFile >> iter && iter <= sampleID_end)
    {
      inFile >> locusID;
      for(unsigned int i=0; i< nGeneCopies; i++)
	{
	  inFile >> id;
	  if(iter >= sampleID_start && iter <= sampleID_end)
	    id_locus.push_back(id);
	}
      if(iter >= sampleID_start && iter <= sampleID_end)
	id_iter.push_back(id_locus);
      id_locus.resize(0);
      if(locusID < n_loci)
	{
	  if(iter >= sampleID_start && iter <= sampleID_end)
	    tipIDs.push_back(id_iter);
	  id_iter.resize(0);	  
	}    
    } 
  // tipIDs.push_back(id_iter); 
  inFile.close();

  return;
}





void Chain::read_listTrees()
{
  ifstream file;
  file.open ("ListTopo.txt");
  string word;

  file >> word;
  file >> word;

  unsigned int iter =0;
  std::string topo_string;
  unsigned int line_counter = 1;
  while(file >> iter)
    {
      file >> topo_string;
      nodeSimple* topo = new nodeSimple;
      topo->convertFromNewick(topo_string,1);
      topo->computeSizes();
      list_trees.push_back(topo);
      line_counter++;

    }  
  file.close();


  return;
}




void Chain::saveKappa2File()
{
	ofstream file;
	file.open ("kappa_MCMCsample.txt");
	file << "iter\tlocus\tkappa\n";
	for(unsigned int i=0; i<kappa.size(); i++)
		for(unsigned int j=0; j<kappa.at(i).size(); j++)
		{
			file << i <<"\t" <<j <<"\t" << kappa.at(i).at(j) <<"\n";
		}
	file.close();
}


void saveKappa2File(std::vector<std::vector<double> > kappa)
{
	ofstream file;
	file.open ("kappa_MCMCsample.txt");
	file << "iter\tlocus\tkappa\n";
	for(unsigned int i=0; i<kappa.size(); i++)
		for(unsigned int j=0; j<kappa.at(i).size(); j++)
		{
			file << i <<"\t" <<j <<"\t" << kappa.at(i).at(j) <<"\n";
		}
	file.close();
}

void Chain::saveMutationScalers2File()
{
	ofstream file;
	file.open ("mutationScaler_MCMCsample.txt");
	file << "iter\tlocus\tmutationScaler\n";
	for(unsigned int i=0; i<mutationScaler.size(); i++)
		for(unsigned int j=0; j<mutationScaler.at(i).size(); j++)
		{
			file << i <<"\t" <<j <<"\t" << mutationScaler.at(i).at(j) <<"\n";
		}
	file.close();
}

void Chain::saveTipIDs()
{
  ofstream file;
  file.open ("MCMCsample_tipIDs.txt");
  file << "iter\tlocus";
  unsigned int nTips = tipIDs.at(0).at(0).size();
  for(unsigned int i=0; i<nTips; i++)
    {
      file << "\ttip" << i;
    }
  file << "\n";
  for(unsigned int i=0; i<tipIDs.size(); i++)
    {
      for(unsigned int j=0; j<tipIDs.at(i).size(); j++)
	{
	  file << i <<"\t" <<j;
	  for(unsigned int k=0; k<nTips; k++)
	    {	    
	      file << "\t" << tipIDs.at(i).at(j).at(k);
	    }
	  file << "\n";
	}
    }
  file.close();
}

void Chain::saveTreeIDs()
{
  ofstream file;
  file.open ("MCMCsample_treeIDs.txt");
  unsigned int nloci = treeIDs.at(0).size();

  file << "iter\tlocus\ttreeID\n";
  for(unsigned int i=0; i<treeIDs.size(); i++)
    {
      for(unsigned int j=0; j<nloci; j++)
	{
	  file << i <<"\t" <<j <<"\t" << treeIDs.at(i).at(j) << "\n";
	}
    }
  file.close();
}


// old version - YC 6/7/2014
void saveMutationScalers2File(std::vector<std::vector<double> > mutationScaler)
{
	ofstream file;
	file.open ("mutationScaler_MCMCsample.txt");
	file << "iter\tlocus\tmutationScaler\n";
	for(unsigned int i=0; i<mutationScaler.size(); i++)
		for(unsigned int j=0; j<mutationScaler.at(i).size(); j++)
		{
			file << i <<"\t" <<j <<"\t" << mutationScaler.at(i).at(j) <<"\n";
		}
	file.close();
}



void Chain::initialize_trees_atPrev(vector<locus> lc)
{
  // std::cout << "In Chain::initialize_trees_atPrev\n";

  trees_atPrev.resize(0);

  for(unsigned int i=0; i< numSubLoci; i++)
    {
      // std::cout << "locus = " << i << "\n";
      if(i!=0)
	lc.at(i).initialize_multiLocusSpecific_mutationRate(lc.at(0).get_multiLocusSpecific_mutationRate());

      unsigned int modelType = lc.at(i).getLikelihoodModel();
      if(modelType == 0) // IS model
	{
	  trees_atPrev.push_back(initial_coalTree_IS(lc.at(i)) );
	}
      else if(modelType == 2) // HKY
	{
	  trees_atPrev.push_back(initial_coalTree_HKY(lc.at(i)) );
	}
      else
	{
	  std::cout << "\n*** Error in  Chain::initialize_trees_atPrev() ***\n";
	  std::cout << "*** The mutation/substitution model is not properly assigned.\n\n";
	}
    }
  return;
}

/**
 * \param n_lc Number of loci
 * \param loci vector of type locus
 * \return NULL
 */
void Chain::InitializeChain(IM im, int nMCMCgen, int c_id, unsigned int processID, unsigned int n_chains, unsigned int numprocesses)
{

  // std::cout << "In Chain::InitializeChain\n";
  
  nChains = n_chains;
  nProcesses = numprocesses;
  chainID = c_id;
  n_loci = im.get_nLoci();
  n_MCMCgen = nMCMCgen;
  popSizeMax = im.get_popSizeMax();
  vector<locus> loci = im.getLoci();
  if(lociInParallel ==1)
    {
      vector<locus> loci_temp;
      for(unsigned int i=locusID_start; i <= locusID_end; i++)
	loci_temp.push_back(loci.at(i));
      loci_temp.at(0).initialize_multiLocusSpecific_mutationRate(loci.at(0).get_multiLocusSpecific_mutationRate());
      loci.resize(numSubLoci);
      loci = loci_temp;      
    }

  priorType = im.get_priorType();
  changeOfRatePoint = im.get_changeOfRatePoint();

  // MCMC savings
  coalTimes.resize(n_MCMCgen);
  tipIDs.resize(n_MCMCgen);
  treeIDs.resize(n_MCMCgen);
  kappa.resize(n_MCMCgen);
  mutationScaler.resize(n_MCMCgen);
  logLikelihood.resize(n_MCMCgen);
  logPrior_trees.resize(n_MCMCgen);
  logJointLikelihood.resize(n_MCMCgen);
  logPrior_each.resize(n_MCMCgen);

  kappa_atPrev.resize(numSubLoci);
  mutationScaler_atPrev.resize(numSubLoci);
  logLikelihood_atPrev.resize(numSubLoci);
  if(lociInParallel ==1) // independent coalescent trees
    logPriorTrees_atPrev_lociParallel.resize(numSubLoci);

  // Added by YC 8/7/2014
  numTries_trees.resize(1,numSubLoci);
  numAccpt_trees.resize(1,numSubLoci);
  numTries_trees.setZero();
  numAccpt_trees.setZero();
  
  //--- Build initial trees ---//
  initialize_trees_atPrev(loci); // 'loci' contains a part of the full data

  //trees_gen.at(0).resize(loci.size());
  //	std::cout << "The initial trees:\n";
  /*
  for(unsigned int i=0; i< numSubLoci; i++)
    {
      std::cout << "Process "<< processID << " Chain "<< chainID <<" ";
      std::cout << "Locus " << i+locusID_start <<": ";
      try{	
	trees_atPrev.at(i)->print_coaltree();
	//coalTimes_atPrev.at(i) = trees_atPrev.at(i)->get_coalescentTimes();
      } catch (std::exception &e) {
	std::cout << "Error in Chain::InitializeChain\n";
	std::cout << "numSubLoci = " << numSubLoci << " but i = " << i <<"\n";
	std::cout << "Can't access elements of trees_gen - index is out of bounds\n";
      }
      
      
    }
  std::cout << "\n";
  */
  
  //std::vector<double> kappa_initial;
  //std::vector<double> mutScaler_initial;
  for(unsigned int i=0; i<numSubLoci; i++)
    {
      kappa_atPrev.at(i) = 2.0; 

	//std::cout << "processID = " << processID 
	//   << " loci.at(0) mutRate = " << loci.at(0).get_multiLocusSpecific_mutationRate() <<"\n";
       if(loci.at(0).get_multiLocusSpecific_mutationRate()==1)
         {
            mutationScaler_atPrev.at(i) = loci.at(i).get_nSites();
	   // std::cout << "i= " <<i <<"locusID_start = " << locusID_start 
           //           << " nSites = " << loci.at(i).get_n_sites_uniq() <<"\n";
         }
       else
         {   
            /// All mutation rate scalers are set to 1.
            mutationScaler_atPrev.at(i) = 1.0;
         }
    }

  // std::cout << "computing the likelihoods..\n";
  for (unsigned int i = 0; i < numSubLoci; i++)
    {
      // std::cout << "locus id = " << i << "\n";
      try 
	{
	  unsigned int likelihoodModel = loci.at(i).getLikelihoodModel();
	  if(likelihoodModel == 0) // IS
	    logLikelihood_atPrev.at(i) = trees_atPrev.at(i)->logLikelihood_IS(loci.at(i), mutationScaler_atPrev.at(i),0);
	  else if(likelihoodModel == 2) // HKY	
	    logLikelihood_atPrev.at(i) = trees_atPrev.at(i)->logLikelihood_HKY(loci.at(i),mutationScaler_atPrev.at(i), kappa_atPrev.at(i));
	  
	} catch (std::exception &e) {
	std::cout << "Error in Chain::InitializeChain\n";
	std::cout << "numSubLoci = " << numSubLoci << " but i = " << i <<"\n";
	std::cout << "Can't access elements of loci - index is out of bounds\n";
      }
    }
  
  // std::cout << "computing the priors..\n";
  if(lociInParallel ==1)
    {
      logPriorTrees_atPrev =1;
      for (unsigned int i = 0; i < numSubLoci; i++)
	{
	  logPriorTrees_atPrev_lociParallel.at(i) 
	    = trees_atPrev.at(i)->compute_prior_indepLoci(loci.at(i), popSizeMax,priorType,changeOfRatePoint);
	  logPriorTrees_atPrev *=logPriorTrees_atPrev_lociParallel.at(i) ;
	}
    }
  else
    {
      compute_totalnum_coalEvents(loci);
      logPriorTrees_atPrev =compute_logPrior_usingTrees_atPrev();
    }
  

  // Make a list of IDs of pairs of mutation scalars to update
  if(nPairsMut == 0)
    {
      nPairsMut = static_cast<unsigned int>( floor(n_loci/2));
    }
  // std::cout << "nPairsMut  = " << nPairsMut <<"\n";

  /*
  unsigned int totalNumIter = im.get_nBurning() + im.get_thinning()*n_MCMCgen;
  pairIDsMut.resize(totalNumIter);
  std::vector<unsigned int> IDs; 
  std::vector<unsigned int> randomIDs; 
  for(unsigned int ll =0; ll<n_loci; ll++)
    IDs.push_back(ll);
  // std::cout << "n_loci = " << n_loci
  //	    << "IDs.size() = " << IDs.size() <<"\n";

  for(unsigned int i =0; i<totalNumIter; i++)
    {
      // std::cout << "i = " << i <<"\n";
      randomIDs = IDs;
      // std::cout << "randomIDs.size() = " << randomIDs.size() <<"\n";
      std::random_shuffle(randomIDs.begin(),randomIDs.end()); // shuffle elemenets in IDs
      pairIDsMut.at(i).resize(2*nPairsMut);
      for(unsigned int j=0; j< 2*nPairsMut; j++)	
	{	  
	  pairIDsMut.at(i).at(j) = randomIDs.at(j);	
	  // std::cout << " randomIDs.at(j) = " << randomIDs.at(j) <<"\n";
	}
    }

  // Share "pariIDsMut" on the CPU with ID 0
  unsigned int id =0;
  for(unsigned int i=0; i<totalNumIter; i++)
    {
      for(unsigned int j=0; j< 2*nPairsMut; j++)
	{
	  if(processID==0)
	    id = pairIDsMut.at(i).at(j);
	  MPI::COMM_WORLD.Barrier();
	  MPI::COMM_WORLD.Bcast(&id, 1, MPI::UNSIGNED, 0);
	  MPI::COMM_WORLD.Barrier();
	  if(processID != 0)
	    pairIDsMut.at(i).at(j) = id;
	}
    }
  */

  // std::cout << "END of Chain::InitializeChain\n";

  return;
} /* End of function Chain::InitializeChain() */





/// Added by YC 6/6/2014
/**
 * \brief Chain::InitializeChain
 * \param n_lc Number of loci
 * \param loci vector of type locus
 * \return NULL
 */
void Chain::InitializeChain_usePriorOnly(IM im, int nMCMCgen, int c_id, unsigned int processID, unsigned int n_chains, unsigned int numprocesses)
{
  nChains = n_chains;
  nProcesses = numprocesses;
  chainID = c_id;
  n_loci = im.get_nLoci();
  n_MCMCgen = nMCMCgen;
  popSizeMax = im.get_popSizeMax();
  vector<locus> loci = im.getLoci();
  coalTimes.resize(n_MCMCgen);
  tipIDs.resize(n_MCMCgen);
  treeIDs.resize(n_MCMCgen);
  kappa.resize(n_MCMCgen);
  kappa_atPrev.resize(n_loci);
  mutationScaler.resize(n_MCMCgen);
  mutationScaler_atPrev.resize(n_loci);
  logLikelihood.resize(n_MCMCgen);
  logLikelihood_atPrev.resize(n_loci);
  logPrior_trees.resize(n_MCMCgen);

  // Added by YC 8/7/2014
  numTries_trees.resize(1,n_loci);
  numAccpt_trees.resize(1,n_loci);
  numTries_trees.setZero();
  numAccpt_trees.setZero();
  
  //--- Build initial trees ---//
  initialize_trees_atPrev(loci);
  //trees_gen.at(0).resize(loci.size());
  //	std::cout << "The initial trees:\n";
  for(unsigned int i=0; i<n_loci; i++)
    {
      // std::cout << "Process "<< processID << " Chain "<< chainID <<" ";
      // std::cout << "Locus " << i <<": ";
      try{
	
	trees_atPrev.at(i)->print_coaltree();
	//coalTimes_atPrev.at(i) = trees_atPrev.at(i)->get_coalescentTimes();
      } catch (std::exception &e) {
	std::cout << "Can't access elements of trees_gen - index is out of bounds\n";
      }
      
      
    }
  std::cout << "\n";

  
  for(unsigned int i=0; i<n_loci; i++)
    {
      kappa_atPrev.at(i) = 2.0;
      
      /// All mutation rate scalers are set to 1.
      mutationScaler_atPrev.at(i) = 1.0;
      //mutScaler_initial.push_back(1.0);
    }

  for (unsigned int i = 0; i < n_loci; i++) 
    {
      try {
	logLikelihood_atPrev.at(i) = 0.0;//trees_atPrev.at(i)->logLikelihood_HKY(loci.at(i),mutationScaler_atPrev.at(i), kappa_atPrev.at(i));
	
      } catch (std::exception &e) {
	std::cout << "Can't access elements of loci - index is out of bounds\n";
      }
    }
  
  compute_totalnum_coalEvents(loci);
  logPriorTrees_atPrev =compute_logPrior_usingTrees_atPrev();
  
  return;
} /* End of function Chain::InitializeChain() */










void Chain::delete_temporaryPara()
{
	for(unsigned int i=0; i<n_loci; i++)
	{
		trees_atPrev.at(i)->deleteCoalTree();
		//delete trees_atPrev.at(i);
	}
	trees_atPrev.resize(0);
	// coalTimes_atPrev.resize(0);
	kappa_atPrev.resize(0);
	mutationScaler_atPrev.resize(0);
	logLikelihood_atPrev.resize(0);
}

void Chain::deleteAll()
{
	deleteTrees();

	coalTimes.resize(0);
	// coalTimes_atPrev.resize(0);
	kappa.resize(0);
	kappa_atPrev.resize(0);
	mutationScaler.resize(0);
	mutationScaler_atPrev.resize(0);
	logLikelihood.resize(0);
	logLikelihood_atPrev.resize(0);
	logPrior_trees.resize(0);
}

void Chain::deleteAll_Lmode()
{
  //deleteTrees();
  deleteListTrees();
  coalTimes.resize(0);
  kappa.resize(0);
	// kappa_atPrev.resize(0);
  mutationScaler.resize(0);
  // mutationScaler_atPrev.resize(0);
	// logLikelihood.resize(0);
	// logLikelihood_atPrev.resize(0);
  logPrior_trees.resize(0);
}


void Chain::deleteListTrees()
{
  for(unsigned int i=0; i<list_trees.size(); i++)
    {
      list_trees.at(i)->deleteTopo();
    }
  list_trees.resize(0);
}

void Chain::deleteTrees()
{
  for(unsigned int i=0; i< numSubLoci; i++)
    {
      trees_atPrev.at(i)->deleteCoalTree();
    }
  trees_atPrev.resize(0);
  deleteListTrees();
}


void Chain::SetChainID(int c_id)
{
	chainID = c_id;
	return;
}

/**
 * YC 3/3/2014
 * Print out the initial state of a chain.
 */
/*
void Chain::print_chain_initial()
{
	cout << "The number of loci is "<< n_loci <<".\n";
	for(unsigned int i=0; i< n_loci; i++)
	{
		cout << "Locus " << i <<": ";
		trees_gen.at(0).at(i)->print_coaltree();
		cout << "log-likelihood: " << logLikelihood.at(0).at(i) <<"\n";
	}
	cout << "log(joint prior of trees): "<< logPrior_trees.at(0) << "\n";
}
*/

/**
 * \function Chain::PrintChain() Printing the current state of a current chain
 * \return NULL
 */
void Chain::PrintChain()
{
	std::cout << "Chain ID is " << chainID << "\n";
	//FIXME: commending for now since temp is not set std::cout << "Temperature is " << temperature << "\n";
	std::cout << "The number of loci is " << n_loci << "\n";
	for (unsigned int i = 0; i < (unsigned) n_loci; i++) {
		std::cout << "Locus " << i + 1 << ": ";
		try {
			// FIXME
			// The following prints out the initial state, but we want to print out the current state
			std::cout << "log-likelihood: " << logLikelihood.at(0).at(i) << "\n";
		} catch (std::exception &e) {
			std::cout << "Can't access likelihoods to pring - vector index out of bounds\n";
		}
	/* AS: Old code */
	//std::cout << "Posterior is " << posterior << "\n";
	std::cout << "Done!\n";
	}
	return;
}/* End of function Chain::PrintChain() */



void Chain::saveTrees2File()
{
	ofstream treefile;
	treefile.open ("treeSample.txt");
	treefile << "iter\tlocus\tTree\n";
	for(unsigned int i=0; i<n_MCMCgen; i++)
		for(unsigned int j=0; j<n_loci; j++)
		{
			treefile << i <<"\t" <<j <<"\t";
			string tree = list_trees.at(treeIDs.at(i).at(j))->convert2Newick(coalTimes.at(i).at(j), tipIDs.at(i).at(j));
			treefile << tree;
			// trees.at(i).at(j)->saveTree(treefile);
		}
	treefile.close();
}


void Chain::saveOriginalTreeTopo2File()
{
  ofstream treefile;
  treefile.open ("originalTreeTopoSample.txt");
  treefile << "iter\tlocus\tTopo\n";
  for(unsigned int i=0; i<n_MCMCgen; i++)
    for(unsigned int j=0; j<n_loci; j++)
      {
	treefile << i <<"\t" <<j <<"\t";
	string tree = list_trees.at(treeIDs.at(i).at(j))->convert2Newick_originalTopo(tipIDs.at(i).at(j));
	treefile << tree;
	// trees.at(i).at(j)->saveTree(treefile);
      }
  treefile.close();
}

void Chain::saveLogEachLikelihoods()
{
  ofstream treefile;
  treefile.open ("MCMCsample_logEachLikelihoods.txt");
  treefile << "iter\tlocus\tlogLik\n";
  for(unsigned int i=0; i<logLikelihood.size(); i++)
    {
      for(unsigned int j=0; j<logLikelihood.at(i).size(); j++)
	{
	  treefile << i <<"\t" << j <<"\t" << logLikelihood.at(i).at(j) <<"\n";
	}
    }
  treefile.close();
}

void Chain::saveEachLogPrior()
{
  ofstream treefile;
  treefile.open ("MCMCsample_eachLogPrior.txt");
  treefile << "iter\tlocus\tlogPrior\n";
  for(unsigned int i=0; i<logPrior_each.size(); i++)
    {
      for(unsigned int j=0; j<logPrior_each.at(i).size(); j++)
	{
	  treefile << i <<"\t" << j <<"\t" << logPrior_each.at(i).at(j) <<"\n";
	}
    }
  treefile.close();
}

void Chain::saveListTopo2File()
{
  // std::cout <<"In Chain::saveListTopo2File()\n";
  // std::cout << "list_trees.size() = " << list_trees.size() <<"\n";

  ofstream treefile;
  treefile.open ("ListTopo.txt");
  treefile << "treeID\ttopo\n";
  string tree;
  if(Forest ==0||Forest ==2)
    {
      for(unsigned int i=0; i<list_trees.size(); i++)
	{
	  treefile << i <<"\t";
	  tree = list_trees.at(i)->convert2Newick_topo();
	  treefile << tree;
	  treefile << "\n";
	  // trees.at(i).at(j)->saveTree(treefile);
	}
    }
  else if(Forest ==1 ||Forest ==3)
    {
      for(unsigned int i=0; i<list_trees.size(); i++)
	{
	  treefile << i <<"\t";
	  tree = list_trees.at(i)->convert2Newick_topo();
	  treefile << tree;
	  treefile << "\n";
	  // trees.at(i).at(j)->saveTree(treefile);
	}
    }
  treefile.close();
}

void Chain::saveLogPriorTrees()
{
  ofstream treefile;
  treefile.open ("MCMCsample_logJointPriorTrees.txt");
  treefile << "iter\tlogPrior\n";
  for(unsigned int i=0; i<logPrior_trees.size(); i++)
    {
      treefile << i <<"\t" << logPrior_trees.at(i) <<"\n";
    }
  treefile.close();
}


void Chain::saveLogJointLikelihoods()
{
  ofstream treefile;
  treefile.open ("MCMCsample_logJointLikelihoods.txt");
  treefile << "iter\tlogJointLik\n";
  for(unsigned int i=0; i<logJointLikelihood.size(); i++)
    {
      treefile << i <<"\t" << logJointLikelihood.at(i) <<"\n";
    }
  treefile.close();
}

void Chain::saveRNormalDist()
{
  ofstream treefile;
  treefile.open ("MCMCsample_RNormalDist.txt");
  treefile << "iter\trNormalDist\n";
  for(unsigned int i=0; i<rNormalDist.size(); i++)
    {
      treefile << i <<"\t" << rNormalDist.at(i) <<"\n";
    }
  treefile.close();
}


void Chain::saveCoalTimes2File()
{
  ofstream treefile;
  treefile.open ("coalescentTimes.txt");
  treefile.precision(10);
  unsigned int nCT = coalTimes.at(0).at(0).size();
  treefile << "iter\tlocus";
  for(unsigned int i=0; i<nCT; i++)
    {
      treefile << "\tct" << i;
    }
  treefile << "\n";
  
  for(unsigned int i=0; i<coalTimes.size(); i++)
    {
      for(unsigned int j=0; j<coalTimes.at(i).size(); j++)
	{
	  treefile << i <<"\t" <<j;
	  //list<double> coaltime = coalTimes.at(i).at(j);
	  //coaltime.sort();
	  coalTimes.at(i).at(j).sort();
	  for(list<double>::iterator it= coalTimes.at(i).at(j).begin() ;it != coalTimes.at(i).at(j).end(); ++it)
	    {
	      treefile << "\t" << *it;
	    }
	  treefile <<"\n";
	}
    }
  treefile.close();
}



void Chain::SetNumLoci(int numloci)
{
	n_loci = numloci;
	return;
}
/* End of function Chain::SetNumLoci() */

// FIXME YC 3/2/2014
// This function is never called. Do you need it?
void Chain::SetNumIterations(int numiterations)
{
	n_MCMCgen = numiterations;
	return;
}


void Chain::print_coalTimes(unsigned int iter, unsigned int id_locus)
{
	std::list<double>::iterator i;
	for(i = coalTimes.at(iter).at(id_locus).begin(); i != coalTimes.at(iter).at(id_locus).end(); i++)
		std::cout << *i <<" ";
	std::cout <<"\n";
}


void Chain::print_tipIDs(unsigned int iter, unsigned int id_locus)
{
	for(unsigned int i = 0; i < tipIDs.at(iter).at(id_locus).size(); i++)
		std::cout << tipIDs.at(iter).at(id_locus).at(i) <<" ";
	std::cout <<"\n";
}

// YC 9/8/2014
void Chain::print_states_atIter(unsigned int iter)
{
  for(unsigned int index_locus = 0; index_locus < numSubLoci; index_locus++)
    {
      try{
	std::cout << temperature << " ";
	if(lociInParallel ==1)
	  std::cout << logPriorTrees_atPrev_lociParallel.at(index_locus) <<" " << index_locus+locusID_start;
	else
	  std::cout <<logPriorTrees_atPrev << " " << index_locus;
	std::cout << " "<< logLikelihood_atPrev.at(index_locus)
		  << " " << kappa_atPrev.at(index_locus)
		  << " " << mutationScaler_atPrev.at(index_locus) <<" ";
	//	if(list_trees.size()>0)
	//	std::cout << treeIDs.at(iter).at(index_locus);
	std::cout <<" ";
	trees_atPrev.at(index_locus)->print_coaltree();
	//if(list_trees.size()>0)
	//{
	//	print_coalTimes(iter,index_locus);
	//	print_tipIDs(iter,index_locus);
	//}
      } catch (std::exception &e) {
	std::cout << "In Chain::print_states_atIter(). "
	  "Can't access elements of Chain - index out of bounds\n";
      }
    }
  
  return;
}

// YC 9/18/2014
std::string Chain::print2string_states_atCrr()
{
  std::string state;
  for(unsigned int index_locus = 0; index_locus < n_loci; index_locus++)
    {
      std::stringstream convert;
      try
	{
	  convert << temperature << " " <<logPriorTrees_atPrev <<" "<< index_locus << " "<< logLikelihood_atPrev.at(index_locus)
		  << " " << kappa_atPrev.at(index_locus)
		  << " " << mutationScaler_atPrev.at(index_locus) <<" ";
	  convert <<" ";
	  state += convert.str();
	  state += trees_atPrev.at(index_locus)->convert2string_newick_root();
	  state += "\n";
	} 
      catch (std::exception &e) 
	{
	  std::cout << "In Chain::print2string_states_atIter. "
	    "Can't access elements of Chain - index out of bounds\n";
	}
    }
  return state;
}

void Chain::print_treeAccptRate()
{
  for(unsigned int i=0; i<n_loci; i++)
    std::cout << numAccpt_trees(0,i)  <<"/"<<numTries_trees(0,i)<< " ";
  std::cout <<"\n";
}


// Old version - YC 9/8/2014
/*
void Chain::print_states_atIter(unsigned int iter)
{
	try{
		std::cout << "temperature: " << temperature
			<< ", log(prior of trees): "<< logPriorTrees_atPrev <<"\n";
	} catch (std::exception &e) {
		std::cout << "In Chain::print_states_atIter(). "
					"Can't access logPrior_trees - index out of bounds\n";
	}
	if(list_trees.size()>0)
		std::cout << "Locus\tlog-lik\tKappa\t\tmutScaler\tTopoID\tTree\n";
	else
		std::cout << "Locus\tlog-lik\tKappa\t\tmutScaler\tTree\n";

	for(unsigned int index_locus = 0; index_locus < n_loci; index_locus++)
	{
		try{
			std::cout << index_locus << "\t"<< logLikelihood_atPrev.at(index_locus)
				<< "\t" << kappa_atPrev.at(index_locus)
				<< "\t" << mutationScaler_atPrev.at(index_locus) <<"\t";
			//	if(list_trees.size()>0)
			//	std::cout << treeIDs.at(iter).at(index_locus);
			std::cout <<"\t";
			trees_atPrev.at(index_locus)->print_coaltree();
			//if(list_trees.size()>0)
			//{
			//	print_coalTimes(iter,index_locus);
			//	print_tipIDs(iter,index_locus);
			//}
		} catch (std::exception &e) {
				std::cout << "In Chain::print_states_atIter(). "
						"Can't access elements of Chain - index out of bounds\n";
		}
	}
	std::cout <<"\n";

	// Added by YC 8/7/2014
	std::cout <<"\n\nNo.accept/No.tries: ";
	for(unsigned int i=0; i<n_loci; i++)
	  std::cout << numAccpt_trees.at(i)  <<"/"<<numTries_trees.at(i)<< " ";
	std::cout <<"\n\n";
	

	print_listTrees();
	std::cout << "\n";

	return;
}
*/


void Chain::print_savedStates(unsigned int iter)
{
	try{
	  std::cout << "log(prior of trees): "<< logPrior_trees.at(iter) <<"\n";
	} catch (std::exception &e) {
		std::cout << "In Chain::print_states_atIter(). "
					"Can't access logPrior_trees - index out of bounds\n";
	}
	if(list_trees.size()>0)
		std::cout << "Locus\tlog-lik\tKappa\t\tmutScaler\tTopoID\tTree\n";
	else
		std::cout << "Locus\tlog-lik\tKappa\t\tmutScaler\tTree\n";

	for(unsigned int index_locus = 0; index_locus < n_loci; index_locus++)
	{
		try{
		  std::cout << index_locus << "\t"<< logLikelihood.at(iter).at(index_locus)
				<< "\t" << kappa.at(iter).at(index_locus)
				<< "\t" << mutationScaler.at(iter).at(index_locus) <<"\t";
		  if(list_trees.size()>0)				  
		    std::cout << treeIDs.at(iter).at(index_locus);
		  if(list_trees.size()>0)
		    {
		      // list_trees.at(treeIDs.at(iter).at(index_locus))->print_topo();
		      std::string ctree = list_trees.at(treeIDs.at(iter).at(index_locus))->convert2Newick(coalTimes.at(iter).at(index_locus), tipIDs.at(iter).at(index_locus));
		      std::cout << ctree <<"\n";
		    }
		  std::cout <<"\t";
		  if(list_trees.size()>0)
			{
				print_coalTimes(iter,index_locus);
				print_tipIDs(iter,index_locus);
			}
		} catch (std::exception &e) {
				std::cout << "In Chain::print_savedStates(). "
						"Can't access elements of Chain - index out of bounds\n";
		}
	}
	std::cout <<"\n";

	print_listTrees();
	std::cout << "\n";

	return;
}


void Chain::fileprint_states_atIter(unsigned int iter, std::ofstream& fp)
{
	try{
		fp << "temperature: " << temperature
			<< ", log(prior of trees): "<< logPrior_trees.at(iter)<<"\n";
	} catch (std::exception &e) {
		fp << "In Chain::print_states_atIter(). "
					"Can't access logPrior_trees - index out of bounds\n";
	}
	if(list_trees.size()>0)
		fp << "Locus\tlog-lik\tKappa\t\tmutScaler\tTopoID\tTree\n";
	else
		fp << "Locus\tlog-lik\tKappa\t\tmutScaler\tTree\n";

	for(unsigned int index_locus = 0; index_locus < n_loci; index_locus++)
	{
		try{
			fp << index_locus << "\t"<< logLikelihood_atPrev.at(index_locus)
				<< "\t" << kappa_atPrev.at(index_locus)
				<< "\t" << mutationScaler_atPrev.at(index_locus) <<"\t";
			if(list_trees.size()>0)
				std::cout << treeIDs.at(iter).at(index_locus);
			fp <<"\t";
			trees_atPrev.at(index_locus)->print_coaltree();
			if(list_trees.size()>0)
			{
				print_coalTimes(iter,index_locus);
				print_tipIDs(iter,index_locus);
			}
		} catch (std::exception &e) {
				fp << "In Chain::print_states_atIter(). "
						"Can't access elements of Chain - index out of bounds\n";
		}
	}
	fp <<"\n";

	/*
	if(list_trees.size()>0)
	{
		std::cout << "Locus\ttreeID\n";
		for(unsigned int index_locus = 0; index_locus < n_loci; index_locus++)
		{
			std::cout <<index_locus << "\t" << treeIDs.at(iter).at(index_locus) <<"\n";
		}
		std::cout << "\n";
	}*/
	fileprint_listTrees(fp);
	fp << "\n";

	return;
}


/*
void Chain::print_trees_lik()
{
	cout << "Gen\tLocus\tTree\n";
	int ngen =0;
	vector<vector<node*> >::iterator iter_gen;
	for(iter_gen = trees_gen.begin(); iter_gen != trees_gen.end(); iter_gen++)
	{
		int nloci = 0;
		cout << ngen << "\t";
		vector<node*>::iterator iter_lc;
		for(iter_lc=iter_gen->begin(); iter_lc != iter_gen->end(); iter_lc++)
		{
			nloci++;
			cout << nloci << "\t";
			(*iter_lc)->print_coaltree();
		}
		ngen++;
	}
}
*/


double Chain::compute_totalLogLik()
{
  double totalLogLik = 0.0;
  for(unsigned int l=0; l<numSubLoci; l++)
    {
      totalLogLik += logLikelihood_atPrev.at(l);
    }
  return totalLogLik;
}

void Chain::compute_rankTemp()
{
  rankTemp = static_cast<int> (1/temperature-1)/10;
  return;
}


void Chain::compute_totalnum_coalEvents(std::vector<locus> loci)
{
	totalNum_coalEvents = 0;

	for(unsigned int i=0; i <n_loci; i++)
	{
	  totalNum_coalEvents += loci.at(i).get_nGeneCopies()-1;
	}

	return;
}


void Chain::compute_coalTimes_tipIDs(unsigned int iter, unsigned int id_locus, node* tree)
{
  if(tree->size_tree() > 1)
    {
      if(coalTimes.at(iter).size() == 0)
	coalTimes.at(iter).resize(n_loci);
      coalTimes.at(iter).at(id_locus).push_back(tree->age);
    }

  if(tree->isTip == 0)
    {
      if(tree->desc[0]->age < tree->desc[1]->age) // younger child is the first
	{
	  compute_coalTimes_tipIDs(iter,id_locus,tree->desc[0]);
	  compute_coalTimes_tipIDs(iter,id_locus,tree->desc[1]);
	}
      else if(tree->desc[0]->age > tree->desc[1]->age) 
	{
	  compute_coalTimes_tipIDs(iter,id_locus,tree->desc[1]);
	  compute_coalTimes_tipIDs(iter,id_locus,tree->desc[0]);
	}
      else if(tree->desc[0]->isTip == 1 && tree->desc[1]->isTip ==1) 
	{
	  if(tree->desc[0]->tipID < tree->desc[1]->tipID) // smaller tipID is the first
	    {
	      compute_coalTimes_tipIDs(iter,id_locus,tree->desc[0]);
	      compute_coalTimes_tipIDs(iter,id_locus,tree->desc[1]);	      
	    }
	  else if(tree->desc[0]->tipID > tree->desc[1]->tipID) 
	    {
	      compute_coalTimes_tipIDs(iter,id_locus,tree->desc[1]);
	      compute_coalTimes_tipIDs(iter,id_locus,tree->desc[0]);
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
	  compute_coalTimes_tipIDs(iter,id_locus,tree->desc[0]);
	  compute_coalTimes_tipIDs(iter,id_locus,tree->desc[1]);	 
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
      if(tipIDs.at(iter).size() == 0)
	tipIDs.at(iter).resize(n_loci);
      tipIDs.at(iter).at(id_locus).push_back(tree->tipID);
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

/***
 *	YC 4/28/2014
 *		It returns the joint prior distribution, on log scale, of trees with
 *	a proposed tree (newTree) at locus "id_locus".
 *	This prior distribution is used to determine whether the proposed
 *	tree is accepted or not at MCMC generation "id_crr_iter".
 *		For example, let n be the id_crr_iter, k be the id_locus,
 *	t* be the proposed tree, and t(i,n) be a coalescent tree
 *	for ith locus at nth generation.
 *	Then this function returns the following probability on log scale:
 *	   Pr( t(0,n), ..., t(k-1,n), t*, t(k+1,n-1), ..., t(n_loci-1,n-1)).
 *
 *		Here we consider a single population and all the trees
 *	lie on the single population under the coalescent. The prior
 *	distribution we compute here is the marginal distribution of trees
 *	when we consider an Uniform distribution on the population size
 *	and integrate out it:
 *	    P( trees) = \int P(trees | population size) P(population size) d(population size)
 *
 */
// YC 2/23/2015
// See Eq (13) from Hey and Nielsen (2007)
// When 'totalCoalRate' is defined as Eq (12) and Eq (16), then
// the correct expression of Eq (13) is
// 2^(n-1) * f^{-(n-2)} * Gamma(n-2, fn/theta_max) / theta_max, where 
// n is the number of sequences.
double Chain::compute_logPrior_usingTrees_atPrev()
{
	double totalCoalRate_allLoci = 0.0;
	double log_prior = 0.0;
	for(unsigned int i=0; i <n_loci; i++)
	{
		trees_atPrev.at(i)->compute_totalCoalescentRate();
		totalCoalRate_allLoci += trees_atPrev.at(i)->get_totalCoalescentRate();
	
	}

	log_prior = (double) totalNum_coalEvents*log(2.0) 
	  - (double) (totalNum_coalEvents-n_loci) * log(totalCoalRate_allLoci)
	  + log_incompleteUpperGamma (totalNum_coalEvents-n_loci, totalCoalRate_allLoci/popSizeMax)
	  - log(popSizeMax);
	

	return log_prior;
}




/***
 *	YC 3/4/2014
 *		It returns the joint prior distribution, on log scale, of trees with
 *	a proposed tree (newTree) at locus "id_locus".
 *	This prior distribution is used to determine whether the proposed
 *	tree is accepted or not at MCMC generation "id_crr_iter".
 *		For example, let n be the id_crr_iter, k be the id_locus,
 *	t* be the proposed tree, and t(i,n) be a coalescent tree
 *	for ith locus at nth generation.
 *	Then this function returns the following probability on log scale:
 *	   Pr( t(0,n), ..., t(k-1,n), t*, t(k+1,n-1), ..., t(n_loci-1,n-1)).
 *
 *		Here we consider a single population and all the trees
 *	lie on the single population under the coalescent. The prior
 *	distribution we compute here is the marginal distribution of trees
 *	when we consider an Uniform distribution on the population size
 *	and integrate out it:
 *	    P( trees) = \int P(trees | population size) P(population size) d(population size)
 *
 */
double Chain::compute_logPriorTrees(node* newTree, unsigned int id_locus)
{
  double totalCoalRate_allLoci = 0.0;
  double log_prior = 0.0;
  for(unsigned int i=0; i <n_loci; i++)
    {
      // REMOVE
      //std::cout << "id_crr_iter = "<< id_crr_iter << ", i= " << i<< ", id_locus = "<< id_locus <<"\n";
      if(i == id_locus)
	{
	  totalCoalRate_allLoci += newTree->get_totalCoalescentRate();
	  
	  // REMOVE
	  //	std::cout << "At locus " <<i <<", coalRate = "<< newTree->get_totalCoalescentRate() << "\n";
	}
      
      else
	{
	  try{
	    totalCoalRate_allLoci +=
	      trees_atPrev.at(i)->get_totalCoalescentRate();
	  } catch (std::exception &e) {
	    std::cout << "In Chain::compute_logPriorTrees(). "
	      "Can't access elements of trees_gen - index out of bounds\n";
	  }
	  
	  // REMOVE
	  // std::cout << "At locus " <<i <<", coalRate = "<< trees_atPrev.at(i)->get_totalCoalescentRate() << "\n";
	}
      
    }

  // The upper bound of population size used in prior is the same as the upper bound of population sizes in L mode
  log_prior = (double) totalNum_coalEvents*log(2.0) - (double) (totalNum_coalEvents-1) * log(totalCoalRate_allLoci)
    + log_incompleteUpperGamma (totalNum_coalEvents-1, totalCoalRate_allLoci/(popSizeMax))
    - log(popSizeMax);
  
  return log_prior;
}



/***
 * Metropolis-Hastings ratio for updating a tree
 */
double Chain::MHratio_trees(unsigned int id_locus, node* newTree, double newlogLik, double logRatioPriors, double slideDist)
{
	// Remove later
	//std::cout <<"in MHratio_trees()\n";
	//std::cout << "old tree: ";
	//trees_gen.at(id_prevIter).at(id_locus)->print_coaltree();
	//std::cout << "new tree: ";
	//newTree->print_coaltree();

	// Remove later
	//std::cout << "The likelihood of the new tree: " << newlogLik <<"\n";


	double logRatioLik = 0.0;
	try{
		logRatioLik = newlogLik - logLikelihood_atPrev.at(id_locus);
	} catch (std::exception &e) {
		std::cout << "In Chain::MHratio_trees(). "
				"Can't access elements of logLikelihood - index out of bounds\n";
	}

	// FIXME YC 3/5/2014
	// Double-check the proposal ratio/
	// May want to see the original paper, Nielsen (2000, Genetics)
	double logProposalRatio = 0.0;
	try{
		if(newTree->age != trees_atPrev.at(id_locus)->age)
		{
			double old_stdv = std::min(trees_atPrev.at(id_locus)->age/3,20.0);
			double new_stdv = std::min(newTree->age/3,20.0);
			// the log of the proposal ratio of
			// f(slideDist; 0, new_stdv) to f(slideDist; 0, old_stdv),
			// where f(x; mu, stdv) is a Normal density at value x
			//       with mean mu and standard deviation stdv.
			logProposalRatio = -log(old_stdv) - pow(slideDist,2)/(2*pow(old_stdv,2))
			  + log(new_stdv) + pow(slideDist,2)/(2*pow(new_stdv,2));  // Modified by YC 6/9/2014
			//logProposalRatio = log(old_stdv)-log(new_stdv)
			//	+ pow(slideDist,2)*( 1/pow(old_stdv,2) - 1/pow(new_stdv,2) )/2;
		}
	} catch (std::exception &e) {
		std::cout << "In Chain::MHratio_trees(). "
				"Can't access elements of trees_gen - index out of bounds\n";
	}

	// Metropolis-Hastings ratio
	double MHratio = exp( temperature *(logRatioLik + logRatioPriors) + logProposalRatio);
	return MHratio;
}

/***
 * Metropolis-Hastings ratio for updating a tree
 */
double Chain::MHratio_trees_newProposal(unsigned int id_locus, node* newTree, double newlogLik, double logRatioPriors, double slideDist)
{

  double logRatioLik = 0.0;
  try{
    logRatioLik = newlogLik - logLikelihood_atPrev.at(id_locus);
  } catch (std::exception &e) {
    std::cout << "In Chain::MHratio_trees(). "
      "Can't access elements of logLikelihood - index out of bounds\n";
  }
  
  // FIXME YC 3/5/2014
  // Double-check the proposal ratio/
  // May want to see the original paper, Nielsen (2000, Genetics)
  double logProposalRatio = 0.0;
	/*
	try{
		if(newTree->age != trees_atPrev.at(id_locus)->age)
		{
			double old_stdv = std::min(trees_atPrev.at(id_locus)->age/3,20.0);
			double new_stdv = std::min(newTree->age/3,20.0);
			// the log of the proposal ratio of
			// f(slideDist; 0, new_stdv) to f(slideDist; 0, old_stdv),
			// where f(x; mu, stdv) is a Normal density at value x
			//       with mean mu and standard deviation stdv.
			logProposalRatio = -log(old_stdv) - pow(slideDist,2)/(2*pow(old_stdv,2))
			  + log(new_stdv) + pow(slideDist,2)/(2*pow(new_stdv,2));  // Modified by YC 6/9/2014
			//logProposalRatio = log(old_stdv)-log(new_stdv)
			//	+ pow(slideDist,2)*( 1/pow(old_stdv,2) - 1/pow(new_stdv,2) )/2;
		}
	} catch (std::exception &e) {
		std::cout << "In Chain::MHratio_trees(). "
				"Can't access elements of trees_gen - index out of bounds\n";
	}
	*/

	// Metropolis-Hastings ratio
	//double MHratio = exp( temperature *(logRatioLik + logRatioPriors) + logProposalRatio);
  double MHratio = exp( temperature *(logRatioLik + logRatioPriors));

  // REMOVE
  /*
  std::cout << "temperature = " << temperature << " and logRatioLik = " << logRatioLik 
	    << " logRatioPriors = " << logRatioPriors << " and MHratio = " << MHratio <<"\n";
  */
  return MHratio;
}


/***
 * Metropolis-Hastings ratio for updating a tree
 */
double Chain::MHratio_trees_usePriorOnly(unsigned int id_locus, node* newTree, double newlogLik, double logRatioPriors, double slideDist)
{
	// Remove later
	//std::cout <<"in MHratio_trees()\n";
	//std::cout << "old tree: ";
	//trees_gen.at(id_prevIter).at(id_locus)->print_coaltree();
	//std::cout << "new tree: ";
	//newTree->print_coaltree();

	// Remove later
	//std::cout << "The likelihood of the new tree: " << newlogLik <<"\n";


	double logRatioLik = 0.0;
	/*
	try{
		logRatioLik = newlogLik - logLikelihood_atPrev.at(id_locus);
	} catch (std::exception &e) {
		std::cout << "In Chain::MHratio_trees(). "
				"Can't access elements of logLikelihood - index out of bounds\n";
	}
	*/
	// FIXME YC 3/5/2014
	// Double-check the proposal ratio/
	// May want to see the original paper, Nielsen (2000, Genetics)
	double logProposalRatio = 0.0;
	try{
		if(newTree->age != trees_atPrev.at(id_locus)->age)
		{
			double old_stdv = std::min(trees_atPrev.at(id_locus)->age/3, 20.0);
			double new_stdv = std::min(newTree->age/3, 20.0);
			// the log of the proposal ratio of
			// f(slideDist; 0, new_stdv) to f(slideDist; 0, old_stdv),
			// where f(x; mu, stdv) is a Normal density at value x
			//       with mean mu and standard deviation stdv.
			// logProposalRatio = log q(old) - log q(new)
			logProposalRatio = -log(old_stdv) - pow(slideDist,2)/(2*pow(old_stdv,2))
			  + log(new_stdv) + pow(slideDist,2)/(2*pow(new_stdv,2));  // Modified by YC 6/9/2014
			//logProposalRatio = log(old_stdv)-log(new_stdv)
			//	+ pow(slideDist,2)*( 1/pow(old_stdv,2) - 1/pow(new_stdv,2) )/2;
		}
	} catch (std::exception &e) {
		std::cout << "In Chain::MHratio_trees(). "
				"Can't access elements of trees_gen - index out of bounds\n";
	}

	// Metropolis-Hastings ratio
	double MHratio = exp( temperature *(logRatioLik + logRatioPriors) + logProposalRatio);
	return MHratio;
}

/***
 * Metropolis-Hastings ratio for updating a tree
 */
double Chain::MHratio_trees_usePriorOnly_newProposal(unsigned int id_locus, node* newTree, double newlogLik, double logRatioPriors, double slideDist)
{
	// Remove later
	//std::cout <<"in MHratio_trees()\n";
	//std::cout << "old tree: ";
	//trees_gen.at(id_prevIter).at(id_locus)->print_coaltree();
	//std::cout << "new tree: ";
	//newTree->print_coaltree();

	// Remove later
	//std::cout << "The likelihood of the new tree: " << newlogLik <<"\n";


	double logRatioLik = 0.0;
	/*
	try{
		logRatioLik = newlogLik - logLikelihood_atPrev.at(id_locus);
	} catch (std::exception &e) {
		std::cout << "In Chain::MHratio_trees(). "
				"Can't access elements of logLikelihood - index out of bounds\n";
	}
	*/
	// FIXME YC 3/5/2014
	// Double-check the proposal ratio/
	// May want to see the original paper, Nielsen (2000, Genetics)
	double logProposalRatio = 0.0;
	/*
	try{
		if(newTree->age != trees_atPrev.at(id_locus)->age)
		{
			double old_stdv = std::min(trees_atPrev.at(id_locus)->age/3, 20.0);
			double new_stdv = std::min(newTree->age/3, 20.0);
			// the log of the proposal ratio of
			// f(slideDist; 0, new_stdv) to f(slideDist; 0, old_stdv),
			// where f(x; mu, stdv) is a Normal density at value x
			//       with mean mu and standard deviation stdv.
			// logProposalRatio = log q(old) - log q(new)
			logProposalRatio = -log(old_stdv) - pow(slideDist,2)/(2*pow(old_stdv,2))
			  + log(new_stdv) + pow(slideDist,2)/(2*pow(new_stdv,2));  // Modified by YC 6/9/2014
			//logProposalRatio = log(old_stdv)-log(new_stdv)
			//	+ pow(slideDist,2)*( 1/pow(old_stdv,2) - 1/pow(new_stdv,2) )/2;
		}
	} catch (std::exception &e) {
		std::cout << "In Chain::MHratio_trees(). "
				"Can't access elements of trees_gen - index out of bounds\n";
	}
	*/

	// Metropolis-Hastings ratio
	//double MHratio = exp( temperature *(logRatioLik + logRatioPriors) + logProposalRatio);
	double MHratio = exp( temperature *(logRatioLik + logRatioPriors));
	return MHratio;
}

// the distance for sliding along a tree is from a Normal distribution with a fixed standard deviation
void Chain::UpdateTrees_newProposal(unsigned int id_crrIter, vector<locus> lc, unsigned int save)
{  

  // REMOVE
  std::cout << "\n Chain::UpdateTrees_newProposal()\n";
  for(unsigned int i=0; i < n_loci; i++)
    {
      std::cout << "locus = " << i <<"\n";
      
      // Added by YC 8/7/2014
      //if(save==1)
      numTries_trees(0,i)++;
      
      node* crrTree = trees_atPrev.at(i);
      node* newTree = new node;
      newTree->deepCopy(crrTree);
      
      // get the distance to move
      // double stdv = 1; // std::min(newTree->age/3,20.0); -- YC 9/25/2014
      double stdv = std::min(splittingTimeMax/3,20.0);
      double slidedist = rNormalDistribution(0.0, stdv);
      
      newTree = newTree->propose_coaltree(slidedist);
      
      
      // Computing the likelihood of the proposed tree, newTree
      unsigned int likelihoodModel = lc.at(i).getLikelihoodModel();
      double newLogLik = 0.0;
      try 
	{
	  if(likelihoodModel == 0) // IS
	    {
	      newLogLik = newTree->logLikelihood_IS(lc.at(i), mutationScaler_atPrev.at(i),0);
	     
	    }
	  else if(likelihoodModel == 2) // HKY
	    newLogLik = newTree->logLikelihood_HKY(lc.at(i), mutationScaler_atPrev.at(i), kappa_atPrev.at(i));
	} 
      catch (std::exception &e) 
	{
	  std::cout << "Can't access lc to calculate likelihood - vector index out of bounds\n";
	}
      
      unsigned int rejectThisTree = 0;
      if(likelihoodModel == 0)
	if(newTree->get_rejectThisTree_IS() == 1)
	 {
	   rejectThisTree = 1;
	   newTree->set_rejectThisTree_IS(0);	 
	 }

       if(rejectThisTree == 0)
	 {
	   // Computing the joint prior distribution of the trees including "newTree"
	   newTree->compute_totalCoalescentRate();
	   double newlogPrior = compute_logPriorTrees(newTree,i);      
	   double logRatioPrior = newlogPrior - logPriorTrees_atPrev;
	   double MHratio = MHratio_trees_newProposal(i, newTree, newLogLik,logRatioPrior, slidedist);
	   
	   //REMOVE
	   // std::cout << "proposed tree is ";
	   //newTree->print_coaltree();
	   // std::cout << "newLogLik = " << newLogLik <<" ";
	   
	   
	   if(runiform() < min(1.0, MHratio))
	     {	  
	       // accept
	       trees_atPrev.at(i)->deleteCoalTree();
	       trees_atPrev.at(i) = new node;
	       trees_atPrev.at(i)->deepCopy(newTree);
	       logLikelihood_atPrev.at(i)=newLogLik;
	       logPriorTrees_atPrev = newlogPrior;
	       
	       // Added by YC 8/7/2014
	       numAccpt_trees(0,i)++;
	       
	     }	   
	   else  // Reject
	     {
	       
	     } 
	 }
       newTree->deleteCoalTree();

       
       
       
    }
  
  
  return;
}


void Chain::collectAllUpdates_Lmode( unsigned int savingID)
{
  treeIDs.at(savingID).resize(numSubLoci);
  coalTimes.at(savingID).resize(numSubLoci);
  // tipIDs.at(savingID).resize(numSubLoci);

  if(lociInParallel ==0) 
    logPrior_trees.at(savingID) = logPriorTrees_atPrev; // which is 0.
 
  for(unsigned int i=0; i < numSubLoci; i++)
    {
      compute_coalTimes_tipIDs(savingID,i,trees_atPrev.at(i));
      unsigned int list_size = list_trees.size();
      if(list_size == 0)
	{
	  nodeSimple *topo = new nodeSimple;
	  std::list<double> ct = coalTimes.at(savingID).at(i);
	  ct.sort();
	  topo->convert(trees_atPrev.at(i), ct, trees_atPrev.at(i)->size_tree()); // Updated by YC 5/8/2014
	  topo->computeSizes();
	  list_trees.push_back(topo);
	  treeIDs.at(savingID).at(i) = 0;	  
	}
      else
	{
	  unsigned int found_sampleTopo = 0;
	  unsigned int count = 0;
	  while(found_sampleTopo == 0 && count < list_size)
	    {
	      found_sampleTopo = list_trees.at(count)->sameTopo(trees_atPrev.at(i));
	      count++;
	    }

	  if(found_sampleTopo == 0)
	    {
	      nodeSimple *topo = new nodeSimple;
	      std::list<double> ct = coalTimes.at(savingID).at(i);
	      ct.sort();
	      topo->convert(trees_atPrev.at(i),ct, trees_atPrev.at(i)->size_tree());// updated by YC 5/8/2014
	      topo->computeSizes();
	      list_trees.push_back(topo);
	      treeIDs.at(savingID).at(i) = list_trees.size()-1;
	    }
	  else if(found_sampleTopo == 1)
	    {
	      treeIDs.at(savingID).at(i) = count -1;
	    }
	  else
	    {
	      std::cout << "\n*** Error in Chain::collectAllUpdates ***\n";
	      std::cout << "found_sampleTopo = " << found_sampleTopo << "\n";
	      std::cout << "trees_atPrev.at(i) = ";
	      trees_atPrev.at(i)->print_coaltree();
	      nodeSimple *topo = new nodeSimple;
	      std::list<double> ct = coalTimes.at(savingID).at(i);
	      ct.sort();
	      topo->convert(trees_atPrev.at(i),ct, trees_atPrev.at(i)->size_tree());
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
    }
  return;
}

void Chain::collectAllUpdates(unsigned int savingID, unsigned int Lmode)
{
  treeIDs.at(savingID).resize(n_loci);
  coalTimes.at(savingID).resize(n_loci);
  tipIDs.at(savingID).resize(n_loci);
  logPrior_each.at(savingID).resize(n_loci);

  if(Lmode == 0) // In M mode
    {
      if(lociInParallel ==1)
	{
	  kappa.at(savingID).resize(n_loci); 
	  if(n_loci >= 2)
	    mutationScaler.at(savingID).resize(n_loci);
	  logLikelihood.at(savingID).resize(n_loci);
	  logPrior_each.at(savingID).resize(n_loci);

	  for(unsigned int i=0; i<numSubLoci; i++)
	    {
	      kappa.at(savingID).at(i) = kappa_atPrev.at(i);
	      if(n_loci >= 2)
		mutationScaler.at(savingID).at(i) = mutationScaler_atPrev.at(i);
	      
	      logLikelihood.at(savingID).at(i) =logLikelihood_atPrev.at(i);
	      logPrior_each.at(savingID).at(i) = logPriorTrees_atPrev_lociParallel.at(i);
	    }
	}
      else
	{
	  kappa.at(savingID) = kappa_atPrev;
	  if(n_loci >= 2)
	    mutationScaler.at(savingID) = mutationScaler_atPrev;
	  
	  logLikelihood.at(savingID) =logLikelihood_atPrev;
	}
    }

  if(lociInParallel ==0) 
    logPrior_trees.at(savingID) = logPriorTrees_atPrev; // which is 0.

  unsigned int nl = n_loci;
  if(Lmode ==0) // In M mode
    {
      if(lociInParallel ==1)
	nl = numSubLoci;
    }

  for(unsigned int i=0; i < nl; i++)
    {
      compute_coalTimes_tipIDs(savingID,i,trees_atPrev.at(i));
      unsigned int list_size = list_trees.size();
      if(list_size == 0)
	{
	  nodeSimple *topo = new nodeSimple;
	  std::list<double> ct = coalTimes.at(savingID).at(i);
	  ct.sort();
	  topo->convert(trees_atPrev.at(i), ct, trees_atPrev.at(i)->size_tree()); // Updated by YC 5/8/2014
	  topo->computeSizes();
	  list_trees.push_back(topo);

	  // REMOVE
	  //trees_atPrev.at(i)->print_coaltree();
	  //topo->print_topo();
	}
      else
	{
	  unsigned int found_sampleTopo = 0;
	  unsigned int count = 0;
	  while(found_sampleTopo == 0 && count < list_size)
	    {
	      found_sampleTopo = list_trees.at(count)->sameTopo(trees_atPrev.at(i));
	      count++;
	    }
	  if(found_sampleTopo == 0)
	    {
	      nodeSimple *topo = new nodeSimple;
	      std::list<double> ct = coalTimes.at(savingID).at(i);
	      ct.sort();
	      topo->convert(trees_atPrev.at(i),ct, trees_atPrev.at(i)->size_tree());// updated by YC 5/8/2014
	      topo->computeSizes();
	      list_trees.push_back(topo);
	      treeIDs.at(savingID).at(i) = list_trees.size()-1;
	    }
	  else if(found_sampleTopo == 1)
	    {
	      treeIDs.at(savingID).at(i) = count -1;
	    }
	  else
	    {
	      std::cout << "\n*** Error in Chain::collectAllUpdates ***\n";
	      std::cout << "found_sampleTopo = " << found_sampleTopo << "\n";
	      std::cout << "trees_atPrev.at(i) = ";
	      trees_atPrev.at(i)->print_coaltree();
	      nodeSimple *topo = new nodeSimple;
	      std::list<double> ct = coalTimes.at(savingID).at(i);
	      ct.sort();
	      topo->convert(trees_atPrev.at(i),ct, trees_atPrev.at(i)->size_tree());
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
    }
  return;
}



// The distance to move along a tree is from a Normal distribution with a standard deviation depending the previous tree height.
void Chain::UpdateTrees(unsigned int id_crrIter, vector<locus> lc, unsigned int save)
{
  //  std::cout << "In UpdateTrees()\n";
	if(save==1)
	{
		treeIDs.at(id_crrIter).resize(n_loci);
		logLikelihood.at(id_crrIter).resize(n_loci);
		coalTimes.at(id_crrIter).resize(n_loci);
		tipIDs.at(id_crrIter).resize(n_loci);
	}

	for(unsigned int i=0; i < n_loci; i++)
	{
		// REMOVE
		//std::cout << "locus " << i << ":\n";

		node* crrTree = trees_atPrev.at(i);
		node* newTree = new node;
		newTree->deepCopy(crrTree);

			// get the distance to move
		double stdv = std::min(newTree->age/3,20.0);
		double slidedist = rNormalDistribution(0.0, stdv);

		newTree = newTree->propose_coaltree(slidedist);

		// REMOVE
		//std::cout << "locus " << i<<":\n";
		//std::cout << "crrTree is ";
		//crrTree->print_coaltree();
		//std::cout << "new tree is ";
		//tr->print_coaltree();


		// Computing the likelihood of the proposed tree, newTree
		double newLogLik = 0.0;
		try {
			newLogLik = newTree->logLikelihood_HKY(lc.at(i), mutationScaler_atPrev.at(i), kappa_atPrev.at(i));
		} catch (std::exception &e) {
			std::cout << "Can't access lc or kappa to calculate likelihood - vector index out of bounds\n";
		}

		// Computing the joint prior distribution of the trees including "newTree"
		newTree->compute_totalCoalescentRate();
		double newlogPrior = compute_logPriorTrees(newTree,i);
		//double oldLogPrior = logPrior_trees.at(id_crrIter);//compute_logPriorTrees(crrTree,id_crrIter,i);

		double logRatioPrior = newlogPrior - logPriorTrees_atPrev;
		double MHratio = MHratio_trees(i, newTree, newLogLik,logRatioPrior, slidedist );

		// REMOVE
		//std::cout << "MHratio = "<< MHratio <<"\n";

		if(runiform() < min(1.0, MHratio) )
		{
			//std::cout << "accept\n";
			//newTree->print_coaltree();

			// accept
			trees_atPrev.at(i)->deleteCoalTree();
			trees_atPrev.at(i) = new node;
			trees_atPrev.at(i)->deepCopy(newTree);
			logLikelihood_atPrev.at(i)=newLogLik;
			logPriorTrees_atPrev = newlogPrior;


		}
		else  // Reject
		{
			// REMOVE
			//std::cout << "reject\n";
			//trees_atPrev.at(i)->print_coaltree();

		}
		newTree->deleteCoalTree();

		if(save==1) // sampling
		{
			compute_coalTimes_tipIDs(id_crrIter,i,trees_atPrev.at(i));
			logLikelihood.at(id_crrIter).at(i) =logLikelihood_atPrev.at(i);

			unsigned int list_size = list_trees.size();
			if(list_size == 0)
			{
				nodeSimple *topo = new nodeSimple;
				std::list<double> ct = coalTimes.at(id_crrIter).at(i);
				ct.sort();
				topo->convert(trees_atPrev.at(i),ct, trees_atPrev.at(i)->size_tree()); // Updated by YC 5/8/2014
				topo->computeSizes();
				list_trees.push_back(topo);
			}
			else
			{
				unsigned int found_sampleTopo = 0;
				unsigned int count = 0;
				while(found_sampleTopo == 0 && count < list_size)
				{
					found_sampleTopo = list_trees.at(count)->sameTopo(trees_atPrev.at(i));
					count++;
				}
				if(found_sampleTopo == 0)
				{
					nodeSimple *topo = new nodeSimple;
					std::list<double> ct = coalTimes.at(id_crrIter).at(i);
					ct.sort();
					topo->convert(trees_atPrev.at(i),ct, trees_atPrev.at(i)->size_tree());// updated by YC 5/8/2014
					topo->computeSizes();
					list_trees.push_back(topo);
					treeIDs.at(id_crrIter).at(i) = list_trees.size()-1;
				}
				else
				{
					treeIDs.at(id_crrIter).at(i) = count -1;
				}
			}

		}

	}

	if(save==1)
		logPrior_trees.at(id_crrIter) = logPriorTrees_atPrev;

	return;
}




void Chain::UpdateTrees_usePriorsOnly(unsigned int id_crrIter, vector<locus> lc, unsigned int save)
{
  //  std::cout << "In UpdateTrees()\n";
	if(save==1)
	{
		treeIDs.at(id_crrIter).resize(n_loci);
		logLikelihood.at(id_crrIter).resize(n_loci);
		coalTimes.at(id_crrIter).resize(n_loci);
		tipIDs.at(id_crrIter).resize(n_loci);
	}

	for(unsigned int i=0; i < n_loci; i++)
	{
		// REMOVE
		//std::cout << "locus " << i << ":\n";

		node* crrTree = trees_atPrev.at(i);
		node* newTree = new node;
		newTree->deepCopy(crrTree);

			// get the distance to move
		double stdv = std::min(newTree->age/3,20.0);
		double slidedist = rNormalDistribution(0.0, stdv);

		// std::cout << "stdv = " << stdv  <<", slidedist = " << slidedist <<"\n";
 

		newTree = newTree->propose_coaltree(slidedist);

		// REMOVE
		//std::cout << "locus " << i<<":\n";
		//std::cout << "crrTree is ";
		//crrTree->print_coaltree();
		//std::cout << "new tree is ";
		//tr->print_coaltree();


		// Computing the likelihood of the proposed tree, newTree
		double newLogLik = 0.0;
		/*
		try {
		  newLogLik = 1.0;// newTree->logLikelihood_HKY(lc.at(i), mutationScaler_atPrev.at(i), kappa_atPrev.at(i));
		} catch (std::exception &e) {
			std::cout << "Can't access lc or kappa to calculate likelihood - vector index out of bounds\n";
		}
		*/
		// Computing the joint prior distribution of the trees including "newTree"
		newTree->compute_totalCoalescentRate();
		double newlogPrior = compute_logPriorTrees(newTree,i);
		//double oldLogPrior = logPrior_trees.at(id_crrIter);//compute_logPriorTrees(crrTree,id_crrIter,i);

		double logRatioPrior = newlogPrior - logPriorTrees_atPrev;
		double MHratio = MHratio_trees_usePriorOnly(i, newTree, newLogLik,logRatioPrior, slidedist );

		// REMOVE
		//std::cout << "logRatioPrior = " << logRatioPrior << ", MHratio = "<< MHratio <<"\n";

		if(runiform() < min(1.0, MHratio) )
		{
			//std::cout << "accept\n";
			//newTree->print_coaltree();

			// accept
			trees_atPrev.at(i)->deleteCoalTree();
			trees_atPrev.at(i) = new node;
			trees_atPrev.at(i)->deepCopy(newTree);
			logLikelihood_atPrev.at(i)= 0.0;// newLogLik;
			logPriorTrees_atPrev = newlogPrior;


		}
		else  // Reject
		{
			// REMOVE
			//std::cout << "reject\n";
			//trees_atPrev.at(i)->print_coaltree();

		}
		newTree->deleteCoalTree();

		if(save==1) // sampling
		{
			compute_coalTimes_tipIDs(id_crrIter,i,trees_atPrev.at(i));
			logLikelihood.at(id_crrIter).at(i) =0.0;//logLikelihood_atPrev.at(i);

			unsigned int list_size = list_trees.size();
			if(list_size == 0)
			{
				nodeSimple *topo = new nodeSimple;
				std::list<double> ct = coalTimes.at(id_crrIter).at(i);
				ct.sort();
				topo->convert(trees_atPrev.at(i),ct, trees_atPrev.at(i)->size_tree()); // Updated by YC 5/8/2014
				topo->computeSizes();
				list_trees.push_back(topo);
			}
			else
			{
				unsigned int found_sampleTopo = 0;
				unsigned int count = 0;
				while(found_sampleTopo == 0 && count < list_size)
				{
					found_sampleTopo = list_trees.at(count)->sameTopo(trees_atPrev.at(i));
					count++;
				}
				if(found_sampleTopo == 0)
				{
					nodeSimple *topo = new nodeSimple;
					std::list<double> ct = coalTimes.at(id_crrIter).at(i);
					ct.sort();
					topo->convert(trees_atPrev.at(i),ct, trees_atPrev.at(i)->size_tree());// updated by YC 5/8/2014
					topo->computeSizes();
					list_trees.push_back(topo);
					treeIDs.at(id_crrIter).at(i) = list_trees.size()-1;
				}
				else
				{
					treeIDs.at(id_crrIter).at(i) = count -1;
				}
			}

		}

	}

	if(save==1)
		logPrior_trees.at(id_crrIter) = logPriorTrees_atPrev;

	return;
}




void Chain::UpdateTrees_usePriorsOnly_newProposal(unsigned int id_crrIter, vector<locus> lc, unsigned int save)
{
  for(unsigned int i=0; i < n_loci; i++)
    {
      // std::cout << "i = " << i <<"\n";
      numTries_trees(0,i)++;

      node* crrTree = trees_atPrev.at(i);
      node* newTree = new node;
      newTree->deepCopy(crrTree);
      
      // get the distance to move
      // double stdv = 1; // std::min(newTree->age/3,20.0); -- YC 9/25/2014
      double stdv = std::min(splittingTimeMax/3,20.0);
      double slidedist = rNormalDistribution(0.0, stdv);
      

      newTree = newTree->propose_coaltree(slidedist);
     

      // Computing the likelihood of the proposed tree, newTree
      double newLogLik = 0.0;
      /*
	try {
	newLogLik = 1.0;// newTree->logLikelihood_HKY(lc.at(i), mutationScaler_atPrev.at(i), kappa_atPrev.at(i));
	} catch (std::exception &e) {
	std::cout << "Can't access lc or kappa to calculate likelihood - vector index out of bounds\n";
	}
      */

      // Computing the joint prior distribution of the trees including "newTree"
      newTree->compute_totalCoalescentRate();
      double newlogPrior = compute_logPriorTrees(newTree,i);
      double logRatioPrior = newlogPrior - logPriorTrees_atPrev;
      double MHratio = MHratio_trees_usePriorOnly_newProposal(i, newTree, newLogLik,logRatioPrior, slidedist );
      
       
      //REMOVE
      /*
      std::cout << "stdv = " << stdv <<" slidedist = "<< slidedist 
		<< " current tree's TMRCA = " << crrTree->age
		<< " Proposed tree's TMRCA = " << newTree->age 
		<<" ";
      */

      if(runiform() < min(1.0, MHratio) )
	{
	  // std::cout << "Accepted\n";
	  // accept
	  trees_atPrev.at(i)->deleteCoalTree();
	  trees_atPrev.at(i) = new node;
	  trees_atPrev.at(i)->deepCopy(newTree);
	  logLikelihood_atPrev.at(i)= 0.0;// newLogLik;
	  logPriorTrees_atPrev = newlogPrior;

	  // Added by YC 8/7/2014
	  numAccpt_trees(0,i)++;

	}
      else
	{	  
	  // std::cout << "Rejected\n";
	}
      newTree->deleteCoalTree();
      
    }
  
  return;
}



/// YC 9/11/2014
void Chain::UpdateTrees_MPI_usePriorsOnly_newProposal(unsigned int id_crrIter, vector<locus> lc, unsigned int save, int currentid, int chainid)
{
  MPI::Status status;
  /*
  double temp = GetTemperature();
  if(save==1 &&  chainid == 0 &&  currentid == 0)
  {
  treeIDs.at(id_crrIter).resize(n_loci);
  }
  logLikelihood.at(id_crrIter).resize(n_loci);
  coalTimes.at(id_crrIter).resize(n_loci);
  tipIDs.at(id_crrIter).resize(n_loci);
  */
	

  for(unsigned int i=0; i < n_loci; i++)
    {
      
      numTries_trees(0,i)++;
      
      node* crrTree = trees_atPrev.at(i);
      node* newTree = new node;
      newTree->deepCopy(crrTree);

      // get the distance to move
      //	double stdv =1;// std::min(newTree->age/3,20.0);
      // double stdv = std::min(splittingTimeMax/6,20.0);
      double stdv = std::min(splittingTimeMax/3,20.0);
      double slidedist = rNormalDistribution(0.0, stdv);
      

      newTree = newTree->propose_coaltree(slidedist);

      // REMOVE
      /*
      std::cout << "id_crrIter = " << id_crrIter 
		<< " slidedist = " << slidedist
		<< " TMRCA(current tree) = " << crrTree->age
		<<" TMRCA( propsed tree) = " << newTree->age <<"\n";
      */

      // Computing the likelihood of the proposed tree, newTree
      double newLogLik = 0.0;
      /*
	try {
	newLogLik = newTree->logLikelihood_HKY(lc.at(i), mutationScaler_atPrev.at(i), kappa_atPrev.at(i));
	} catch (std::exception &e) {
	std::cout << "Can't access lc or kappa to calculate likelihood - vector index out of bounds\n";
	}
      */
      
      // Computing the joint prior distribution of the trees including "newTree"
      newTree->compute_totalCoalescentRate();
      double newlogPrior = compute_logPriorTrees(newTree,i);
      //double oldLogPrior = logPrior_trees.at(id_crrIter);//compute_logPriorTrees(crrTree,id_crrIter,i);      
      double logRatioPrior = newlogPrior - logPriorTrees_atPrev;
      // double MHratio = MHratio_trees_usePriorOnly(i, newTree, newLogLik,logRatioPrior, slidedist);
      double MHratio = MHratio_trees_usePriorOnly_newProposal(i, newTree, newLogLik,logRatioPrior, slidedist);

      if(runiform() < min(1.0, MHratio) )
	{	  
	  // accept
	  trees_atPrev.at(i)->deleteCoalTree();
	  trees_atPrev.at(i) = new node;
	  trees_atPrev.at(i)->deepCopy(newTree);
	  logLikelihood_atPrev.at(i)=0.0;//newLogLik;
	  logPriorTrees_atPrev = newlogPrior;
	  
	  // Added by YC 9/11/2014
	  numAccpt_trees(0,i)++;
	  
	}
      else  // Reject
	{
	  // REMOVE
	  // std::cout << "reject\n";
	  //trees_atPrev.at(i)->print_coaltree();	  
	}
      newTree->deleteCoalTree();
      
    }
  
  return;
}


/// YC 6/6/2014 -- old tree proposal
int Chain::UpdateTrees_MPI_usePriorsOnly(unsigned int id_crrIter, vector<locus> lc, unsigned int save, int currentid, int chainid)
{
	MPI::Status status;
	double temp = GetTemperature();
	if(save==1 &&  chainid == 0 &&  currentid == 0)
	{
		treeIDs.at(id_crrIter).resize(n_loci);
	}
		logLikelihood.at(id_crrIter).resize(n_loci);
	

	coalTimes.at(id_crrIter).resize(n_loci);
	tipIDs.at(id_crrIter).resize(n_loci);
	

	for(unsigned int i=0; i < n_loci; i++)
	{
	  numTries_trees(0,i)++;

		node* crrTree = trees_atPrev.at(i);
		node* newTree = new node;
		newTree->deepCopy(crrTree);

			// get the distance to move
		// double stdv = std::min(newTree->age/3,20.0);
		double stdv = std::min(splittingTimeMax/3,20.0);
		double slidedist = rNormalDistribution(0.0, stdv);

		// std::cout << "stdv = " << stdv  <<", slidedist = " << slidedist <<"\n";

		newTree = newTree->propose_coaltree(slidedist);

		// REMOVE
		//std::cout << "locus " << i<<":\n";
		//std::cout << "crrTree is ";
		//crrTree->print_coaltree();
		//std::cout << "new tree is ";
		//tr->print_coaltree();


		// Computing the likelihood of the proposed tree, newTree
		double newLogLik = 0.0;
		/*
		try {
			newLogLik = newTree->logLikelihood_HKY(lc.at(i), mutationScaler_atPrev.at(i), kappa_atPrev.at(i));
		} catch (std::exception &e) {
			std::cout << "Can't access lc or kappa to calculate likelihood - vector index out of bounds\n";
		}
		*/

		// Computing the joint prior distribution of the trees including "newTree"
		newTree->compute_totalCoalescentRate();
		double newlogPrior = compute_logPriorTrees(newTree,i);
		//double oldLogPrior = logPrior_trees.at(id_crrIter);//compute_logPriorTrees(crrTree,id_crrIter,i);

		double logRatioPrior = newlogPrior - logPriorTrees_atPrev;
		double MHratio = MHratio_trees_usePriorOnly(i, newTree, newLogLik,logRatioPrior, slidedist);

		// REMOVE
		//std::cout << "MHratio = "<< MHratio <<"\n";

		if(runiform() < min(1.0, MHratio) )
		{
			//std::cout << "accept\n";
			//newTree->print_coaltree();

			// accept
			trees_atPrev.at(i)->deleteCoalTree();
			trees_atPrev.at(i) = new node;
			trees_atPrev.at(i)->deepCopy(newTree);
			logLikelihood_atPrev.at(i)=0.0;//newLogLik;
			logPriorTrees_atPrev = newlogPrior;

			// Added by YC 9/11/2014
			numAccpt_trees(0,i)++;

		}
		else  // Reject
		{
			// REMOVE
			//std::cout << "reject\n";
			//trees_atPrev.at(i)->print_coaltree();

		}
		newTree->deleteCoalTree();

	
		double logl = 0.0;
		double logprior = 0.0;
		unsigned int ctsize = 0;
		unsigned int list_size = 0;
		int coldprocess = 0;
		int sendingsignal = 0;
		//AS: this is the case when the head node, chain 0 is still cold
		if(save==1 && currentid == 0  && chainid == 0 && temp == 1.0) // sampling
		{
			compute_coalTimes_tipIDs(id_crrIter,i,trees_atPrev.at(i));
			logLikelihood.at(id_crrIter).at(i) =0.0;//logLikelihood_atPrev.at(i);

			unsigned int list_size = list_trees.size();
			if(list_size == 0)
			{
				nodeSimple *topo = new nodeSimple;
				std::list<double> ct = coalTimes.at(id_crrIter).at(i);
				ct.sort();
				topo->convert(trees_atPrev.at(i),ct, trees_atPrev.at(i)->size_tree()); // Updated by YC 5/8/2014
				topo->computeSizes();
				list_trees.push_back(topo);
			}
			else
			{
				unsigned int found_sampleTopo = 0;
				unsigned int count = 0;
				while(found_sampleTopo == 0 && count < list_size)
				{
					found_sampleTopo = list_trees.at(count)->sameTopo(trees_atPrev.at(i));
					count++;
				}
				if(found_sampleTopo == 0)
				{
					nodeSimple *topo = new nodeSimple;
					std::list<double> ct = coalTimes.at(id_crrIter).at(i);
					ct.sort();
					topo->convert(trees_atPrev.at(i),ct, trees_atPrev.at(i)->size_tree());// updated by YC 5/8/2014
					topo->computeSizes();
					list_trees.push_back(topo);
					treeIDs.at(id_crrIter).at(i) = list_trees.size()-1;
				}
				else
				{
					treeIDs.at(id_crrIter).at(i) = count -1;
				}
				// std::cout << "new tree added to current list!\n";
			}


		}
		//AS: how do we do this?? Thu May 29 19:10:01 EDT 2014
		//case where the cold chain is on the head node, but on a different chain than the 0 chain
		/*
		if(save==1 && currentid == 0  && chainid != 0 && temp == 1.0) // sampling
		{
			compute_coalTimes_tipIDs(id_crrIter,i,trees_atPrev.at(i));
			logLikelihood.at(id_crrIter).at(i) =logLikelihood_atPrev.at(i);
				

			unsigned int list_size = list_trees.size();
			if(list_size == 0)
			{
				nodeSimple *topo = new nodeSimple;
				std::list<double> ct = coalTimes.at(id_crrIter).at(i);
				ct.sort();
				topo->convert(trees_atPrev.at(i),ct, trees_atPrev.at(i)->size_tree()); // Updated by YC 5/8/2014
				topo->computeSizes();
				list_trees.push_back(topo);
				std::cout << "new list is formed!\n";
			}
			else
			{
				unsigned int found_sampleTopo = 0;
				unsigned int count = 0;
				while(found_sampleTopo == 0 && count < list_size)
				{
					found_sampleTopo = list_trees.at(count)->sameTopo(trees_atPrev.at(i));
					count++;
				}
				if(found_sampleTopo == 0)
				{
					nodeSimple *topo = new nodeSimple;
					std::list<double> ct = coalTimes.at(id_crrIter).at(i);
					ct.sort();
					topo->convert(trees_atPrev.at(i),ct, trees_atPrev.at(i)->size_tree());// updated by YC 5/8/2014
					topo->computeSizes();
					list_trees.push_back(topo);
					treeIDs.at(id_crrIter).at(i) = list_trees.size()-1;
				}
				else
				{
					treeIDs.at(id_crrIter).at(i) = count -1;
				}
				std::cout << "new tree added to current list!\n";
			}


		}
		*/



		if(save == 1) // sampling
		{
	/*		if (currentid !=0 && temp == 1.0) {
				sendingsignal = 1;
				MPI::COMM_WORLD.Send(&sendingsignal, 1, MPI::INT, 0, 4321);
				std::cout << "sendingsignal from process " << currentid << " in iteration " << id_crrIter << "\n";

			}
			if (currentid == 0 && temp != 1.0 && chainid == 0) {

				MPI::COMM_WORLD.Recv(&sendingsignal, 1, MPI::INT, MPI_ANY_SOURCE, 4321);
				std::cout << "sendsignal received\n";
			}
	*/	
			if (currentid !=0 && temp == 1.0) {
				//AS: coal times, tip IDs are computed on the cold process
				compute_coalTimes_tipIDs(id_crrIter,i,trees_atPrev.at(i));
				logLikelihood.at(id_crrIter).at(i) =0.0;//logLikelihood_atPrev.at(i);
				coldprocess = currentid;
				//AS: send cold process ID to the head node
				int cid = currentid;
				MPI::COMM_WORLD.Send(&cid, 1, MPI::INT, 0, 14);
				//AS: send log likelihood
				double llhere = 0.0;//logLikelihood_atPrev.at(i);
				try {
					MPI::COMM_WORLD.Send(&llhere, 1, MPI::DOUBLE, 0, 356);
				} catch (MPI::Exception e) {
					std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
					return -1;				
				}
			} 
			if (currentid == 0 && temp != 1.0 && chainid == 0) {
				//AS: receieve cold process ID on head node
				MPI::COMM_WORLD.Recv(&coldprocess, 1, MPI::INT, MPI_ANY_SOURCE, 14, status);
				//AS: receieve log likelihood on head node*/
				double llhere = 0.0;
				try {
					MPI::COMM_WORLD.Recv(&llhere, 1, MPI::DOUBLE, MPI_ANY_SOURCE, 356, status);
				} catch (MPI::Exception e) {
					std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
					return -1;
				}
				logLikelihood.at(id_crrIter).at(i) = 0.0;//llhere;

				list_size = list_trees.size();
				//std::cout << "Receiving ctsize\n";
				MPI::COMM_WORLD.Recv(&ctsize, 1, MPI::UNSIGNED, coldprocess, 198, status);
				std::list<double> ct;
				//ct.resize(ctsize - 1);
				double ctvalue = 0.0;
				for (unsigned int y = 0; y < ctsize; y++) {
					MPI::COMM_WORLD.Recv(&ctvalue, 1, MPI::DOUBLE, coldprocess, 210, status);
					ct.push_back(ctvalue);
					coalTimes.at(id_crrIter).at(i).push_back(ctvalue);
					ctvalue = 0.0;
				}
				unsigned int ti = 0;
				for (unsigned int y = 0; y < ctsize+1; y++) {
					MPI::COMM_WORLD.Recv(&ti, 1, MPI::UNSIGNED, coldprocess, 211, status);
					tipIDs.at(id_crrIter).at(i).push_back(ti);
					ti = 0;
				}	
				ct.sort();
				MPI::COMM_WORLD.Send(&list_size, 1, MPI::UNSIGNED, coldprocess, 1234);
				nodeSimple *topo = new nodeSimple;
				//AS: On cold process, send the coalescent tree
				topo->MPIreceive_coaltree(ct, coldprocess);
				topo->print_topo();

				//AS: since list size is 0, just compute sizes needed for topology and insert
				topo->computeSizes();
				if (list_size == 0) {
					list_trees.push_back(topo);
				}
				else //AS: list size is not 0! which means there are already topologies existing
				{
					unsigned int found_sampleTopo = 0;
					unsigned int count = 0;
					//AS: Recursively receive topology info and see if same topology is found or not
					while(found_sampleTopo == 0 && count < list_size) {
						found_sampleTopo = list_trees.at(count)->sameTopo(topo);
						//std::cout << "Am I here at all?? receive_topoinfo...\n";
						//found_sampleTopo = list_trees.at(count)->receive_topoinfo(coldprocess);
						count++;
					} 
					//AS: if same topology was not found, then push back current topology
					if(found_sampleTopo == 0)
					{
						list_trees.push_back(topo);
						treeIDs.at(id_crrIter).at(i) = list_trees.size()-1;
					}
					else //AS: same topology was found!
					{
						treeIDs.at(id_crrIter).at(i) = count -1;
					}
				}
			}
			//AS: On cold process, send the size of tree and coalescent times
			if (currentid != 0 && temp == 1.0) {
				std::list<double> ct = coalTimes.at(id_crrIter).at(i);
				
				ctsize = ct.size();
				MPI::COMM_WORLD.Send(&ctsize, 1, MPI::UNSIGNED, 0, 198);
				trees_atPrev.at(i)->MPIsend_coaltimes_tipIDs();
			
				unsigned int list_size = 0;
				unsigned int count = 0;
				MPI::COMM_WORLD.Recv(&list_size, 1, MPI::UNSIGNED, 0, 1234, status);
				trees_atPrev.at(i)->MPIsend_coaltree();
			}
		}
	}

	if(save == 1 && currentid == 0 && temp == 1.0 && chainid == 0) {
		logPrior_trees.at(id_crrIter) = logPriorTrees_atPrev;
	}
	//AS: on cold process, send the logpriortrees
	if (save == 1 &&  temp == 1.0 && currentid != 0) {
			logPrior_trees.at(id_crrIter) = logPriorTrees_atPrev;
		try {
			MPI::COMM_WORLD.Send(&logPriorTrees_atPrev, 1, MPI::DOUBLE, 0, 175);
		} catch (MPI::Exception e) {
			std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
			return -1;
		}
	}
	//AS: on head node, receive it
	if (save == 1 &&  temp != 1.0 && currentid == 0 && chainid == 0) {
		double logprior = 0.0;
		try {
			MPI::COMM_WORLD.Recv(&logprior, 1, MPI::DOUBLE, MPI_ANY_SOURCE, 175, status);
		} catch (MPI::Exception e) {
			std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
			return -1;
		}
		logPrior_trees.at(id_crrIter) = logprior;
	}

	return 1;
}


void Chain::compute_jointPrior_indepLoci()
{
  logPriorTrees_atPrev =0;
  for(unsigned int i=0; i< numSubLoci; i++)
    logPriorTrees_atPrev += logPriorTrees_atPrev_lociParallel.at(i);
  return;
}

void Chain::compute_logLikelihood_indepLoci()
{
  logJointLikelihood_atPrev_lociParallel =0;
  for(unsigned int i=0; i< numSubLoci; i++)
    logJointLikelihood_atPrev_lociParallel += logLikelihood_atPrev.at(i);
  return;
}


void Chain::UpdateTrees_LociInMPI_newProposal(unsigned int id_crrIter, vector<locus> lc, unsigned int save, int currentid, int chainid)
{

  // std::cout << "In Chain::UpdateTrees_LociInMPI_newProposal()\n";
  
  MPI::Status status;

  for(unsigned int i=0; i < numSubLoci; i++)
    {
      
      numTries_trees(0,i)++;

      node* crrTree = trees_atPrev.at(i);
      node* newTree = new node;
      newTree->deepCopy(crrTree);
      
      // get the distance to move
      // double stdv = 1; // std::min(newTree->age/3,20.0); -- YC 9/25/2014
      //double stdv = runiform()*std::min(splittingTimeMax/3,20.0);
      //double stdv = runiform()*5;
      double nSegSites = (double) lc.at(i+locusID_start).get_n_sites_uniq();
      double nSites = (double) lc.at(i+locusID_start).get_nSites();
      double nSeq = (double) lc.at(i+locusID_start).get_nGeneCopies();
      //double stdv = std::max(1.0/3, (double)(nSegSites/(2*nSeq-2))/2 );
      double stdv = std::min(20.0, std::max(1/(2*nSeq-2), nSegSites/(2*nSeq-2) ));
      if(lc.at(0).get_multiLocusSpecific_mutationRate() ==1) // mutation rate per site
	{
	  stdv = std::min(20.0, std::max(0.1/3, nSegSites/nSites/3 ));
	  // stdv = std::min(20.0, std::max(0.1/(2*nSeq-2), nSegSites/nSites/(2*nSeq-2) ));
	}

      double slidedist = rNormalDistribution(0.0, stdv);
      // rNormalDist.push_back(slidedist); // for debugging only
      newTree = newTree->propose_coaltree(slidedist);
      
      /*
      if(i+locusID_start ==0)
	{
	  std::cout <<" locusID =" << i+locusID_start << " stdv = " << stdv 
		    <<" slidedist = " << slidedist <<"\n";
	  // std::cout << "nSegsites = " << nSegSites << " nSites = " <<nSites << " nSeq = " << nSeq <<"\n";
	  std::cout << "old tree = ";
	  crrTree->print_coaltree();
	  std::cout << "new tree = ";
	  newTree->print_coaltree();
	}
      */

      // Computing the likelihood of the proposed tree, newTree
      unsigned int likelihoodModel = lc.at(i+locusID_start).getLikelihoodModel();
      double newLogLik = 0.0;
      try 
	{
	  if(likelihoodModel == 0) // IS
	    newLogLik = newTree->logLikelihood_IS(lc.at(i+locusID_start), mutationScaler_atPrev.at(i),0);
	  else if(likelihoodModel == 2) // HKY
	    newLogLik = newTree->logLikelihood_HKY(lc.at(i+locusID_start), mutationScaler_atPrev.at(i), kappa_atPrev.at(i));
	} 
      catch (std::exception &e) 
	{
	  std::cout << "Can't access lc or kappa to calculate likelihood - vector index out of bounds\n";
	}
      
      
      unsigned int rejectThisTree = 0;
      if(likelihoodModel == 0) // IS model
	if(newTree->get_rejectThisTree_IS() == 1 && std::isnan(newLogLik) == 1)
	 {
	   rejectThisTree = 1;
	   newTree->set_rejectThisTree_IS(0);	 
	 }

       if(rejectThisTree == 0)
	 {
	   // Computing the joint prior distribution of the trees including "newTree"
	   newTree->compute_totalCoalescentRate();
	   double newlogPrior = 0.0;
	   double logRatioPrior = 0.0;
	   if(lociInParallel==1)
	     {
	       // newlogPrior = newTree->compute_prior_indepLoci(lc.at(i+locusID_start),popSizeMax,splittingTimeMax,priorType,changeOfRatePoint);
	       newlogPrior = newTree->compute_prior_indepLoci(lc.at(i+locusID_start),popSizeMax,priorType,changeOfRatePoint);
	       logRatioPrior = newlogPrior - logPriorTrees_atPrev_lociParallel.at(i);
	     }
	   else
	     {
	       newlogPrior = compute_logPriorTrees(newTree,i);
	       logRatioPrior = newlogPrior - logPriorTrees_atPrev;
	     }
	   double MHratio = MHratio_trees_newProposal(i, newTree, newLogLik,logRatioPrior, slidedist );

	   
	   if(runiform() < std::min(1.0, MHratio) )
	     {
	       // accept
	       trees_atPrev.at(i)->deleteCoalTree();
	       trees_atPrev.at(i) = new node;
	       trees_atPrev.at(i)->deepCopy(newTree);
	       logLikelihood_atPrev.at(i)=newLogLik;
	       if(lociInParallel==1)
		 logPriorTrees_atPrev_lociParallel.at(i) = newlogPrior;
	       else
		 logPriorTrees_atPrev = newlogPrior;
	       
	       numAccpt_trees(0,i)++;
	       
	       /*
	       if(i+locusID_start ==0)
		 {
		   // REMOVE
		   std::cout << "Accepted\n";
		 }
	       */
	     }
	   else  // Reject
	     {
	       // REMOVE
	       //std::cout << "reject\n";
	       //trees_atPrev.at(i)->print_coaltree();
	       
	     }
	 }
       newTree->deleteCoalTree();
    }
  return;
}


// YC 9/11/2014
void Chain::UpdateTrees_MPI_newProposal(unsigned int id_crrIter, vector<locus> lc, unsigned int save, int currentid, int chainid)
{
  
  MPI::Status status;
  // double temp = GetTemperature(); -- Redundant! Class "chain" has member "temperature" YC 9/25/2014

  /*
  if(save==1 &&  chainid == 0 &&  currentid == 0)
    {
      treeIDs.at(id_crrIter).resize(n_loci);
      coalTimes.at(id_crrIter).resize(n_loci);
      tipIDs.at(id_crrIter).resize(n_loci);
    }
  */
  // logLikelihood.at(id_crrIter).resize(n_loci);
 
  /* YC 9/25/2014
  if(save ==1 &&  currentid == 0)
    {
      coalTimes.at(id_crrIter).resize(n_loci);
      tipIDs.at(id_crrIter).resize(n_loci);
    }
  */

  for(unsigned int i=0; i < n_loci; i++)
    {
      
      numTries_trees(0,i)++;

      node* crrTree = trees_atPrev.at(i);
      node* newTree = new node;
      newTree->deepCopy(crrTree);
      
      // get the distance to move
      // double stdv = 1; // std::min(newTree->age/3,20.0); -- YC 9/25/2014
      double stdv = std::min(splittingTimeMax/3,20.0);
      double slidedist = rNormalDistribution(0.0, stdv);

      newTree = newTree->propose_coaltree(slidedist);

      // Computing the likelihood of the proposed tree, newTree
      unsigned int likelihoodModel = lc.at(i).getLikelihoodModel();
      double newLogLik = 0.0;
      try 
	{
	  if(likelihoodModel == 0) // IS
	    newLogLik = newTree->logLikelihood_IS(lc.at(i), mutationScaler_atPrev.at(i),0);
	  else if(likelihoodModel == 2) // HKY
	    newLogLik = newTree->logLikelihood_HKY(lc.at(i), mutationScaler_atPrev.at(i), kappa_atPrev.at(i));
	} 
      catch (std::exception &e) 
	{
	  std::cout << "Can't access lc or kappa to calculate likelihood - vector index out of bounds\n";
	}
      
      
      unsigned int rejectThisTree = 0;
      if(likelihoodModel == 0)
	if(newTree->get_rejectThisTree_IS() == 1)
	 {
	   rejectThisTree = 1;
	   newTree->set_rejectThisTree_IS(0);	 
	 }

       if(rejectThisTree == 0)
	 {
	   // Computing the joint prior distribution of the trees including "newTree"
	   newTree->compute_totalCoalescentRate();
	   double newlogPrior = compute_logPriorTrees(newTree,i);
	   double logRatioPrior = newlogPrior - logPriorTrees_atPrev;
	   double MHratio = MHratio_trees_newProposal(i, newTree, newLogLik,logRatioPrior, slidedist );
	   
	   if(runiform() < min(1.0, MHratio) )
	     {
	       // accept
	       trees_atPrev.at(i)->deleteCoalTree();
	       trees_atPrev.at(i) = new node;
	       trees_atPrev.at(i)->deepCopy(newTree);
	       logLikelihood_atPrev.at(i)=newLogLik;
	       logPriorTrees_atPrev = newlogPrior;
	       
	       numAccpt_trees(0,i)++;
	       
	     }
	   else  // Reject
	     {
	       // REMOVE
	       //std::cout << "reject\n";
	       //trees_atPrev.at(i)->print_coaltree();
	       
	     }
	 }
       newTree->deleteCoalTree();
    }
  return;
}


/*
void Chain::samplingUpdatedTrees(unsigned int id_crrIter)
{
	unsigned int nloci = chains.at(chainid).GetNumLoci();

	for (unsigned int i = 0; i < nloci; i++) {
		chains.at(0).compute_coalTimes_tipIDs(id_crrIter, i, chains.at(chainid).GetTreesAtPrev(i));
		chains.at(0).SetLogLikelihood(id_crrIter, i, chains.at(chainid).GetLogLikelihoodAtPrev(i));
		unsigned int list_size = chains.at(0).GetListTrees().size();		
		if (list_size == 0) {
			nodeSimple *topo = new nodeSimple;
			std::list<double> ct = chains.at(0).GetCoalTimes(id_crrIter, i);
			ct.sort();
			topo->convert(chains.at(chainid).GetTreesAtPrev(i),ct, chains.at(chainid).GetTreesAtPrev(i)->size_tree());
			topo->computeSizes();
			chains.at(0).SetListTrees(topo);
		} else {
			unsigned int found_sampleTopo = 0;
			unsigned int count = 0;
			while(found_sampleTopo == 0 && count < list_size) {
				found_sampleTopo = chains.at(0).GetListTrees().at(count)->sameTopo(chains.at(chainid).GetTreesAtPrev(i));
				count++;
			}
			if (found_sampleTopo == 0) {
				nodeSimple *topo = new nodeSimple;
				std::list<double> ct = chains.at(0).GetCoalTimes(id_crrIter, i);
				ct.sort();
				topo->convert(chains.at(chainid).GetTreesAtPrev(i), ct, chains.at(chainid).GetTreesAtPrev(i)->size_tree());
				topo->computeSizes();
				chains.at(0).SetListTrees(topo);
				chains.at(0).SetTreeIDs(id_crrIter, i, chains.at(0).GetListTrees().size() - 1);
			} else {
				chains.at(0).SetTreeIDs(id_crrIter, i, count - 1);
			}	
		}
	}	
	chains.at(0).SetLogPriorsAtPrev(id_crrIter, chains.at(chainid).GetLogPriorPrev());
	chains.at(0).SetKappaAtPrev(id_crrIter, chains.at(chainid).GetKappaAtPrev());
	if (nloci > 1) {
		chains.at(0).SetMutationScalerAtPrev(id_crrIter, chains.at(chainid).GetMutationScalerAtPrev());
	}
	
	return;
	}*/



// YC 9/11/2014
int Chain::UpdateTrees_MPI_newProposal_old(unsigned int id_crrIter, vector<locus> lc, unsigned int save, int currentid, int chainid, int numprocesses, int coldchain)
{
  
  MPI::Status status;
  double temp = GetTemperature(); // Redundant! Class "chain" has member "temperature" YC 9/25/2014

  if(save==1 &&  chainid == 0 &&  currentid == 0)
    {
      treeIDs.at(id_crrIter).resize(n_loci);
      coalTimes.at(id_crrIter).resize(n_loci);
      tipIDs.at(id_crrIter).resize(n_loci);
    }
  logLikelihood.at(id_crrIter).resize(n_loci);
 
  /* YC 9/25/2014
  if(save ==1 &&  currentid == 0)
    {
      coalTimes.at(id_crrIter).resize(n_loci);
      tipIDs.at(id_crrIter).resize(n_loci);
    }
  */

  for(unsigned int i=0; i < n_loci; i++)
    {
      
      numTries_trees(0,i)++;

      node* crrTree = trees_atPrev.at(i);
      node* newTree = new node;
      newTree->deepCopy(crrTree);
      
      // get the distance to move
      // double stdv = 1; // std::min(newTree->age/3,20.0); -- YC 9/25/2014
      double stdv = std::min(splittingTimeMax/3,20.0);
      double slidedist = rNormalDistribution(0.0, stdv);

      newTree = newTree->propose_coaltree(slidedist);

      // Computing the likelihood of the proposed tree, newTree
      double newLogLik = 0.0;
      try 
	{
	  newLogLik = newTree->logLikelihood_HKY(lc.at(i), mutationScaler_atPrev.at(i), kappa_atPrev.at(i));
	} 
      catch (std::exception &e) 
	{
	  std::cout << "Can't access lc or kappa to calculate likelihood - vector index out of bounds\n";
	}
      
      // Computing the joint prior distribution of the trees including "newTree"
      newTree->compute_totalCoalescentRate();
      double newlogPrior = compute_logPriorTrees(newTree,i);
      double logRatioPrior = newlogPrior - logPriorTrees_atPrev;
      double MHratio = MHratio_trees(i, newTree, newLogLik,logRatioPrior, slidedist );

      if(runiform() < min(1.0, MHratio) )
	{
	  // accept
	  trees_atPrev.at(i)->deleteCoalTree();
	  trees_atPrev.at(i) = new node;
	  trees_atPrev.at(i)->deepCopy(newTree);
	  logLikelihood_atPrev.at(i)=newLogLik;
	  logPriorTrees_atPrev = newlogPrior;
	  
	  numAccpt_trees(0,i)++;
	  
	}
      else  // Reject
	{
		  // REMOVE
	  //std::cout << "reject\n";
	  //trees_atPrev.at(i)->print_coaltree();
	  
	}
      newTree->deleteCoalTree();
      
	
		double logl = 0.0;
		double logprior = 0.0;
		unsigned int ctsize = 0;
		unsigned int list_size = 0;
		int coldprocess = 0;
		int sendingsignal = 0;
		temp = GetTemperature();
		//AS: this is the case when the head node, chain 0 is still cold
		if(save==1 && currentid == 0  && chainid == 0 && temp == 1.0) // sampling
		{
			compute_coalTimes_tipIDs(id_crrIter,i,trees_atPrev.at(i));
			logLikelihood.at(id_crrIter).at(i) =logLikelihood_atPrev.at(i);

			unsigned int list_size = list_trees.size();
			if(list_size == 0)
			{
				nodeSimple *topo = new nodeSimple;
				std::list<double> ct = coalTimes.at(id_crrIter).at(i);
				ct.sort();
				topo->convert(trees_atPrev.at(i),ct, trees_atPrev.at(i)->size_tree()); // Updated by YC 5/8/2014
				topo->computeSizes();
				list_trees.push_back(topo);
				//std::cout << "new list is formed!\n";
			}
			else
			{
				unsigned int found_sampleTopo = 0;
				unsigned int count = 0;
				while(found_sampleTopo == 0 && count < list_size)
				{
					found_sampleTopo = list_trees.at(count)->sameTopo(trees_atPrev.at(i));
					count++;
				}
				if(found_sampleTopo == 0)
				{
					nodeSimple *topo = new nodeSimple;
					std::list<double> ct = coalTimes.at(id_crrIter).at(i);
					ct.sort();
					topo->convert(trees_atPrev.at(i),ct, trees_atPrev.at(i)->size_tree());// updated by YC 5/8/2014
					topo->computeSizes();
					list_trees.push_back(topo);
					treeIDs.at(id_crrIter).at(i) = list_trees.size()-1;
				}
				else
				{
					treeIDs.at(id_crrIter).at(i) = count -1;
				}
				//std::cout << "new tree added to current list!\n";
			}


		}

	// Head node does not have the cold chain
		if(save == 1) // sampling
		{
		// For non-head node but with the cold chain
			if (currentid !=0 && temp == 1.0) {
				//AS: coal times, tip IDs are computed on the cold process
				//compute_coalTimes_tipIDs(id_crrIter,i,trees_atPrev.at(i));
				logLikelihood.at(id_crrIter).at(i) =logLikelihood_atPrev.at(i);
				coldprocess = currentid;
				//AS: send cold process ID to the head node
				int cid = currentid;
				MPI::COMM_WORLD.Send(&cid, 1, MPI::INT, 0, 1427);
				//AS: send log likelihood
				double llhere = logLikelihood_atPrev.at(i);
				try {
					MPI::COMM_WORLD.Send(&llhere, 1, MPI::DOUBLE, 0, 356);
				} catch (MPI::Exception e) {
					std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
					return -1;				
				}
			} 
	// For head node 
			if (currentid == 0 && temp != 1.0 && chainid == 0 && numprocesses > 1 && coldchain < 0) {
				//AS: receieve cold process ID on head node
				MPI::COMM_WORLD.Recv(&coldprocess, 1, MPI::INT, MPI_ANY_SOURCE, 1427, status);
				//AS: receieve log likelihood on head node*/
				double llhere = 0.0;
				try {
					MPI::COMM_WORLD.Recv(&llhere, 1, MPI::DOUBLE, MPI_ANY_SOURCE, 356, status);
				} catch (MPI::Exception e) {
					std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
					return -1;
				}
				logLikelihood.at(id_crrIter).at(i) = llhere;

				list_size = list_trees.size();
				//std::cout << "Receiving ctsize\n";
				MPI::COMM_WORLD.Recv(&ctsize, 1, MPI::UNSIGNED, coldprocess, 198, status);
				std::list<double> ct;
				//ct.resize(ctsize - 1);
				double ctvalue = 0.0;
				for (unsigned int y = 0; y < ctsize; y++) {
					MPI::COMM_WORLD.Recv(&ctvalue, 1, MPI::DOUBLE, coldprocess, 298, status);					
					ct.push_back(ctvalue);
					coalTimes.at(id_crrIter).at(i).push_back(ctvalue);
					ctvalue = 0.0;
				}
				unsigned int ti = 0;
				for (unsigned int y = 0; y < ctsize+1; y++) {
					MPI::COMM_WORLD.Recv(&ti, 1, MPI::UNSIGNED, coldprocess, 211, status);
					tipIDs.at(id_crrIter).at(i).push_back(ti);
					ti = 0;
				}	
				ct.sort();

				nodeSimple *topo = new nodeSimple;
				//AS: On cold process, send the coalescent tree
				topo->MPIreceive_coaltree(ct, coldprocess);
				topo->print_topo();

				//AS: since list size is 0, just compute sizes needed for topology and insert
				topo->computeSizes();
				if (list_size == 0) {
					list_trees.push_back(topo);
				}
				else //AS: list size is not 0! which means there are already topologies existing
				{
					unsigned int found_sampleTopo = 0;
					unsigned int count = 0;
					//AS: Recursively receive topology info and see if same topology is found or not
					while(found_sampleTopo == 0 && count < list_size) {
						found_sampleTopo = list_trees.at(count)->sameTopo(topo);
						count++;
					} 
					//AS: if same topology was not found, then push back current topology
					if(found_sampleTopo == 0)
					{
						list_trees.push_back(topo);
						treeIDs.at(id_crrIter).at(i) = list_trees.size()-1;
					}
					else //AS: same topology was found!
					{
						treeIDs.at(id_crrIter).at(i) = count -1;
					}
				}
			}
			//AS: On cold process, send the size of tree and coalescent times
			if (currentid != 0 && temp == 1.0) {
				ctsize = trees_atPrev.at(i)->size_tree()-1;
				MPI::COMM_WORLD.Send(&ctsize, 1, MPI::UNSIGNED, 0, 198);
				trees_atPrev.at(i)->MPIsend_coaltimes_tipIDs();
			
				//AS: Send coal_tree	
				trees_atPrev.at(i)->MPIsend_coaltree();
				//trees_atPrev.at(i)->send_tipids();
			}
		}
		//	std::cout << "Done with locus " << i << " on process " << currentid<< " \n\n";
	}

	if(save == 1 && currentid == 0 && temp == 1.0 && chainid == 0 && coldchain >= 0) {
		logPrior_trees.at(id_crrIter) = logPriorTrees_atPrev;
	}
	//AS: on cold process, send the logpriortrees
	if (save == 1 &&  temp == 1.0 && currentid != 0) {
			logPrior_trees.at(id_crrIter) = logPriorTrees_atPrev;
		try {
			MPI::COMM_WORLD.Send(&logPriorTrees_atPrev, 1, MPI::DOUBLE, 0, 1765);
		} catch (MPI::Exception e) {
			std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
			return -1;
		}
	}
	//AS: on head node, receive it
	if (save == 1 &&  temp != 1.0 && currentid == 0 && chainid == 0 && numprocesses > 1 && coldchain < 0) {
		double logprior = 0.0;
		try {
			MPI::COMM_WORLD.Recv(&logprior, 1, MPI::DOUBLE, MPI_ANY_SOURCE, 1765, status);
		} catch (MPI::Exception e) {
			std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
			return -1;
		}
		logPrior_trees.at(id_crrIter) = logprior;
	}

	return 1;
}

/***
 * Update the ratio between two mutation rate scalars.
 * The upate is drawn from a uniform log scale between 1/maxratio and maxratio.
 * maxratio is set to be 3 times the maximum value of the individual scalars
 */
void Chain::update_randomPairsMutationScalars_lociInParallel(unsigned int id_crrIter, std::vector<locus> lc, unsigned int crrProcID, unsigned int nProc)
{
  // std::cout << "\nIn Chain::update_randomPairsMutationScalars_lociInParallel\n";


  pairIDsMut.resize(2*nPairsMut);
  std::vector<unsigned int> IDs; 
  for(unsigned int ll =0; ll<n_loci; ll++)
    IDs.push_back(ll);
  std::random_shuffle(IDs.begin(),IDs.end()); // shuffle elemenets in IDs
  unsigned int id =0;
  for(unsigned int j=0; j< 2*nPairsMut; j++)	
    {	  
      if(crrProcID==0)
	id = IDs.at(j);
      MPI::COMM_WORLD.Barrier();
      MPI::COMM_WORLD.Bcast(&id, 1, MPI::UNSIGNED, 0);
      MPI::COMM_WORLD.Barrier();      
      pairIDsMut.at(j) = id;
    }



  // static int start = 0;
  /* windowsize and maxratio are on log scales
     all mutation rate scalars have the same priors and windowsize */
  // if (start == 0)
  // {
  /* 5/27/2011 changed the range expansion factor from 3 to 2,
		     *  based on modeling that showed 2 is sufficient */
  /* 11/6/2011  changed back to 3.  no point changing to 2,
   *  and makes new results inconsistent with old results */
  
  // Modified by YC 6/9/2014
  // See IMa2 codes("update_mc_params.cpp") for how to determine the maximum ratio and window size.
  double maxRatio = 3.0*log(10000); //3.0 * L[0].u_rec[0].pr.max; // the upper bound of the mutation rate - YC 3/6/2014
  double windowSize = log(10000);
  //double windowSize = log(10000)/n_loci;
		//    if (calcoptions[MUTATIONPRIORRANGE])
		//      windowsize = L[0].u_rec[0].pr.max;
		//    else
		//      windowsize = L[0].u_rec[0].win; // window width, may be used for updating - YC 3/6/2014
		//    start = 1;
	//}

  /*
  std::cout << "nPairsMut = " << nPairsMut 
	    << " pairIDsMut.size() = " << pairIDsMut.size()
	    << " pairIDsMut.at(0).size() = " << pairIDsMut.at(0).size() <<"\n";
  */

  unsigned int id1, id2;
  for(unsigned int i=0; i<nPairsMut; i++)
    {
      // two loci to update
      id1 = pairIDsMut.at(2*i);
      id2 = pairIDsMut.at(2*i+1);
      // Note: after each iteration, the first vector, pairIDsMut.at(0), is removed.
      // Therefore, the vector containing ids to update at the next iteration is pairIDsMut.at(0).

      //std::cout << " crrProcID = " << crrProcID
      //		<< "id1 = " << id1 << " id2 = " << id2 <<"\n";

       // determine the CPU IDs to update
      unsigned int procID1_toUpdate = 0;
      unsigned int procID2_toUpdate = 0;
      unsigned int nSubLoci_smallest = static_cast<unsigned int>(n_loci/nProc);
      unsigned int nRemainder= n_loci -nProc*nSubLoci_smallest;
      unsigned int nProcs_withSmallNumSubLoci = nProc-nRemainder;
      unsigned int nLoci_onProcsWithFewerLoci = nSubLoci_smallest*nProcs_withSmallNumSubLoci;
      if(id1 <nLoci_onProcsWithFewerLoci)
	{
	  procID1_toUpdate = static_cast<unsigned int>(id1/nSubLoci_smallest);
	}
      else
	{
	  procID1_toUpdate =  nProcs_withSmallNumSubLoci+static_cast<unsigned int>((id1-nLoci_onProcsWithFewerLoci)/(nSubLoci_smallest+1));
	}
      if(id2 < nLoci_onProcsWithFewerLoci)
	{
	  procID2_toUpdate = static_cast<unsigned int>(id2/nSubLoci_smallest);
	}
      else
	{
	  procID2_toUpdate =  nProcs_withSmallNumSubLoci+static_cast<unsigned int>((id2-nLoci_onProcsWithFewerLoci)/(nSubLoci_smallest+1));
	}

      //std::cout << " id_crrIter = " << id_crrIter
      //		<< " procID1_toUpdate = "<< procID1_toUpdate 
      //		<<" procID2_toUpdate = " << procID2_toUpdate <<"\n";

      // update mutation scalars
      if(procID1_toUpdate == procID2_toUpdate)
	{
	  if(crrProcID == procID1_toUpdate)
	    {
	      double oldScaler1 = mutationScaler_atPrev.at(id1 - locusID_start); 
	      double oldScaler2 = mutationScaler_atPrev.at(id2 - locusID_start);
	      double logRatio = log (oldScaler1 / oldScaler2);	
	      // Propose a new logRatio
	      double U = runiform ();
	      double newlogRatio = 0.0;
	      if (U > 0.5)
		newlogRatio = logRatio + (2.0 * U - 1.0) * windowSize;
	      else
		newlogRatio = logRatio - windowSize * U * 2.0;
	  
	      if (newlogRatio > maxRatio)
		newlogRatio = 2.0 * maxRatio - newlogRatio;
	      else if (newlogRatio < -maxRatio)
		newlogRatio = 2.0 * (-maxRatio) - newlogRatio;
	      
	      double multiplier = exp ((newlogRatio - logRatio) / 2);
	      std::vector<double> newScalers;
	      newScalers.push_back(oldScaler1 * multiplier);
	      newScalers.push_back(oldScaler2 / multiplier);
	      
	      vector<double>  newlogLik;
	      newlogLik.resize(2);
	      // Compute the new likelihood values
	      newlogLik.at(0) = trees_atPrev.at(id1 -locusID_start)
		->logLikelihood_IS(lc.at(id1), newScalers.at(0),1);
	      newlogLik.at(1) = trees_atPrev.at(id2 -locusID_start)
		->logLikelihood_IS(lc.at(id2), newScalers.at(1),1);
	      double diff_logLik = newlogLik.at(0) - logLikelihood_atPrev.at(id1-locusID_start)
		+ newlogLik.at(1) - logLikelihood_atPrev.at(id2-locusID_start);

	      /// Metropolis-Hastings ratio
	      double MHratio = exp(temperature * diff_logLik);
	      
	      if (runiform() < MHratio)
		{
		  // std::cout << "id1 = " <<id1 << " id2 = " << id2 <<"\n";
		  // Accept
		  mutationScaler_atPrev.at(id1-locusID_start) = newScalers.at(0);
		  logLikelihood_atPrev.at(id1-locusID_start) = newlogLik.at(0);
		  mutationScaler_atPrev.at(id2-locusID_start) = newScalers.at(1);
		  logLikelihood_atPrev.at(id2-locusID_start) = newlogLik.at(1);
		  /*
		  for(unsigned int i=0; i<2; i++)
		    {
		      mutationScaler_atPrev.at(IDs.at(i)-locusID_start) = newScalers.at(i);
		      logLikelihood_atPrev.at(IDs.at(i)-locusID_start) = newlogLik.at(i);
		    }
		  */
		}
	    }// END of if(crrProcID == procID1_toUpdate)
	}// END of if(procID1_toUpdate == procID2_toUpdate)
      else
	{
	  double oldScaler1 =0; double oldScaler2 =0;
	  double oldLogLik1 = 0.0; double oldLogLik2 = 0.0;
	  double newScalar1 = 0.0; double newScalar2 = 0.0;
	  double newLogLik1 = 0.0; double newLogLik2 = 0.0;
	  unsigned int accept =0;
	  if(crrProcID == procID1_toUpdate) // sender
	    {
	      oldScaler1 =mutationScaler_atPrev.at(id1 - locusID_start); 
	      MPI::COMM_WORLD.Send(&oldScaler1, 1, MPI::DOUBLE, procID2_toUpdate, 843718);
	      // std::cout << "sender crrProcID = " << crrProcID << " odlScaler1 = " << oldScaler1 <<"\n";
	      
	      oldLogLik1 = logLikelihood_atPrev.at(id1-locusID_start);
	      MPI::COMM_WORLD.Send(&oldLogLik1, 1, MPI::DOUBLE, procID2_toUpdate, 1798);
	    }
	  else if(crrProcID == procID2_toUpdate) // receiver
	    {
	      oldScaler2 =mutationScaler_atPrev.at(id2 - locusID_start); 
	      MPI::COMM_WORLD.Recv(&oldScaler1, 1, MPI::DOUBLE, procID1_toUpdate, 843718);
	      // std::cout << "receiver crrProcID = " << crrProcID << " odlScaler1 = " << oldScaler1 <<"\n";
	      
	      oldLogLik2 = logLikelihood_atPrev.at(id2-locusID_start);
	      MPI::COMM_WORLD.Recv(&oldLogLik1, 1, MPI::DOUBLE, procID1_toUpdate, 1798);
	      
	      double logRatio = log (oldScaler1 / oldScaler2);	
	      // Propose a new logRatio
	      double U = runiform ();
	      double newlogRatio = 0.0;
	      if (U > 0.5)
		newlogRatio = logRatio + (2.0 * U - 1.0) * windowSize;
	      else
		newlogRatio = logRatio - windowSize * U * 2.0;
	      
	      if (newlogRatio > maxRatio)
		newlogRatio = 2.0 * maxRatio - newlogRatio;
	      else if (newlogRatio < -maxRatio)
		newlogRatio = 2.0 * (-maxRatio) - newlogRatio;
	      
	      double multiplier = exp ((newlogRatio - logRatio) / 2);
	      newScalar1 = oldScaler1 * multiplier;
	      newScalar2 = oldScaler2 / multiplier;
	      
	      // std::cout << "newSaclar1= " << newScalar1 <<" newScalar2= " << newScalar2  <<"\n";
	      
	      // Sending newScalar1
	      MPI::COMM_WORLD.Send(&newScalar1, 1, MPI::DOUBLE, procID1_toUpdate,37551798);  
	    }
	  MPI::COMM_WORLD.Barrier();
	  
	  // computing new likelihoods
	  if(crrProcID == procID1_toUpdate) 
	    {	  
	      // receiving newScalar1
	      MPI::COMM_WORLD.Recv(&newScalar1, 1, MPI::DOUBLE, procID2_toUpdate,37551798);
	      
	      //std::cout << "crrProcID = " << crrProcID << " computing lik.."
	      //	    << "IDs.at(0) = " << IDs.at(0) << "locusID_start = "<< locusID_start <<"\n";
	      newLogLik1 = trees_atPrev.at(id1 -locusID_start)
		->logLikelihood_IS(lc.at(id1), newScalar1,1);
	      //std::cout << "crrProcID = " << crrProcID << "newloglik1=" << newLogLik1 <<"\n";
	      
	      // sending newLogLik1
	      MPI::COMM_WORLD.Send(&newLogLik1, 1, MPI::DOUBLE, procID2_toUpdate,30537);  
	    }
	  else if(crrProcID == procID2_toUpdate)
	    {
	      // receiving newLogLik1
	      MPI::COMM_WORLD.Recv(&newLogLik1, 1, MPI::DOUBLE, procID1_toUpdate, 30537);
	      
	      //std::cout << "crrProcID = " << crrProcID << " computing lik.."
	      //	    << "IDs.at(1) = " << IDs.at(1) << "locusID_start = "<< locusID_start <<"\n";
	      newLogLik2 = trees_atPrev.at(id2 -locusID_start)
		->logLikelihood_IS(lc.at(id2), newScalar2,1);
	      // std::cout << "crrProcID = " << crrProcID << "newloglik2=" << newLogLik2 <<"\n";
	      double diff_logLik = newLogLik1 - oldLogLik1 + newLogLik2 - oldLogLik2;
	      
	      double MHratio = exp(temperature * diff_logLik);
	      if (runiform() < MHratio)
		{
		  // Accept
		  accept = 1;
		}
	      MPI::COMM_WORLD.Send(&accept, 1, MPI::UNSIGNED, procID1_toUpdate, 7252015);
	      
	      if(accept ==1)
		{
		  // std::cout << "id1 = " <<id1 << " id2 = " << id2 <<"\n";

		  mutationScaler_atPrev.at(id2-locusID_start) = newScalar2;
		  logLikelihood_atPrev.at(id2-locusID_start) = newLogLik2;	     
		}	 
	    }
	  MPI::COMM_WORLD.Barrier();
	  
	  // MH update result
	  if(crrProcID == procID1_toUpdate) 
	    {
	      // receive accept
	      MPI::COMM_WORLD.Recv(&accept, 1, MPI::UNSIGNED, procID2_toUpdate, 7252015);
	      
	      if(accept==1)
		{
		  mutationScaler_atPrev.at(id1-locusID_start) = newScalar1;
		  logLikelihood_atPrev.at(id1-locusID_start) = newLogLik1;	      
		}
	    }
	}// END of if(procID1_toUpdate != procID2_toUpdate)
      MPI::COMM_WORLD.Barrier();
      
      
    }

  // Remove the vector of ids that have updated.
  // pairIDsMut.erase(pairIDsMut.begin());
  pairIDsMut.resize(0);

  return;
}

/***
 * Update the ratio between two mutation rate scalars.
 * The upate is drawn from a uniform log scale between 1/maxratio and maxratio.
 * maxratio is set to be 3 times the maximum value of the individual scalars
 */
void Chain::update_mutationScalar_lociInParallel(unsigned int id_crrIter, std::vector<locus> lc, unsigned int crrProcID, unsigned int nProc)
{
  // std::cout << "In Chain::update_mutationScaler_Kappa_lociInParallel() \n";

  /// Randomly pick two loci to update
  vector<unsigned int> IDs;
  unsigned int id;
  if(crrProcID ==0)
    {
      id = runiform_discrete(n_loci);
    }
  MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Bcast(&id, 1, MPI::UNSIGNED, 0);
  MPI::COMM_WORLD.Barrier();
  IDs.push_back(id);
  if(crrProcID ==0)
    {
      id = runiform_discrete(n_loci);
      while (id == IDs.at(0))
	id = runiform_discrete(n_loci);
    }
  MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Bcast(&id, 1, MPI::UNSIGNED, 0);
  MPI::COMM_WORLD.Barrier();
  IDs.push_back(id);
  
  // std::cout << "crrProcID = " << crrProcID << " IDs.at(0) = " << IDs.at(0)<<" IDs.at(1)=" << IDs.at(1) <<"\n";
  
  
  // static int start = 0;
  /* windowsize and maxratio are on log scales
     all mutation rate scalars have the same priors and windowsize */
  // if (start == 0)
  // {
  /* 5/27/2011 changed the range expansion factor from 3 to 2,
		     *  based on modeling that showed 2 is sufficient */
  /* 11/6/2011  changed back to 3.  no point changing to 2,
   *  and makes new results inconsistent with old results */
  
  // Modified by YC 6/9/2014
  // See IMa2 codes("update_mc_params.cpp") for how to determine the maximum ratio and window size.
  double maxRatio = 3.0*log(10000); //3.0 * L[0].u_rec[0].pr.max; // the upper bound of the mutation rate - YC 3/6/2014
  double windowSize = log(10000);
  //double windowSize = log(10000)/n_loci;
		//    if (calcoptions[MUTATIONPRIORRANGE])
		//      windowsize = L[0].u_rec[0].pr.max;
		//    else
		//      windowsize = L[0].u_rec[0].win; // window width, may be used for updating - YC 3/6/2014
		//    start = 1;
	//}

  // determine the CPU IDs to update
  unsigned int procID1_toUpdate = 0;
  unsigned int procID2_toUpdate = 0;
  unsigned int nSubLoci_smallest = static_cast<unsigned int>(n_loci/nProc);
  unsigned int nRemainder= n_loci -nProc*nSubLoci_smallest;
  unsigned int nProcs_withSmallNumSubLoci = nProc-nRemainder;
  unsigned int nLoci_onProcsWithFewerLoci = nSubLoci_smallest*nProcs_withSmallNumSubLoci;
  if(IDs.at(0) <nLoci_onProcsWithFewerLoci)
    {
      procID1_toUpdate = static_cast<unsigned int>(IDs.at(0)/nSubLoci_smallest);
    }
  else
    {
      procID1_toUpdate =  nProcs_withSmallNumSubLoci+static_cast<unsigned int>((IDs.at(0)-nLoci_onProcsWithFewerLoci)/(nSubLoci_smallest+1));
    }
  if(IDs.at(1) < nLoci_onProcsWithFewerLoci)
    {
      procID2_toUpdate = static_cast<unsigned int>(IDs.at(1)/nSubLoci_smallest);
    }
  else
    {
      procID2_toUpdate =  nProcs_withSmallNumSubLoci+static_cast<unsigned int>((IDs.at(1)-nLoci_onProcsWithFewerLoci)/(nSubLoci_smallest+1));
    }

  // std::cout << "crrProcID = " << crrProcID 
  //	    << " procID1_toUpdate = " << procID1_toUpdate <<" procID2_toUpdate=" << procID2_toUpdate <<"\n";

  MPI::COMM_WORLD.Barrier();
  if(procID1_toUpdate == procID2_toUpdate)
    {
      if(crrProcID == procID1_toUpdate)
	{
	  double oldScaler1 = mutationScaler_atPrev.at(IDs.at(0) - locusID_start); 
	  double oldScaler2 = mutationScaler_atPrev.at(IDs.at(1) - locusID_start);
	  double logRatio = log (oldScaler1 / oldScaler2);	
	  // Propose a new logRatio
	  double U = runiform ();
	  double newlogRatio = 0.0;
	  if (U > 0.5)
	    newlogRatio = logRatio + (2.0 * U - 1.0) * windowSize;
	  else
	    newlogRatio = logRatio - windowSize * U * 2.0;
	  
	  if (newlogRatio > maxRatio)
	    newlogRatio = 2.0 * maxRatio - newlogRatio;
	  else if (newlogRatio < -maxRatio)
	    newlogRatio = 2.0 * (-maxRatio) - newlogRatio;
	  
	  double multiplier = exp ((newlogRatio - logRatio) / 2);
	  std::vector<double> newScalers;
	  newScalers.push_back(oldScaler1 * multiplier);
	  newScalers.push_back(oldScaler2 / multiplier);
 
	  vector<double>  newlogLik;
	  newlogLik.resize(2);
	  double diff_logLik = 0.0;
	  for(unsigned int i=0; i<2; i++)
	    {
	      // Compute the new likelihood values
	      newlogLik.at(i) = trees_atPrev.at(IDs.at(i) -locusID_start)
		->logLikelihood_IS(lc.at(IDs.at(i)), newScalers.at(i),1);
	      diff_logLik += newlogLik.at(i) - logLikelihood_atPrev.at(IDs.at(i)-locusID_start);
	    }
	  /// Metropolis-Hastings ratio
	  double MHratio = exp(temperature * diff_logLik);
	  
	  if (runiform() < MHratio)
	    {
	      // Accept
	      for(unsigned int i=0; i<2; i++)
		{
		  mutationScaler_atPrev.at(IDs.at(i)-locusID_start) = newScalers.at(i);
		  logLikelihood_atPrev.at(IDs.at(i)-locusID_start) = newlogLik.at(i);
		}
	    }
	}// END of if(crrProcID == procID1_toUpdate)
    }// END of if(procID1_toUpdate == procID2_toUpdate)
  else
    {
      double oldScaler1 =0; double oldScaler2 =0;
      double oldLogLik1 = 0.0; double oldLogLik2 = 0.0;
      double newScalar1 = 0.0; double newScalar2 = 0.0;
      double newLogLik1 = 0.0; double newLogLik2 = 0.0;
      unsigned int accept =0;
      if(crrProcID == procID1_toUpdate) // sender
	{
	  oldScaler1 =mutationScaler_atPrev.at(IDs.at(0) - locusID_start); 
	  MPI::COMM_WORLD.Send(&oldScaler1, 1, MPI::DOUBLE, procID2_toUpdate, 843718);
	  // std::cout << "sender crrProcID = " << crrProcID << " odlScaler1 = " << oldScaler1 <<"\n";

	  oldLogLik1 = logLikelihood_atPrev.at(IDs.at(0)-locusID_start);
	  MPI::COMM_WORLD.Send(&oldLogLik1, 1, MPI::DOUBLE, procID2_toUpdate, 1798);
	}
      else if(crrProcID == procID2_toUpdate) // receiver
	{
	  oldScaler2 =mutationScaler_atPrev.at(IDs.at(1) - locusID_start); 
	  MPI::COMM_WORLD.Recv(&oldScaler1, 1, MPI::DOUBLE, procID1_toUpdate, 843718);
	  // std::cout << "receiver crrProcID = " << crrProcID << " odlScaler1 = " << oldScaler1 <<"\n";

	  oldLogLik2 = logLikelihood_atPrev.at(IDs.at(1)-locusID_start);
	  MPI::COMM_WORLD.Recv(&oldLogLik1, 1, MPI::DOUBLE, procID1_toUpdate, 1798);
	  
	  double logRatio = log (oldScaler1 / oldScaler2);	
	  // Propose a new logRatio
	  double U = runiform ();
	  double newlogRatio = 0.0;
	  if (U > 0.5)
	    newlogRatio = logRatio + (2.0 * U - 1.0) * windowSize;
	  else
	    newlogRatio = logRatio - windowSize * U * 2.0;
	  
	  if (newlogRatio > maxRatio)
	    newlogRatio = 2.0 * maxRatio - newlogRatio;
	  else if (newlogRatio < -maxRatio)
	    newlogRatio = 2.0 * (-maxRatio) - newlogRatio;
	  
	  double multiplier = exp ((newlogRatio - logRatio) / 2);
	  newScalar1 = oldScaler1 * multiplier;
	  newScalar2 = oldScaler2 / multiplier;

	  // std::cout << "newSaclar1= " << newScalar1 <<" newScalar2= " << newScalar2  <<"\n";
	  
	  // Sending newScalar1
	  MPI::COMM_WORLD.Send(&newScalar1, 1, MPI::DOUBLE, procID1_toUpdate,37551798);  
	}
      MPI::COMM_WORLD.Barrier();

      // computing new likelihoods
      if(crrProcID == procID1_toUpdate) 
	{	  
	  // receiving newScalar1
	  MPI::COMM_WORLD.Recv(&newScalar1, 1, MPI::DOUBLE, procID2_toUpdate,37551798);

	  //std::cout << "crrProcID = " << crrProcID << " computing lik.."
	  //	    << "IDs.at(0) = " << IDs.at(0) << "locusID_start = "<< locusID_start <<"\n";
	  newLogLik1 = trees_atPrev.at(IDs.at(0) -locusID_start)
	    ->logLikelihood_IS(lc.at(IDs.at(0)), newScalar1,1);
	  //std::cout << "crrProcID = " << crrProcID << "newloglik1=" << newLogLik1 <<"\n";

	  // sending newLogLik1
	  MPI::COMM_WORLD.Send(&newLogLik1, 1, MPI::DOUBLE, procID2_toUpdate,30537);  
	}
      else if(crrProcID == procID2_toUpdate)
	{
	  // receiving newLogLik1
	  MPI::COMM_WORLD.Recv(&newLogLik1, 1, MPI::DOUBLE, procID1_toUpdate, 30537);

	  //std::cout << "crrProcID = " << crrProcID << " computing lik.."
	  //	    << "IDs.at(1) = " << IDs.at(1) << "locusID_start = "<< locusID_start <<"\n";
	  newLogLik2 = trees_atPrev.at(IDs.at(1) -locusID_start)
	    ->logLikelihood_IS(lc.at(IDs.at(1)), newScalar2,1);
	  // std::cout << "crrProcID = " << crrProcID << "newloglik2=" << newLogLik2 <<"\n";
	  double diff_logLik = newLogLik1 - oldLogLik1 + newLogLik2 - oldLogLik2;
	  
	  double MHratio = exp(temperature * diff_logLik);
	 if (runiform() < MHratio)
	   {
	     // Accept
	     accept = 1;
	   }
	 MPI::COMM_WORLD.Send(&accept, 1, MPI::UNSIGNED, procID1_toUpdate, 7252015);

	 if(accept ==1)
	   {
	     mutationScaler_atPrev.at(IDs.at(1)-locusID_start) = newScalar2;
	     logLikelihood_atPrev.at(IDs.at(1)-locusID_start) = newLogLik2;	     
	   }	 
	}
      MPI::COMM_WORLD.Barrier();
      
      // MH update result
      if(crrProcID == procID1_toUpdate) 
	{
	  // receive accept
	  MPI::COMM_WORLD.Recv(&accept, 1, MPI::UNSIGNED, procID2_toUpdate, 7252015);
	  
	  if(accept==1)
	    {
	      mutationScaler_atPrev.at(IDs.at(0)-locusID_start) = newScalar1;
	      logLikelihood_atPrev.at(IDs.at(0)-locusID_start) = newLogLik1;	      
	    }
	}
    }// END of if(procID1_toUpdate != procID2_toUpdate)
  MPI::COMM_WORLD.Barrier();
  
  // std::cout << "End of Chain::update_mutationScaler_Kappa_lociInParallel() \n";

  return;
}


/***
 * Update the ratio between two mutation rate scalars.
 * The upate is drawn from a uniform log scale between 1/maxratio and maxratio.
 * maxratio is set to be 3 times the maximum value of the individual scalars
 */
void Chain::update_mutationScaler_Kappa(unsigned int id_crrIter, std::vector<locus> lc, unsigned int crrProcID, unsigned int nProc)
{
  // std::cout << "In Chain::update_mutationScaler_Kappa() \n";

  if(lociInParallel ==1)
    {
      update_randomPairsMutationScalars_lociInParallel(id_crrIter,lc,crrProcID, nProc);
      // update_mutationScalar_lociInParallel(id_crrIter,lc,crrProcID, nProc);
    }
  else
    {
      /// Randomly pick two loci to update
      vector<unsigned int> IDs;
      IDs.push_back(runiform_discrete(n_loci));
      unsigned int id = runiform_discrete(n_loci);
      while (id == IDs.at(0))
	id = runiform_discrete(n_loci);
      IDs.push_back(id);
      
      // static int start = 0;
      /* windowsize and maxratio are on log scales
	 all mutation rate scalars have the same priors and windowsize */
      // if (start == 0)
      // {
      /* 5/27/2011 changed the range expansion factor from 3 to 2,
       *  based on modeling that showed 2 is sufficient */
      /* 11/6/2011  changed back to 3.  no point changing to 2,
       *  and makes new results inconsistent with old results */
      
      // Modified by YC 6/9/2014
      // See IMa2 codes("update_mc_params.cpp") for how to determine the maximum ratio and window size.
      double maxRatio = 3.0*log(10000); //3.0 * L[0].u_rec[0].pr.max; // the upper bound of the mutation rate - YC 3/6/2014
      double windowSize = log(10000)/n_loci;
      //    if (calcoptions[MUTATIONPRIORRANGE])
      //      windowsize = L[0].u_rec[0].pr.max;
      //    else
      //      windowsize = L[0].u_rec[0].win; // window width, may be used for updating - YC 3/6/2014
      //    start = 1;
      //}
      
      double oldScaler1 = 0.0, oldScaler2 = 0.0;
      try{
	oldScaler1 = mutationScaler_atPrev.at(IDs.at(0)); 
	oldScaler2 = mutationScaler_atPrev.at(IDs.at(1));
      } catch (std::exception &e) {
	std::cout << "In Chain::update_mutationScaler_Kappa(). "
	  "Can't access a member of mutationScaler - vector index out of bounds\n";
      }
      double logRatio = log (oldScaler1 / oldScaler2);
      
  // Propose a new logRatio
      double U = runiform ();
      double newlogRatio = 0.0;
      if (U > 0.5)
    newlogRatio = logRatio + (2.0 * U - 1.0) * windowSize;
      else
	newlogRatio = logRatio - windowSize * U * 2.0;
      
      if (newlogRatio > maxRatio)
	newlogRatio = 2.0 * maxRatio - newlogRatio;
      else if (newlogRatio < -maxRatio)
	newlogRatio = 2.0 * (-maxRatio) - newlogRatio;
      
      double multiplier = exp ((newlogRatio - logRatio) / 2);
      std::vector<double> newScalers;
      newScalers.push_back(oldScaler1 * multiplier);
      newScalers.push_back(oldScaler2 / multiplier);
      
      unsigned int updateKappa = 0;
      if(lc.at(IDs.at(0)).getLikelihoodModel() == 2 && lc.at(IDs.at(1)).getLikelihoodModel() == 2) // HKY
	{
      updateKappa =1;
	}
      /// Propose new kappa values for the two loci
      vector<double> newKappa, newlogLik;
      newKappa.resize(2); newlogLik.resize(2);
      // Modified by YC 6/9/2014
      // Refer IMa2 codes ("initialize.cpp" and "update_mc_param.cpp") to see how to determine the window size and the upper bound.
      double windowSize_kappa = 2.0;
      double upperBound_kappa = 100.0;//2.0;
      double diff_logLik = 0.0;
      for(unsigned int i=0; i<2; i++)
	{
	  if(updateKappa ==1 )
	    {
	      double U = runiform ();
	      if (U > 0.5)
		{
		  try{
		    newKappa.at(i) = kappa_atPrev.at(IDs.at(i)) + (2.0 * U - 1.0) * windowSize_kappa; // L[li].kappa_rec->win;
		  } catch (std::exception &e) {
		    std::cout << "In Chain::update_mutationScaler_Kappa(). "
		      "Can't access a member of kappa - vector index out of bounds\n";
		  }
		  if (newKappa.at(i) > upperBound_kappa)// L[li].kappa_rec->pr.max)
		    newKappa.at(i) = 2.0 * upperBound_kappa - newKappa.at(i);// L[li].kappa_rec->pr.max - newkappa[i];
		}
	      else
		{
		  try{
		    newKappa.at(i) = kappa_atPrev.at(IDs.at(i)) - windowSize_kappa*U*2.0; //L[li].kappa_rec->win * U * 2.0;
		  } catch (std::exception &e) {
		    std::cout << "In Chain::update_mutationScaler_Kappa(). "
		      "Can't access a member of kappa - vector index out of bounds\n";
		  }
		  if (newKappa.at(i) < 0)
		    newKappa.at(i) = -newKappa.at(i);
		}
	      
	    }
	  // Compute the new likelihood values
	  if(updateKappa == 0) // other models like IS
	    newlogLik.at(i) = trees_atPrev.at(IDs.at(i))
	      ->logLikelihood_IS(lc.at(IDs.at(i)),newScalers.at(i),1);
	  if(updateKappa == 1) // HKY
	    newlogLik.at(i) = trees_atPrev.at(IDs.at(i))
	      ->logLikelihood_HKY(lc.at(IDs.at(i)),newScalers.at(i),newKappa.at(i));
	  diff_logLik += newlogLik.at(i) - logLikelihood_atPrev.at(IDs.at(i));
	}
      /// Metropolis-Hastings ratio
      double MHratio = exp(temperature * diff_logLik);
      
      if (runiform() < MHratio)
	{
	  // Accept
	  for(unsigned int i=0; i<2; i++)
	    {
	      if(updateKappa ==1 ) // HKY
		kappa_atPrev.at(IDs.at(i)) = newKappa.at(i);
	      mutationScaler_atPrev.at(IDs.at(i)) = newScalers.at(i);
	      logLikelihood_atPrev.at(IDs.at(i)) = newlogLik.at(i);
	    }
	}
      else
	{
      // Reject
      //std::cout << "Reject the proposed mutation scalers and kappa values\n";
	}
    }
  
  return;
}


/***
 * Update the ratio between two mutation rate scalars.
 * The upate is drawn from a uniform log scale between 1/maxratio and maxratio.
 * maxratio is set to be 3 times the maximum value of the individual scalars
 */
void Chain::update_mutationScaler_Kappa_usePriorsOnly(unsigned int id_crrIter, std::vector<locus> lc) //(unsigned int id_crrIter, std::vector<locus> lc,unsigned int save, int currentid, int doisendinfo, int chainid)
{
	/// Randomly pick two loci to update
	vector<unsigned int> IDs;
	IDs.push_back(runiform_discrete(n_loci));
	unsigned int id = runiform_discrete(n_loci);
	while (id == IDs.at(0))
		id = runiform_discrete(n_loci);
	IDs.push_back(id);

	// static int start = 0;
	/* windowsize and maxratio are on log scales
	   all mutation rate scalars have the same priors and windowsize */
	// if (start == 0)
	// {
		    /* 5/27/2011 changed the range expansion factor from 3 to 2,
		     *  based on modeling that showed 2 is sufficient */
		    /* 11/6/2011  changed back to 3.  no point changing to 2,
		     *  and makes new results inconsistent with old results */
	// Modified by YC 6/9/2014
	// See IMa2 codes("update_mc_params.cpp"). How to determine the maximum ratio and window size.
	double maxRatio =  3.0*log(10000); //3.0 * L[0].u_rec[0].pr.max; // the upper bound of the mutation rate - YC 3/6/2014
	double windowSize = log(10000)/n_loci;
		//    if (calcoptions[MUTATIONPRIORRANGE])
		//      windowsize = L[0].u_rec[0].pr.max;
		//    else
		//      windowsize = L[0].u_rec[0].win; // window width, may be used for updating - YC 3/6/2014
		//    start = 1;
	//}

	double oldScaler1 = 0.0, oldScaler2 = 0.0;
	try{
		oldScaler1 = mutationScaler_atPrev.at(IDs.at(0)); // C[ci]->G[lj].uvals[aj];
		oldScaler2 = mutationScaler_atPrev.at(IDs.at(1));
	} catch (std::exception &e) {
		std::cout << "In Chain::update_mutationScaler_Kappa(). "
				"Can't access a member of mutationScaler - vector index out of bounds\n";
	}
	double logRatio = log (oldScaler1 / oldScaler2);

	// Propose a new logRatio
	double U = runiform ();
	double newlogRatio = 0.0;
    if (U > 0.5)
    	newlogRatio = logRatio + (2.0 * U - 1.0) * windowSize;
	else
	    newlogRatio = logRatio - windowSize * U * 2.0;

	if (newlogRatio > maxRatio)
		newlogRatio = 2.0 * maxRatio - newlogRatio;
	else if (newlogRatio < -maxRatio)
		newlogRatio = 2.0 * (-maxRatio) - newlogRatio;

	double multiplier = exp ((newlogRatio - logRatio) / 2);
	std::vector<double> newScalers;
	newScalers.push_back(oldScaler1 * multiplier);
	newScalers.push_back(oldScaler2 / multiplier);

  unsigned int updateKappa = 0;
  if(lc.at(IDs.at(0)).getLikelihoodModel() == 2 && lc.at(IDs.at(1)).getLikelihoodModel() == 2) // HKY
    {
      updateKappa =1;
    }
 
  /// Propose new kappa values for the two loci
  vector<double> newKappa, newlogLik;
  newKappa.resize(2); newlogLik.resize(2);
  // Modified by YC 6/9/2014
  // Refer IMa2 codes to see how to determine the window size and the upper bound.
  double windowSize_kappa = 2.0;// 1.0;
  double upperBound_kappa = 100.0; //2.0;
  double diff_logLik = 0.0;
  for(unsigned int i=0; i<2; i++)
    {
      if(updateKappa ==1 )
	{
	  double U = runiform ();
	  if (U > 0.5)
	    {
	      newKappa.at(i) = kappa_atPrev.at(IDs.at(i)) + (2.0 * U - 1.0) * windowSize_kappa; // L[li].kappa_rec->win;
	      if (newKappa.at(i) > upperBound_kappa)// L[li].kappa_rec->pr.max)
		newKappa.at(i) = 2.0 * upperBound_kappa - newKappa.at(i);
	    }
	  else
	    {
	      newKappa.at(i) = kappa_atPrev.at(IDs.at(i)) - windowSize_kappa*U*2.0; 	 
	      if (newKappa.at(i) < 0)
		newKappa.at(i) = -newKappa.at(i);
	    }
	}
	  // Compute the new likelihood values
      newlogLik.at(i) = 0.0;//trees_atPrev.at(IDs.at(i))
      //->logLikelihood_HKY(lc.at(IDs.at(i)),newScalers.at(i),newKappa.at(i));
      diff_logLik =0.0;// += newlogLik.at(i) - logLikelihood_atPrev.at(IDs.at(i));
    }
  
  /// Metropolis-Hastings ratio
  double MHratio = exp(temperature * diff_logLik);
  

  if (runiform() < MHratio)
    {
      // Accept
      for(unsigned int i=0; i<2; i++)
	{
	  if(updateKappa ==1 ) // HKY
	    kappa_atPrev.at(IDs.at(i)) = newKappa.at(i);
	  mutationScaler_atPrev.at(IDs.at(i)) = newScalers.at(i);
	  logLikelihood_atPrev.at(IDs.at(i)) = 0.0;// newlogLik.at(i);
	}
    }

  return;
}





void Chain::updateKappa(unsigned int id_crrIter, locus lc,unsigned int save, int currentid)
{
	double U = runiform();
	/*
	double oldKappa = 0.0;
	try{
		oldKappa = kappa_atPrev.at(0);
	} catch (std::exception &e) {
			std::cout << "In Chain::updateKappa\n"
					"\tCan't access kappa - vector index out of bounds\n";
	}
	*/
	double newKappa = 0.0;
	// FIXME YC 3/11/2014
	double windowSize = 1.0; //L[li].kappa_rec->win
	double kappa_upperBound = 10; // L[li].kappa_rec->pr.max
	if(U > 0.5)
	{
		newKappa = kappa_atPrev.at(0) + (2.0 * U - 1.0) * windowSize;
		if (newKappa > kappa_upperBound)
		      newKappa = 2.0 * kappa_upperBound - newKappa;
	}
	else
	{
	    newKappa = kappa_atPrev.at(0) - windowSize * U * 2.0;
	    if (newKappa < 0)
	      newKappa = -newKappa;

	}

/*	double oldlogLik = 0.0;
	try{
		oldlogLik = logLikelihood.at(id_crrIter).at(0);
	}catch (std::exception &e) {
		std::cout << "In Chain::updateKappa\n"
				"\tCan't access logLikelihood - vector index out of bounds\n";
	}*/

	// Note that mutation rate scaler =1 on one locus case - YC 3/11/2014
	double newlogLik = trees_atPrev.at(0)->logLikelihood_HKY(lc, 1.0,newKappa);
	double MHratio = exp(temperature*(newlogLik - logLikelihood_atPrev.at(0)));

	if (runiform() < MHratio)
	{
		// accept
		kappa_atPrev.at(0) = newKappa;
		//std::vector<double> tmpKappa;
		//tmpKappa.push_back(newKappa);
		//kappa.at(id_crrIter) = tmpKappa;

		logLikelihood_atPrev.at(0) = newlogLik;
	}
	else  // reject
	{
		// kappa.at(id_crrIter) = kappa.at(id_prevIter);
	}

	return;
}

/// Added by YC 6/6/2014
void Chain::updateKappa_usePriorsOnly(unsigned int id_crrIter, locus lc,unsigned int save, int currentid) 
//(unsigned int id_crrIter, locus lc,unsigned int save, int currentid, int doisendinfo, int chainid)
{
  double U = runiform();

  double newKappa = 0.0;
	// FIXME YC 3/11/2014
  double windowSize = 2.0; // 1.0; //L[li].kappa_rec->win
  double kappa_upperBound = 100.0; // L[li].kappa_rec->pr.max
  if(U > 0.5)
    {
      newKappa = kappa_atPrev.at(0) + (2.0 * U - 1.0) * windowSize;
      if (newKappa > kappa_upperBound)
	newKappa = 2.0 * kappa_upperBound - newKappa;
    }
  else
    {
      newKappa = kappa_atPrev.at(0) - windowSize * U * 2.0;
      if (newKappa < 0)
	newKappa = -newKappa;
      
    }
  
	// Note that mutation rate scaler =1 on one locus case - YC 3/11/2014
  double newlogLik = 0.0;//trees_atPrev.at(0)->logLikelihood_HKY(lc, 1.0,newKappa);
	//double MHratio = exp(temperature*(newlogLik - logLikelihood_atPrev.at(0)));
  double MHratio = 1.0;

  if (runiform() < MHratio)
    {
      // accept
      kappa_atPrev.at(0) = newKappa;
      //std::vector<double> tmpKappa;
      //tmpKappa.push_back(newKappa);
      //kappa.at(id_crrIter) = tmpKappa;
      
      logLikelihood_atPrev.at(0) = 0.0;//newlogLik;
    }
  
  
  return;
}




/**
 * \function Chain::UpdateChain Updates the chain with newly calculated posterior, and other updates
 * \param n_current_iteration Current generation of the chain
 * \return NULL
 */
void Chain::UpdateChain(unsigned int id_crrIter, vector<locus> lc,unsigned int save, int currentid, int numprocesses, int chainid, unsigned int n_chains)
{
  
  // Update coalescent trees
  if (lociInParallel==0 && numprocesses == 1 && n_chains == 1) 
    {
      UpdateTrees_newProposal(id_crrIter, lc, save);
    } 
  else 
    {
      if(lociInParallel==1)
	{
	  //std::cout << "calling UpdateTrees_LociInMPI_newProposal\n";
	  UpdateTrees_LociInMPI_newProposal(id_crrIter, lc, save, currentid, chainid);
	  //std::cout << "Done in process = " << currentid <<"\n";
	}
      else
	{
	  try 
	    {
	      UpdateTrees_MPI_newProposal(id_crrIter, lc, save, currentid, chainid);
	    } 
	  catch (std::exception &e) 
	    {
	      std::cout << "In UpdateChain caught exception:" << e.what() << "\n";
	  }
	}
    } 

  // std::cout << "lc.at(0).get_multiLocusSpecific_mutationRate() = " << lc.at(0).get_multiLocusSpecific_mutationRate() <<"\n";
  
  if(n_loci == 1)
    {
      if(lc.at(0).getLikelihoodModel() == 2) // HKY
	{
	  try 
	    {
	      updateKappa(id_crrIter, lc.at(0), save, currentid);//, doisendinfo, chainid, coldchain);
	      //updateKappa(id_crrIter, lc.at(0), save, currentid, doisendinfo, chainid, coldchain);
	    } catch (std::exception &e) {
	    std::cout << "In UpdateKappa caught exception: " << e.what() << "\n";
	  }
	}
    }
  else if(lc.at(0).get_multiLocusSpecific_mutationRate() ==2)
    {	
      try {
	update_mutationScaler_Kappa(id_crrIter, lc,currentid, numprocesses);
      } catch (std::exception &e) {
	std::cout << "In updateMutationscaler caught exception: " << e.what() << "\n";
      }
    }
  
  //std::cout <<"End of UpdateChain() on process \n";
  

  current_iteration = id_crrIter;
  return;
}
/* End of function Chain::UpdateChain() */




/**
 * \function Chain::UpdateChain Updates the chain with newly calculated posterior, and other updates
 * \param n_current_iteration Current generation of the chain
 * \return NULL
 */
// void Chain::UpdateChain(unsigned int id_crrIter, vector<locus> lc,unsigned int save, int currentid, int coldchain)
//void Chain::UpdateChain_usePriorsOnly(unsigned int id_crrIter, vector<locus> lc,unsigned int save, int currentid, int coldchain, int numprocesses, int chainid)
void Chain::UpdateChain_usePriorsOnly(unsigned int id_crrIter, vector<locus> lc,unsigned int save, int currentid, int numprocesses, int chainid, unsigned int n_chains)
{
  //std::cout << "Updating Chains\n";
  //std::cout << "numprocesses = "<< numprocesses << " n_chains = " << n_chains <<"\n";

  // Update coalescent trees
  if (numprocesses == 1 && n_chains == 1) 
    {
      //std::cout << "Updating Trees\n";
      UpdateTrees_usePriorsOnly_newProposal(id_crrIter, lc, save);
    } 
  else 
    {
      UpdateTrees_MPI_usePriorsOnly_newProposal(id_crrIter, lc, save, currentid, chainid);
    }

  if(n_loci == 1)
    {
      if(lc.at(0).getLikelihoodModel() == 2) // HKY
	{
	  updateKappa_usePriorsOnly(id_crrIter, lc.at(0), save, currentid); //(id_crrIter, lc.at(0), save, currentid, doisendinfo, chainid);
	  
	}
    }
  else if(lc.at(0).get_multiLocusSpecific_mutationRate() ==2)
    {	
      update_mutationScaler_Kappa_usePriorsOnly(id_crrIter, lc); //(id_crrIter, lc, save, currentid, doisendinfo, chainid);
    }
  
  current_iteration = id_crrIter;
  return;
}
/* End of function Chain::UpdateChain() */




/**
 * \function Chain::SetTemperature Setting the temperature of a chain
 * \param n_temperature Temperature of the current chain
 * \return NULL
 */
void Chain::SetTemperature(double n_temperature)
{
	temperature = n_temperature;
	return;
}
/* End of function Chain::SetTemperature() */


double Chain::compute_sumOfInversePrior()
{
  unsigned int nSamples = logPrior_trees.size();
  double sumInverse = 0.0;
  for(unsigned int i=0; i<nSamples; i++)
    {
      sumInverse += 1/exp(logPrior_trees.at(i));
    }
  sumInverse /= nSamples;
  return sumInverse;
}
