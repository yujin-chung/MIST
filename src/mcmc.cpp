/* MIST copyright 2016 by Yujin Chung and Jody Hey */



#include <iostream>
#include "Chain.hpp"
#include <fstream>

using namespace std;

void MCMC::addSwapChainsTry(unsigned int chainID1, unsigned int chainID2, unsigned int successSwap)
{
  if(chainID1 > chainID2)
    {
      AccepRate_swapChains(chainID1, chainID2)++;
      if(successSwap == 1)
	AccepRate_swapChains(chainID2, chainID1)++;
    }
  else
    {
      AccepRate_swapChains(chainID2, chainID1)++;
      if(successSwap == 1)
	AccepRate_swapChains(chainID1, chainID2)++;
    }
  return;
}

void MCMC::update_orderGlobalChainIDs(unsigned int chainID1, unsigned int chainID2)
{
  unsigned int order_chainID1 = order_globalChainIDs.at(chainID1);
  order_globalChainIDs.at(chainID1) = order_globalChainIDs.at(chainID2);
  order_globalChainIDs.at(chainID2) = order_chainID1;
  return;
}

void MCMC::computeGlobalTreeAcceptanceRate()
{  
#ifdef MPI_ENABLED
  MPI::COMM_WORLD.Barrier();
	      
  unsigned int nrow = numAccpt_trees_allChains_local.rows();
  unsigned int ncol = numAccpt_trees_allChains_local.cols();
  for(unsigned int r =0; r < nrow; r++)
    {
      for(unsigned int c =0; c< ncol; c++)
	{
	  MPI::COMM_WORLD.Barrier();
	  double accpt_local = numAccpt_trees_allChains_local(r,c);
	  unsigned int tries_local = numTries_trees_allChains_local(r,c);
	  double accpt_global =0;
	  unsigned int tries_global =0;
	  MPI::COMM_WORLD.Reduce(&accpt_local, &accpt_global, 1, MPI_DOUBLE, MPI_SUM, 0);
	  MPI::COMM_WORLD.Reduce(&tries_local, &tries_global, 1, MPI_INT, MPI_SUM, 0);	  
	  MPI::COMM_WORLD.Barrier();

	  if(process_id == 0)
	    {
	      // Replace the numbers
	      numAccpt_trees_allChains_global(r,c) = accpt_global;
	      numTries_trees_allChains_global(r,c) = tries_global; 
	    }
	}
    }
#endif
  return;
}

void MCMC::compute_globalAcceptanceRate_chainSwap()
{  
#ifdef MPI_ENABLED
  MPI::COMM_WORLD.Barrier();
	      
  unsigned int nrow = AccepRate_swapChains_local.rows();
  unsigned int ncol = AccepRate_swapChains_local.cols();
  for(unsigned int r =0; r < nrow; r++)
    {
      for(unsigned int c =0; c< ncol; c++)
	{
	  MPI::COMM_WORLD.Barrier();
	  unsigned int local = AccepRate_swapChains_local(r,c);
	  unsigned int global =0;
	  MPI::COMM_WORLD.Reduce(&local, &global, 1, MPI_INT, MPI_SUM, 0); 
	  MPI::COMM_WORLD.Barrier();

	  if(process_id == 0)
	    {
	      // Replace the numbers
	      AccepRate_swapChains_global(r,c) = global;
	    }
	}
    }
#endif
  return;
}


/**
 * \function MCMC::InitializeMCMC Initialize the MCMC with chains, temperatures, etc
 * \param n_process_ID ID of the process, for temperature assignment (initial)
 * \param IM im - this is the master locus object
 * \return NULL
 */
void MCMC::InitializeMCMC(unsigned int crrProcID, IM im, unsigned int nProcs)
{
  unsigned int n_loci = im.get_nLoci();
  nProcesses = nProcs;
  n_chains = im.get_nChains();
  n_MCMCgen = im.get_nSampledTrees();
  process_id = crrProcID;
  samplingFromPriorOnly = im.get_samplingFromPriorOnly();
  nPairsMut = im.get_nPairsMut();
  // coldChainID = 0;

  lociInParallel = im.get_lociInParallel();
  if(lociInParallel == 1)
    {
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
      /*
      std::cout << "On process " << crrProcID
		<< ": reading data from locusID = " << locusID_start
		<< " to locusID = " << locusID_end
		<< " (the number of loci = " << numSubLoci << ")"
		<<"\n";
      */
    }
  else
    {
      numSubLoci = n_loci;
      locusID_start = 0;
      locusID_end = numSubLoci-1;

    }

  unsigned int totalN_chains = 0;
  if(lociInParallel ==1)
    totalN_chains = n_chains;
  else
    totalN_chains = nProcesses * n_chains;

  
  chains.resize(n_chains);
  for(unsigned int chainID = 0; chainID < n_chains; chainID++)
    {
      //	std::cout << "Initialize Chain " << chainID << "....\n";
      chains.at(chainID).set_lociInParallel(lociInParallel);
      chains.at(chainID).set_locusID_start(locusID_start);
      chains.at(chainID).set_locusID_end(locusID_end);
      chains.at(chainID).set_numSubLoci(numSubLoci);
      chains.at(chainID).set_nPairsMut(nPairsMut);

      chains.at(chainID).InitializeChain(im, n_MCMCgen, chainID, crrProcID,n_chains,nProcesses);


      // Assign temperatures of chains for Metropolis-coupled MCMC
      // "beta" is called "heat value" in MCMCMC literature,
      // but in this coding, called "temperature".
      // 
      double beta = 0.0;
      if(lociInParallel ==1 ) // loci in parallel
	beta = 1/(1+10.0 * ((double) chainID));
      else // chains in parallel
	beta = 1/(1+10.0 *(double) ((crrProcID+1)*n_chains+chainID));
	
      chains.at(chainID).SetTemperature(beta);
      chains.at(chainID).compute_rankTemp();
      chains.at(chainID).set_splittingTimeMax(im.get_splittingTimeMax());
    }

  // FIXME - YC 2/16/2015
  // This is not going to work when loci are in parallel. 
  // Need to extend it for the case
  if(crrProcID == 0 & lociInParallel ==0)
    {
      AccepRate_swapChains.resize(totalN_chains,totalN_chains);
      AccepRate_swapChains.setZero();
      for(unsigned int i=0; i< totalN_chains; i++)
	order_globalChainIDs.push_back(i);
    }

  
  // Added by YC 9/12/2014
  
  numAccpt_trees_allChains_local.resize(totalN_chains, numSubLoci);
  numTries_trees_allChains_local.resize(totalN_chains, numSubLoci);
  numAccpt_trees_allChains_local.setZero();
  numTries_trees_allChains_local.setZero();
  
  AccepRate_swapChains_local.resize(totalN_chains,totalN_chains);
  AccepRate_swapChains_local.setZero();

  local_totalLogLik.resize(1,totalN_chains);
  local_totalLogLik.setZero();
  
  if(crrProcID == 0)
    {      
      numAccpt_trees_allChains_global.resize(totalN_chains, n_loci);
      numTries_trees_allChains_global.resize(totalN_chains, n_loci);
      numAccpt_trees_allChains_global.setZero();
      numTries_trees_allChains_global.setZero();

      AccepRate_swapChains_global.resize(totalN_chains,totalN_chains);
      AccepRate_swapChains_global.setZero();
      
      global_totalLogLik.resize(1,totalN_chains);
      global_totalLogLik.setZero();
    }
  

  return;
} // END of MCMC::InitializeMCMC



/// Added by YC 6/6/2014
/**
 * \function MCMC::InitializeMCMC Initialize the MCMC with chains, temperatures, etc
 * \param n_process_ID ID of the process, for temperature assignment (initial)
 * \param IM im - this is the master locus object
 * \return NULL
 */
void MCMC::InitializeMCMC_usePriorOnly(unsigned int n_process_id, IM im, unsigned int nProcs)
{
  unsigned int n_loci = im.get_nLoci();
  nProcesses = nProcs;
  n_chains = im.get_nChains();
  n_MCMCgen = im.get_nSampledTrees();
  process_id = n_process_id;
  samplingFromPriorOnly = im.get_samplingFromPriorOnly();
  
  chains.resize(n_chains);
  for(unsigned int chainID = 0; chainID < n_chains; chainID++)
    {
      //	std::cout << "Initialize Chain " << chainID << "....\n";
      chains.at(chainID).InitializeChain_usePriorOnly(im, n_MCMCgen, chainID,n_process_id,n_chains,nProcesses);
      
      // Assign temperatures of chains for Metropolis-coupled MCMC
      // "beta" is called "temperature"
      double beta = 1/(1+10.0 * (n_process_id * n_chains + (double) chainID));
      chains.at(chainID).SetTemperature(beta);
      chains.at(chainID).compute_rankTemp();
      chains.at(chainID).set_splittingTimeMax(im.get_splittingTimeMax());
    }
  
  unsigned int totalN_chains = nProcesses * n_chains;
  if(n_process_id == 0)
    {
      // unsigned int totalN_chains = nProcesses * n_chains;
      AccepRate_swapChains.resize(totalN_chains,totalN_chains);
      AccepRate_swapChains.setZero();
      for(unsigned int i=0; i< totalN_chains; i++)
	order_globalChainIDs.push_back(i);
    }

  
  // Added by YC 9/12/2014
  
  numAccpt_trees_allChains_local.resize(totalN_chains, n_loci);
  numTries_trees_allChains_local.resize(totalN_chains, n_loci);
  numAccpt_trees_allChains_local.setZero();
  numTries_trees_allChains_local.setZero();
  
  AccepRate_swapChains_local.resize(totalN_chains,totalN_chains);
  AccepRate_swapChains_local.setZero();

  local_totalLogLik.resize(1,totalN_chains);
  local_totalLogLik.setZero();
  
  if(n_process_id == 0)
    {
      numAccpt_trees_allChains_global.resize(nProcesses*n_chains, n_loci);
      numTries_trees_allChains_global.resize(nProcesses*n_chains, n_loci);
      numAccpt_trees_allChains_global.setZero();
      numTries_trees_allChains_global.setZero();

      AccepRate_swapChains_global.resize(totalN_chains,totalN_chains);
      AccepRate_swapChains_global.setZero();
      
      global_totalLogLik.resize(1,totalN_chains);
      global_totalLogLik.setZero();
    }
  


  return;
}



void MCMC::deleteHotChains()
{
	for(unsigned int i=0; i<n_chains; i++)
	{
		if(chains.at(i).GetTemperature() != 1.0)
		{
			chains.at(i).deleteAll();
		}
		//else
		//	chains.at(i).delete_temporaryPara();
	}

}

void MCMC::deleteTrees()
{
	for(unsigned int i=0; i<n_chains; i++)
		chains.at(i).deleteTrees();
}


void MCMC::print_globalTreeAcceptanceRate()
{
  if(lociInParallel ==1)
    {
      double nIters = (double) (get_numTriesTrees_chain(0))(0,0);       
      Eigen::MatrixXd treeAccptRate =get_numAccptTrees_chain(0)*100/nIters;
      std::cout << "locusID = " << locusID_start << "-" << locusID_end <<": "
		<< treeAccptRate <<"\n";
    }
  else
    {
      double nIters = (double) numTries_trees_allChains_global(0,0); 
      Eigen::MatrixXd treeAccptRate = numAccpt_trees_allChains_global*100/nIters;
      std::cout << treeAccptRate <<"\n\n";
    }
  return;
}

void MCMC::print_AccepRate_swapChains_global()
{ 
  Eigen::MatrixXd swapChainsAcceptRate = AccepRate_swapChains_global;
  unsigned int ncols = AccepRate_swapChains_global.cols();  
  for(unsigned int i=0; i<ncols-1; i++)
    {
      for(unsigned int j=i+1; j<ncols; j++)
	{
	  double nTries = swapChainsAcceptRate(j,i);
	  if( nTries != 0.0)
	    swapChainsAcceptRate(i,j) *= 100/swapChainsAcceptRate(j,i);
	}  
    }
  std::cout << swapChainsAcceptRate <<"\n"; 
  return;
}

void MCMC::print_mcmc_initial()
{
	cout << "The number of chains (per CPU) is "<< n_chains <<".\n";
	for(unsigned int i=0; i< n_chains; i++)
	{
		cout << "Chain "<< i+1 <<": \n\n";
		chains.at(i).print_states_atIter(0);
	}
	cout << "\n";
}

void MCMC::print_savedTrees(unsigned int upto_IterID)
{
  unsigned int nloci = chains.at(0).GetNumLoci();

  std::cout << "The saved coalescent times and tree IDs\n";
  for(unsigned int iter =0; iter<=upto_IterID; iter++)
    {
      chains.at(0).print_savedStates(iter);
    }
}

void MCMC::updateTreeAcceptanceRate()
{
  for(unsigned int chainID = 0; chainID <n_chains; chainID++)
    {
      unsigned int rankTemp = getRankTempOfChain(chainID);
      Eigen::MatrixXi numTries = get_numTriesTrees_chain(chainID);
      Eigen::MatrixXd numAccpt = get_numAccptTrees_chain(chainID);
      numTries_trees_allChains_local.row(rankTemp) += numTries.row(0);
      numAccpt_trees_allChains_local.row(rankTemp) += numAccpt.row(0);
      setZero_numTriesTrees_chain(chainID);
      setZero_numAccptTrees_chain(chainID);
    }
  return;
}

/**
 * \function MCMC::UpdateAllChains() Updates all chains in the MCMC space
 * \param n_current_iteration Current generation of the MCMC
 * \return NULL
 */
void MCMC::UpdateAllChains(unsigned int id_crrIter, vector<locus> lc, unsigned int save, unsigned int swapper, unsigned int swappee)
{
  for (unsigned int  i = 0; i < n_chains; i++)
    {
      // std::cout << "Chain i = " << i <<"\n";
      try 
	{
	  if(samplingFromPriorOnly == 1)
	    {
	      chains.at(i).UpdateChain_usePriorsOnly(id_crrIter, lc,save, process_id, nProcesses, i, n_chains);  	     // (id_crrIter, lc,save, process_id, coldchain,nProcesses, i);  
	    }
	  else
	    {
	      chains.at(i).UpdateChain(id_crrIter, lc,save, process_id, nProcesses, i, n_chains);  	     
	    }	 
	}	    
      catch (std::exception &e) 
	{
	  std::cout << "In MCMC::UpdateAllChains():\n";
	  std::cout << "Can't access elements of MCMC to update - index out of bounds\n";
	}
    } // END of for (unsigned int  i = 0; i < n_chains; i++)


  if(lociInParallel==0)
    {
      // Added by YC 9/12/2014
      updateTreeAcceptanceRate();

      
      // swap chains
      swapChains(swapper, swappee);
    }
  
  return;
}     /* End of function MCMC::UpdateAllChains() */


void MCMC::collectUpdates(unsigned int savingID)
{
  if(lociInParallel ==1)
    collectUpdates_lociInParallel(savingID);
  else
    collectUpdates_chainsInParallel(savingID);

  return;
}

void MCMC::collectUpdates_lociInParallel(unsigned int savingID)
{
  if(process_id == 0)
    {
      chains.at(0).collectAllUpdates(savingID, 0);
    }

  MPI::COMM_WORLD.Barrier();
  //std::cout << "calling collectUpdates_bwProcs_lociInParallel() on procID = "<< process_id << "\n";
  collectUpdates_bwProcs_lociInParallel(savingID);
  // std::cout << "Done with collectUpdates_bwProcs_lociInParallel() on procID = "<< process_id << "\n";
  
  MPI::COMM_WORLD.Barrier();

  //std::cout << "calling collectUpdates_logPrior_lociInParallel() on procID = "<< process_id << "\n";
  collectUpdates_logPrior_lociInParallel(savingID);
  collectUpdates_logLikelihood_lociInParallel(savingID);
  //std::cout << "Done with collectUpdates_logPrior_lociInParallel() on procID = "<< process_id << "\n";

  MPI::COMM_WORLD.Barrier();
  return;
}


void MCMC::collectUpdates_logLikelihood_lociInParallel(unsigned int savingID)
{
  double global_loglik = 0.0;
  chains.at(0).compute_logLikelihood_indepLoci();
  double local_loglik = chains.at(0).getLogJointLikelihood_atPrev_lociParallel();
  MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Reduce(&local_loglik, &global_loglik, 1, MPI_DOUBLE, MPI_SUM, 0);
  MPI::COMM_WORLD.Barrier();
  if(process_id ==0)
    chains.at(0).SetLogJointLike(savingID,global_loglik);
  
  return;
}


void MCMC::collectUpdates_logPrior_lociInParallel(unsigned int savingID)
{
  double global_logPrior = 0.0;
  chains.at(0).compute_jointPrior_indepLoci();
  double local_logPrior = chains.at(0).GetLogPriorAtPrev();
  MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Reduce(&local_logPrior, &global_logPrior, 1, MPI_DOUBLE, MPI_SUM, 0);
  MPI::COMM_WORLD.Barrier();
  if(process_id ==0)
    chains.at(0).SetLogPriors(savingID,global_logPrior);
  
  return;
}

void MCMC::collectUpdates_bwProcs_lociInParallel(unsigned int savingID)
{  
  for(unsigned int p=1; p< nProcesses; p++)
    {
      MPI::COMM_WORLD.Barrier();
      unsigned int sender_locusID_start = 0;
      unsigned int sender_numSubLoci = 0;
      if(process_id ==p) // sender
	{
	  sender_locusID_start = locusID_start;
	  MPI::COMM_WORLD.Send(&sender_locusID_start, 1, MPI::INT, 0, 1497);
	 
	  sender_numSubLoci = numSubLoci;
	  MPI::COMM_WORLD.Send(&sender_numSubLoci, 1, MPI::INT, 0, 69462);
	  
	}
      else if(process_id ==0) // reciever
	{	  
	  MPI::COMM_WORLD.Recv(&sender_locusID_start, 1, MPI::INT, p, 1497);
	  MPI::COMM_WORLD.Recv(&sender_numSubLoci, 1, MPI::INT, p, 69462);
	}

      MPI::COMM_WORLD.Barrier();
      for(unsigned int i=0; i<sender_numSubLoci; i++)
	{
	  double loglik = 0.0;
	  double kappa = 0.0;
	  double mutationScaler = 0.0;
	  double logprior=0.0;
	  if(process_id ==p) // sender
	    {
	      //std::cout << "p = " << p << "sending loglik\n";
	      loglik = chains.at(0).GetLogLikelihoodAtPrev(i);
	      MPI::COMM_WORLD.Send(&loglik, 1, MPI::DOUBLE, 0, 26480);
	      //std::cout << "DONE: p = " << p << "sending loglik\n";
	      
	      //std::cout << "p = " << p << " sending kappa\n";
	      kappa = chains.at(0).getKappaAtPrev(i);
	      MPI::COMM_WORLD.Send(&kappa, 1, MPI::DOUBLE, 0,26491);
	      //std::cout << "DONE: p = " << p << " sending kappa\n";

	      logprior = chains.at(0).getLogPriorTrees_atPrev_lociParallel(i);
	      MPI::COMM_WORLD.Send(&logprior, 1, MPI::DOUBLE, 0, 28542321);

	      if( chains.at(0).GetNumLoci() > 1)
		{
		  mutationScaler = chains.at(0).getMutationScalerAtPrev(i);
		  MPI::COMM_WORLD.Send(&mutationScaler, 1, MPI::DOUBLE, 0, 149711);
		}
	  
	      node *tree2save = chains.at(0).GetTreesAtPrev(i);
	      tree2save->MPIsend_coaltimes_tipIDs();
	      tree2save->MPIsend_coaltree();
	    }
	  else if(process_id ==0)// reciever
	    {
	      MPI::COMM_WORLD.Recv(&loglik, 1, MPI::DOUBLE, p, 26480);
	      chains.at(0).SetLogLikelihood_lociMPI(savingID, sender_locusID_start+i, loglik);

	      MPI::COMM_WORLD.Recv(&logprior, 1, MPI::DOUBLE, p, 28542321);
	      chains.at(0).setEachLogPrior_lociMPI(savingID, sender_locusID_start+i, logprior);


	      MPI::COMM_WORLD.Recv(&kappa, 1, MPI::DOUBLE, p, 26491);
	      chains.at(0).setKappa(savingID, sender_locusID_start+i, kappa);
	  

	      if( chains.at(0).GetNumLoci() > 1)
		{
		  MPI::COMM_WORLD.Recv(&mutationScaler, 1, MPI::DOUBLE, p, 149711);
		  chains.at(0).setMutationScaler(savingID, sender_locusID_start+i, mutationScaler);
		}

	      // Receive the coalescent times
	      unsigned int ctsize = chains.at(0).GetTreesAtPrev(0)->size_tree()-1;
	      std::list<double> ct;
	      double ctvalue = 0.0;
	      for (unsigned int y = 0; y < ctsize; y++) 
		{
		  MPI::COMM_WORLD.Recv(&ctvalue, 1, MPI::DOUBLE, p, 298);	      
		  ct.push_back(ctvalue);
		  ctvalue = 0.0;
		}
	      chains.at(0).setCoalTimes(savingID, sender_locusID_start+i, ct);
	      
	      // Recieve the tipIDs
	      unsigned int ti = 0;
	      std::vector<unsigned int> newTipIDs;
	      for (unsigned int y = 0; y < ctsize+1; y++) 
		{
		  MPI::COMM_WORLD.Recv(&ti, 1, MPI::UNSIGNED, p, 211);
		  newTipIDs.push_back(ti);
		  ti = 0;
		}	
	      chains.at(0).setTipIDs(savingID, sender_locusID_start+i, newTipIDs);
	      
	      // Receive tree topology
	      nodeSimple *topo = new nodeSimple;
	      ct.sort();
	      topo->MPIreceive_coaltree(ct, p);
	      topo->computeSizes();
	  
	      unsigned int list_size = chains.at(0).GetListTrees().size();
	      if (list_size == 0) 
		{
		  chains.at(0).SetListTrees(topo);
		  chains.at(0).SetTreeIDs(savingID,sender_locusID_start+i, chains.at(0).GetListTrees().size()-1); // added by YC 9/25/2014
		} 
	      else 
		{
		  unsigned int found_sampleTopo = 0;
		  unsigned int count = 0;
		  while(found_sampleTopo == 0 && count < list_size) 
		    {
		      found_sampleTopo = chains.at(0).GetListTrees().at(count)->sameTopo(topo);
		      count++;
		    }
		  if (found_sampleTopo == 0) 
		    {
		      chains.at(0).SetListTrees(topo);
		      chains.at(0).SetTreeIDs(savingID, sender_locusID_start+i, chains.at(0).GetListTrees().size()-1);
		    } 
		  else 
		    {
		      chains.at(0).SetTreeIDs(savingID, sender_locusID_start+i, count-1);
		    }	
		}  
	      
	      
	    } // END of recieving 
	} // END of for(unsigned int i=0; i<numSubLoci; i++)
    }
  return;
}

void MCMC::collectUpdates_chainsInParallel(unsigned int savingID)
{  
  if(nProcesses == 1)
    {
      collectUpdates_onSameProc(savingID);
    }
  else // MPI enabled
    {
      // Determine whether the current process has the cold chain or not.
      // Also find the chain ID of cold chain on its process.
      unsigned int flag_coldChain = 0;
      unsigned int count_chains = 0;
      unsigned int coldChainID = 0;
      while(count_chains < n_chains && flag_coldChain == 0)
	{
	  if(chains.at(count_chains).getRankTemp() == 0)
	    {
	      coldChainID = count_chains;
	      flag_coldChain = 1;
	    }
	  else
	    count_chains++;
	}
      
      if(flag_coldChain == 1) // The current process has the cold chain
	{
	  if(process_id == 0) // The head node has the cold chain
	    {
	      collectUpdates_onSameProc(savingID);
	    }
	  else
	    {          
	      unsigned int procID_wColdChain = process_id;
	      MPI::COMM_WORLD.Send(&procID_wColdChain, 1, MPI::INT, 0, 32348);
	      collectUpdates_fromColdChain_onDiffProcs(savingID, procID_wColdChain, coldChainID);	      
	    }
	}
      else if(process_id == 0) // The head node does not have the cold chain
	{      
	  unsigned int procID_wColdChain = 0;
	  MPI::COMM_WORLD.Recv(&procID_wColdChain, 1, MPI::INT, MPI_ANY_SOURCE, 32348);
	  collectUpdates_fromColdChain_onDiffProcs(savingID, procID_wColdChain, coldChainID);
	}
    }
  return;
}


void MCMC::collectUpdates_fromColdChain_onDiffProcs(unsigned int savingID, unsigned int procID_wColdChain,unsigned int coldChainID)
{
  unsigned int nloci = chains.at(0).GetNumLoci();

  double logTreePrior = 0.0;
  if(process_id == procID_wColdChain)
    {
      logTreePrior = chains.at(coldChainID).GetLogPriorAtPrev();
      MPI::COMM_WORLD.Send(&logTreePrior, 1, MPI::DOUBLE, 0, 30380);

    }
  else if(process_id == 0)
    {	  
      MPI::COMM_WORLD.Recv(&logTreePrior, 1, MPI::DOUBLE, procID_wColdChain, 30380);
      chains.at(0).SetLogPriors(savingID, logTreePrior);
    }
  else
    {
      std::cout << "*** ERROR in samplingFromColdChain_onDiffProcs() ***\n";
    }
  

  for(unsigned int l=0; l< nloci; l++)
    { 
      double loglik = 0.0;
      double kappa = 0.0;
      double mutationScaler = 0.0;
      if(process_id == procID_wColdChain)
	{
	  loglik = chains.at(coldChainID).GetLogLikelihoodAtPrev(l);
	  MPI::COMM_WORLD.Send(&loglik, 1, MPI::DOUBLE, 0, 32376);

	  kappa = chains.at(coldChainID).getKappaAtPrev(l);
	  MPI::COMM_WORLD.Send(&kappa, 1, MPI::DOUBLE, 0, 31400);
	  
	  if(nloci > 1)
	    {
	      mutationScaler = chains.at(coldChainID).getMutationScalerAtPrev(l);
	      MPI::COMM_WORLD.Send(&mutationScaler, 1, MPI::DOUBLE, 0, 314051303);
	    }
	  
	  node *tree2save = chains.at(coldChainID).GetTreesAtPrev(l);
	  tree2save->MPIsend_coaltimes_tipIDs();
	  tree2save->MPIsend_coaltree();
	  
	}
      else if(process_id == 0)
	{
	  chains.at(0).resizeLogLikelihood(savingID);
	  chains.at(0).resizeTreeIDs(savingID);
	  chains.at(0).resizeCoalTimes(savingID);
	  chains.at(0).resizeTipIDs(savingID);
	  chains.at(0).resizeKappa(savingID);
	  chains.at(0).resizeMutationScaler(savingID);
  
  
	  MPI::COMM_WORLD.Recv(&loglik, 1, MPI::DOUBLE, procID_wColdChain, 32376);
	  chains.at(0).SetLogLikelihood(savingID, l, loglik);
	  
	  MPI::COMM_WORLD.Recv(&kappa, 1, MPI::DOUBLE, procID_wColdChain, 31400);
	  chains.at(0).setKappa(savingID, l, kappa);
	  
	  if(nloci > 1)
	    {
	      MPI::COMM_WORLD.Recv(&mutationScaler, 1, MPI::DOUBLE, procID_wColdChain, 314051303);
	      chains.at(0).setMutationScaler(savingID, l, mutationScaler);
	    }

	  // Receive the coalescent times
	  unsigned int ctsize = chains.at(0).GetTreesAtPrev(l)->size_tree()-1;
	  std::list<double> ct;
	  double ctvalue = 0.0;
	  for (unsigned int y = 0; y < ctsize; y++) {
	    MPI::COMM_WORLD.Recv(&ctvalue, 1, MPI::DOUBLE, procID_wColdChain, 298);					
	    ct.push_back(ctvalue);
	    // coalTimes.at(id_crrIter).at(i).push_back(ctvalue);
	    ctvalue = 0.0;
	  }
	  chains.at(0).setCoalTimes(savingID, l, ct);

	  // Receive the and tipIDs
	  unsigned int ti = 0;
	  std::vector<unsigned int> newTipIDs;
	  for (unsigned int y = 0; y < ctsize+1; y++) 
	    {
	      MPI::COMM_WORLD.Recv(&ti, 1, MPI::UNSIGNED, procID_wColdChain, 211);
	      newTipIDs.push_back(ti);
		//tipIDs.at(id_crrIter).at(i).push_back(ti);
	      ti = 0;
	    }	
	  chains.at(0).setTipIDs(savingID, l, newTipIDs);

	  // Receive tree topology
	  nodeSimple *topo = new nodeSimple;
	  ct.sort();
	  topo->MPIreceive_coaltree(ct, procID_wColdChain);
	  topo->computeSizes();
	  
	  unsigned int list_size = chains.at(0).GetListTrees().size();
	  if (list_size == 0) 
	    {
	      chains.at(0).SetListTrees(topo);
	      chains.at(0).SetTreeIDs(savingID, l, chains.at(0).GetListTrees().size()-1); // added by YC 9/25/2014
	    } 
	  else 
	    {
	      unsigned int found_sampleTopo = 0;
	      unsigned int count = 0;
	      while(found_sampleTopo == 0 && count < list_size) 
		{
		  found_sampleTopo = chains.at(0).GetListTrees().at(count)->sameTopo(topo);
		  count++;
		}
	      if (found_sampleTopo == 0) 
		{
		  chains.at(0).SetListTrees(topo);
		  chains.at(0).SetTreeIDs(savingID, l, chains.at(0).GetListTrees().size()-1);
		} 
	      else 
		{
		  chains.at(0).SetTreeIDs(savingID, l, count-1);
		}	
	    }  
	} // END of else if(process_id == 0)
      else
	{
	  std::cout << "*** ERROR in samplingFromColdChain_onDiffProcs() ***\n";
	}
    }
  
  return;
}


void MCMC::collectUpdates_onSameProc(unsigned int savingID)
{
  // std::cout << "In MCMC::collectUpdates_onSameProc(): n_chains = " << n_chains <<"\n";

  // Find the chain ID with the cold chain
  unsigned int coldChainID = 0;
  if(n_chains > 1)
    {
      unsigned int id = 0;
      unsigned int foundColdChain = 0;
      while(id < n_chains && foundColdChain == 0)
	{
	  // std::cout << "chains.at(id).getRankTemp() = " << chains.at(id).getRankTemp() <<"\n";
	  if(chains.at(id).getRankTemp() == 0) // the cold chain has rank 0.
	    {
	      coldChainID = id;
	      foundColdChain =1;
	    }
	  id++;
	}
      if(foundColdChain == 0)
	{
	  std::cout << "*** Error in MCMC::collectUpdates() ***\n";
	  std::cout << "*** Can't find a cold chain ***\n";
	}
    }

  if(coldChainID == 0)
    {
      chains.at(0).collectAllUpdates(savingID, 0);
    }
  else
    {
      collectUpdates_fromColdChain_onSameProc(savingID, coldChainID);
    }

  return;
}



void MCMC::collectUpdates_fromColdChain_onSameProc(unsigned int savingID, unsigned int coldChainID)
{
  /*
  std::cout << "In MCMC::collectUpdates_fromColdChain_onSameProc() savingID = " << savingID
	    << " coldChainID  = " << coldChainID 
	    <<"\n";
  */

  unsigned int nloci = chains.at(coldChainID).GetNumLoci();

  chains.at(0).SetLogPriors(savingID, chains.at(coldChainID).GetLogPriorAtPrev());
  chains.at(0).SetKappa(savingID, chains.at(coldChainID).GetKappaAtPrev());
  chains.at(0).SetLogLikelihood(savingID, chains.at(coldChainID).getLogLikelihoodAtPrev());

  if (nloci > 1) 
    chains.at(0).SetMutationScaler(savingID, chains.at(coldChainID).GetMutationScalerAtPrev());

  chains.at(0).resizeTreeIDs(savingID);

  for (unsigned int i = 0; i < nloci; i++) 
    {
      node *tree2save = chains.at(coldChainID).GetTreesAtPrev(i);
      chains.at(0).compute_coalTimes_tipIDs(savingID, i, tree2save);
  
      unsigned int list_size = chains.at(0).GetListTrees().size();
  
      if (list_size == 0) 
	{
	  nodeSimple *topo = new nodeSimple;
	  std::list<double> ct = chains.at(0).GetCoalTimes(savingID, i);
	  ct.sort();
	  topo->convert(chains.at(coldChainID).GetTreesAtPrev(i), ct, tree2save->size_tree());
	  topo->computeSizes();
	  chains.at(0).SetListTrees(topo);
	  chains.at(0).SetTreeIDs(savingID, i, chains.at(0).GetListTrees().size()-1); // added by YC 9/25/2014
	} 
      else 
	{
	  unsigned int found_sampleTopo = 0;
	  unsigned int count = 0;
	  while(found_sampleTopo == 0 && count < list_size) 
	    {
	      found_sampleTopo = chains.at(0).GetListTrees().at(count)->sameTopo(tree2save);
	      count++;
	    }
	  if (found_sampleTopo == 0) 
	    {
	      nodeSimple *topo = new nodeSimple;
	      std::list<double> ct = chains.at(0).GetCoalTimes(savingID, i);
	      ct.sort();
	      topo->convert(chains.at(coldChainID).GetTreesAtPrev(i), ct, tree2save->size_tree());
	      topo->computeSizes();
	      chains.at(0).SetListTrees(topo);
	      chains.at(0).SetTreeIDs(savingID, i, chains.at(0).GetListTrees().size()-1);
	    } 
	  else 
	    {
	      chains.at(0).SetTreeIDs(savingID, i, count-1);
	    }	
	}
    }	  
  return;
}





void MCMC::swapChains(unsigned int swapper_globalChainID, unsigned int swappee_globalChainID)
{

  unsigned int procID_swapper = 0;
  unsigned int procID_swappee = 0;
  unsigned int swapperID_onProc = swapper_globalChainID;
  unsigned int swappeeID_onProc = swappee_globalChainID;
  if(nProcesses > 1)
    {
      procID_swapper = static_cast<int>(swapper_globalChainID/n_chains);
      procID_swappee = static_cast<int>(swappee_globalChainID/n_chains);
      swapperID_onProc = swapper_globalChainID - n_chains * procID_swapper;
      swappeeID_onProc = swappee_globalChainID - n_chains * procID_swappee;
    }

  if (nProcesses > 1 || n_chains>1) 
    {
      if(process_id == procID_swapper || process_id == procID_swappee)
	{
	  if (procID_swapper == procID_swappee) 
	    {
	      SwapChainHeats_WithinProcess(swapperID_onProc, swappeeID_onProc);
	    } 
	  else
	    {
	      SwapChainHeats_BwProcesses(swapper_globalChainID, swappee_globalChainID, n_chains, process_id);
	    }	      
	}
    }
  return;
}


/*
void MCMC::samplingFromColdChain(unsigned int id_crrIter)
{
  unsigned int flag_coldChain = 0;
  unsigned int count_chains = 0;
  unsigned int coldChainID = 0;
  while(count_chains < n_chains && flag_coldChain == 0)
    {
      if(chains.at(count_chains).getRankTemp() == 0)
	{
	  coldChainID = count_chains;
	  flag_coldChain = 1;
	}
      else
	count_chains++;
    }
  
  if(flag_coldChain == 1) // cold chain
    {
      if(process_id == 0) // The head node has the cold chain
	{
	  samplingFromColdChain_onTheSameProc(id_crrIter, coldChainID);
	}
      else
	{
	  unsigned int procID_wColdChain = process_id;      
	  MPI::COMM_WORLD.Send(&procID_wColdChain, 1, MPI::INT, 0, 32348);
	  samplingFromColdChain_onDiffProcs(id_crrIter, procID_wColdChain, coldChainID);
	}
    }
  else if(process_id == 0) // The head node does not have the cold chain
    {
      unsigned int procID_wColdChain = 0;
      MPI::COMM_WORLD.Recv(&procID_wColdChain, 1, MPI::INT, MPI_ANY_SOURCE, 32348);
      samplingFromColdChain_onDiffProcs(id_crrIter, procID_wColdChain, coldChainID);
    }

  return;
}

void MCMC::samplingFromColdChain_onDiffProcs(unsigned int id_crrIter, unsigned int procID_wColdChain,unsigned int coldChainID)
{
  unsigned int nloci = chains.at(0).GetNumLoci();

  double logTreePrior = 0.0;
  if(process_id == procID_wColdChain)
    {
      logTreePrior = chains.at(coldChainID).GetLogPriorAtPrev();
      MPI::COMM_WORLD.Send(&logTreePrior, 1, MPI::DOUBLE, 0, 30380);

    }
  else if(process_id == 0)
    {	  
      MPI::COMM_WORLD.Recv(&logTreePrior, 1, MPI::DOUBLE, procID_wColdChain, 30380);
      chains.at(0).SetLogPriors(id_crrIter, logTreePrior);
    }
  else
    {
      std::cout << "*** ERROR in samplingFromColdChain_onDiffProcs() ***\n";
    }
  

  for(unsigned int l=0; l< nloci; l++)
    { 
      double loglik = 0.0;
      double kappa = 0.0;
      double mutationScaler = 0.0;
      if(process_id == procID_wColdChain)
	{
	  loglik = chains.at(coldChainID).GetLogLikelihoodAtPrev(l);
	  MPI::COMM_WORLD.Send(&loglik, 1, MPI::DOUBLE, 0, 32376);

	  kappa = chains.at(coldChainID).getKappaAtPrev(l);
	  MPI::COMM_WORLD.Send(&kappa, 1, MPI::DOUBLE, 0, 31400);
	  
	  if(nloci > 1)
	    {
	      mutationScaler = chains.at(coldChainID).getMutationScalerAtPrev(l);
	      MPI::COMM_WORLD.Send(&mutationScaler, 1, MPI::DOUBLE, 0, 314051303);
	    }
	  
	  node *tree2save = chains.at(coldChainID).GetTreesAtPrev(l);
	  tree2save->MPIsend_coaltimes_tipIDs();
	  tree2save->MPIsend_coaltree();
	  
	}
      else if(process_id == 0)
	{
	  MPI::COMM_WORLD.Recv(&loglik, 1, MPI::DOUBLE, procID_wColdChain, 32376);
	  chains.at(0).SetLogLikelihood(id_crrIter, l, loglik);
	  
	  MPI::COMM_WORLD.Recv(&kappa, 1, MPI::DOUBLE, procID_wColdChain, 31400);
	  chains.at(0).setKappa(id_crrIter, l, kappa);
	  
	  if(nloci > 1)
	    {
	      MPI::COMM_WORLD.Recv(&mutationScaler, 1, MPI::DOUBLE, procID_wColdChain, 314051303);
	      chains.at(0).setMutationScaler(id_crrIter, l, mutationScaler);
	    }

	  // Receive the coalescent times
	  unsigned int ctsize = chains.at(0).GetTreesAtPrev(l)->size_tree()-1;
	  std::list<double> ct;
	  double ctvalue = 0.0;
	  for (unsigned int y = 0; y < ctsize; y++) {
	    MPI::COMM_WORLD.Recv(&ctvalue, 1, MPI::DOUBLE, procID_wColdChain, 298);					
	    ct.push_back(ctvalue);
	    // coalTimes.at(id_crrIter).at(i).push_back(ctvalue);
	    ctvalue = 0.0;
	  }
	  chains.at(0).setCoalTimes(id_crrIter, l, ct);

	  // Receive the and tipIDs
	  unsigned int ti = 0;
	  std::vector<unsigned int> newTipIDs;
	  for (unsigned int y = 0; y < ctsize+1; y++) 
	    {
	      MPI::COMM_WORLD.Recv(&ti, 1, MPI::UNSIGNED, procID_wColdChain, 211);
	      newTipIDs.push_back(ti);
		//tipIDs.at(id_crrIter).at(i).push_back(ti);
	      ti = 0;
	    }	
	  chains.at(0).setTipIDs(id_crrIter, l, newTipIDs);

	  // Receive tree topology
	  nodeSimple *topo = new nodeSimple;
	  ct.sort();
	  topo->MPIreceive_coaltree(ct, procID_wColdChain);
	  topo->computeSizes();
	  
	  unsigned int list_size = chains.at(0).GetListTrees().size();
	  if (list_size == 0) 
	    {
	      chains.at(0).SetListTrees(topo);
	      chains.at(0).SetTreeIDs(id_crrIter, l, chains.at(0).GetListTrees().size()-1); // added by YC 9/25/2014
	    } 
	  else 
	    {
	      unsigned int found_sampleTopo = 0;
	      unsigned int count = 0;
	      while(found_sampleTopo == 0 && count < list_size) 
		{
		  found_sampleTopo = chains.at(0).GetListTrees().at(count)->sameTopo(topo);
		  count++;
		}
	      if (found_sampleTopo == 0) 
		{
		  chains.at(0).SetListTrees(topo);
		  chains.at(0).SetTreeIDs(id_crrIter, l, chains.at(0).GetListTrees().size()-1);
		} 
	      else 
		{
		  chains.at(0).SetTreeIDs(id_crrIter, l, count-1);
		}	
	    }  
	} // END of else if(process_id == 0)
      else
	{
	  std::cout << "*** ERROR in samplingFromColdChain_onDiffProcs() ***\n";
	}
    }
  return;
}

void MCMC::samplingFromColdChain_onTheSameProc(unsigned int id_crrIter,unsigned int coldChainID)
{
  unsigned int nloci = chains.at(coldChainID).GetNumLoci();

  chains.at(0).SetLogPriors(id_crrIter, chains.at(coldChainID).GetLogPriorAtPrev());
  chains.at(0).SetKappa(id_crrIter, chains.at(coldChainID).GetKappaAtPrev());

  if (nloci > 1) 
    {
      chains.at(0).SetMutationScaler(id_crrIter, chains.at(coldChainID).GetMutationScalerAtPrev());
    }

  for (unsigned int i = 0; i < nloci; i++) 
    {
      node *tree2save = chains.at(coldChainID).GetTreesAtPrev(i);
      chains.at(0).compute_coalTimes_tipIDs(id_crrIter, i, tree2save);
      chains.at(0).SetLogLikelihood(id_crrIter, i, chains.at(coldChainID).GetLogLikelihoodAtPrev(i));
      unsigned int list_size = chains.at(0).GetListTrees().size();		
      if (list_size == 0) 
	{
	  nodeSimple *topo = new nodeSimple;
	  std::list<double> ct = chains.at(0).GetCoalTimes(id_crrIter, i);
	  ct.sort();
	  topo->convert(chains.at(coldChainID).GetTreesAtPrev(i), ct, tree2save->size_tree());
	  topo->computeSizes();
	  chains.at(0).SetListTrees(topo);
	  chains.at(0).SetTreeIDs(id_crrIter, i, chains.at(0).GetListTrees().size()-1); // added by YC 9/25/2014
	} 
      else 
	{
	  unsigned int found_sampleTopo = 0;
	  unsigned int count = 0;
	  while(found_sampleTopo == 0 && count < list_size) 
	    {
	      found_sampleTopo = chains.at(0).GetListTrees().at(count)->sameTopo(tree2save);
	      count++;
	    }
	  if (found_sampleTopo == 0) 
	    {
	      nodeSimple *topo = new nodeSimple;
	      std::list<double> ct = chains.at(0).GetCoalTimes(id_crrIter, i);
	      ct.sort();
	      topo->convert(chains.at(coldChainID).GetTreesAtPrev(i), ct, tree2save->size_tree());
	      topo->computeSizes();
	      chains.at(0).SetListTrees(topo);
	      chains.at(0).SetTreeIDs(id_crrIter, i, chains.at(0).GetListTrees().size()-1);
	    } 
	  else 
	    {
	      chains.at(0).SetTreeIDs(id_crrIter, i, count-1);
	    }	
	}
    }	  
  return;
}
*/


/**
 * \function MCMC::PrintMCMCState() Returns/prints the current MCMC state
 * \param NULL
 * \return NULL
 */

void MCMC::PrintMCMCState()
{
	unsigned int i = 0;
	for( i = 0; i < chains.size(); i++) {
		try {
			chains.at(i).PrintChain();
		} catch (std::exception &e) {
			std::cout << "Can't access elements of MCMC to print - index out of bounds\n";
		}
	}
	return;
}
/* End of function MCMC::PrintMCMCState() */

//AS: this function has to do all that the sending from chain identified by chainid
//to chain 0 on the head node

void MCMC::sharecoaltree(unsigned int chainid, unsigned int id_crrIter)
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
	chains.at(0).SetLogPriors(id_crrIter, chains.at(chainid).GetLogPriorAtPrev());
	chains.at(0).SetKappa(id_crrIter, chains.at(chainid).GetKappaAtPrev());
	if (nloci > 1) {
		chains.at(0).SetMutationScaler(id_crrIter, chains.at(chainid).GetMutationScalerAtPrev());
	}
	
	return;
}



/**
 * \function MCMC::GetChainTemperature Gets the temperature of the chain identified by chain_ID
 * \param chain_ID Identifier of the chain
 * \return double Temperature (beta)
 */

double MCMC::GetChainTemperature(unsigned int chain_ID)
{
  double temp = 0.0;
  try 
    {
      temp = chains.at(chain_ID).GetTemperature();
    }
  catch (std::exception &e) 
    {
      std::cout << "*** Error in MCMC::GetChainTemperature() ***\n";
      std::cout << "\ttried to get the temperature of chain with ID = "<< chain_ID << " on process " << process_id << "\n";
    }	
  return temp;
}
/* End of function MCMC::GetChainTemperature() */
/**
 * \function MCMC::GetChainPosterior Returns the posterior density of the current chain
 * identified by chain_ID
 * \param chain_ID Identifier of the chain
 * \return double posterior density
 */
double MCMC::GetChainPosterior(int chain_ID)
{
	double post = 0.0;
	try {
		post = chains.at(chain_ID).GetPosterior();
	} catch (std::exception &e) {
		std::cout << "Can't access posterior of this chain!! Index is out of bounds\n";
	}
	return post;
}
/* End of function MCMC::GetChainPosterior() */


Chain MCMC::getColdChain()
{
	unsigned int id = 0;
	// Cold chain's temperature (beta) is 1.
	while(GetChainTemperature(id)!= 1.0 && id < n_chains)
	{
		// REMOVE
		//std::cout << "id = " <<id<< " temperature = " << GetChainTemperature(id) <<"\n";

		id++;
	}
	return chains.at(id);
}

int MCMC::getColdChain(int process_id)
{
	unsigned int id = 0;
	// Cold chain's temperature (beta) is 1.
	for (unsigned int i = 0; i < n_chains; i++)
	{
		// REMOVE
		//std::cout << "id = " <<id<< " temperature = " << GetChainTemperature(id) << " on process " << process_id << "\n";
		if (GetChainTemperature(i) == 1.0) {
			return i;
		}

	}
	//std::cout << "No cold chains on this process!";
	return -1;
}



/**
 * \function MCMC::fileprintAllChains() printing the state of all chains at the given iteration,
 * \param id_iter Id of the iteration in MCMC
 * \return NULL
 */

void MCMC::fileprintAllChains(unsigned int id_iter, std::ofstream& fp)
{
	for(unsigned int i=0; i<n_chains; i++)
	{
		try {
			chains.at(i).fileprint_states_atIter(id_iter, fp);
		} catch (std::exception &e) {
			std::cout << "Can't access elements of this chain to print - index out of bounds\n";
		}
	}
}
/* End of function MCMC::fileprintAllChains() */



/**
 * \function MCMC::printAllChains() printing the state of all chains at the given iteration,
 * \param id_iter Id of the iteration in MCMC
 * \return NULL
 */

void MCMC::printAllChains(unsigned int id_iter)
{
	for(unsigned int i=0; i<n_chains; i++)
	{
		try {
			chains.at(i).print_states_atIter(id_iter);
		} catch (std::exception &e) {
			std::cout << "Can't access elements of this chain to print - index out of bounds\n";
		}
	}
}
/* End of function MCMC::printAllChains() */


void MCMC::print_coldChainState()
{
  if(lociInParallel ==1)
  {
	chains.at(0).print_states_atIter(0);
  }
  else
  {
    for(unsigned int c = 0; c<n_chains; c++)
      {
	if(chains.at(c).getRankTemp() == 0)
	  chains.at(c).print_states_atIter(0);  	
      }
  }	
  return;
}

void MCMC::printMCMCstate_eachChain(unsigned int chainID, unsigned int crrIter)
{
  chains.at(chainID).print_states_atIter(crrIter);
  return;
}

void MCMC::print_globalMCMCstate()
{
  unsigned int nrow = global_mcmcStates.rows();
  for(unsigned int r=0; r < nrow; r++)
    {
      std::cout << global_mcmcStates.row(r) << " " << global_trees.at(r) <<"\n";
      // std::cout << global_mcmcStates.at(i) <<"\n";
    }
  return;
}

void MCMC::print_allLogLik()
{
  if(lociInParallel==1)
    {
      //std::cout << "printing alllogLik\n";
      //std::cout << "n_chains = " << n_chains <<"\n";
      unsigned int chID =0;
      // for(unsigned int chID = 0; chID <n_chains; chID++)
      {
	double loglik = chains.at(chID).compute_totalLogLik();
	// unsigned int order = chains.at(chID).getRankTemp();
	// std::cout << "loglik = " << loglik  << " and order = "<< order <<"\n";
	if(std::isnan(loglik) ==0.0)
	  local_totalLogLik(0, 0) = loglik;
      }
      // std::cout << "local_totalLogLik = " << local_totalLogLik <<"\n";
      MPI::COMM_WORLD.Barrier(); 
      
      for(unsigned int id = 0; id< n_chains; id++)
	{    
	  MPI::COMM_WORLD.Barrier(); 
	  double local = local_totalLogLik(0,id);
	  double global =0;
	  MPI::COMM_WORLD.Reduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, 0);
	  if(process_id ==0)
	    {
	      global_totalLogLik(0,id) = global;
	      // std::cout << "global = " << global <<"\n";
	    }
	  MPI::COMM_WORLD.Barrier();
	}
    }
  else
    {
      for(unsigned int chID = 0; chID <n_chains; chID++)
	{
	  double loglik = chains.at(chID).compute_totalLogLik();
	  unsigned int order = chains.at(chID).getRankTemp();
	  local_totalLogLik(0, order) = loglik;
	}	 
      MPI::COMM_WORLD.Barrier(); 
      
      for(unsigned int id = 0; id< nProcesses*n_chains; id++)
	{    
	  MPI::COMM_WORLD.Barrier(); 
	  double local = local_totalLogLik(0,id);
	  double global =0;
	  MPI::COMM_WORLD.Reduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, 0);
	  if(process_id ==0)
	    global_totalLogLik(0,id) = global;
	  MPI::COMM_WORLD.Barrier();
	}
    }
      
  MPI::COMM_WORLD.Barrier(); 
  if(process_id == 0)
    {
      std::cout << global_totalLogLik <<"\n";
    }
  local_totalLogLik.setZero();
  return;
}

void MCMC::collect_localCurrentMCMCstates()
{
  // Relace 'global_mcmcStates' by the sum of local_mcmcStates
  unsigned int nrow = local_mcmcStates.rows();
  unsigned int ncol = local_mcmcStates.cols();
  for(unsigned int r = 0; r<nrow; r++)
    {
      for(unsigned int c =0; c<ncol; c++)
	{	 
	  MPI::COMM_WORLD.Barrier(); 
	  double local = local_mcmcStates(r,c);
	  double global =0;
	  MPI::COMM_WORLD.Reduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, 0);
	  global_mcmcStates(r,c) = global;
	  MPI::COMM_WORLD.Barrier();
	}
    }
  /*
  // Broadcasting the tree states
  for(unsigned int procID =0; procID <nProcesses; procID++)
    {
#ifdef MPI_ENABLED
      MPI::COMM_WORLD.Barrier();
#endif
      if(procID == process_id)
	{
	  unsigned int 


	  if(process_id == 0)
	    {
	      std::cout << "in collect_localCurrentMCMCstates()\n";
	      for(unsigned int chID =0; chID<n_chains; chID++)
		{
		  std::string state = chains.at(chID).print2string_states_atCrr();
		  std::cout << state;
		}
	    }
	}
#ifdef MPI_ENABLED
      MPI::COMM_WORLD.Barrier();
#endif
    }
  */
  return;
}


/**
 * \function MCMC::ChangeMCMCState() Updating the MCMC state
 * \param chain_ID Identifier of which chain to update
 * \param new_temperature Temperature of the chain after update/swap
 * \return NULL
 */

void MCMC::ChangeMCMCState(int chain_ID, double new_temperature)
{
	int pid = GetProcessID();
	try {
		chains.at(chain_ID).SetTemperature(new_temperature);
	} catch (std::exception &e) {
		std::cout << "Can't set new temperature in chain " << chain_ID << " on process " << pid << "\n";
	}
}

/* End of function MCMC::ChangeMCMCState() */
/**
 * \function MCMC::SwapChainHeats Perform swap of heats between two chains
 * identified by chain1 and chain2
 * \param chain1 identifier of chain1
 * \param chain2 identifier of chain2
 * \return NULL
 */

double MCMC::GetMHSwapTerm(int chain_id)
{
  double mhswapterm = 0.0;
  for (int l = 0; l < chains.at(chain_id).GetNumLoci(); l++) 
    {
      mhswapterm = mhswapterm + chains.at(chain_id).GetLogLikelihoodAtPrev(l);
      // By YC - 12/9/2014
      // Since parameters for mutation/substitution model are from uniform, 
      // their prior distribution may not be in the MH term.
      // mhswapterm = mhswapterm + log(1 / GetUpperKappa());
      
    }
  mhswapterm = mhswapterm + chains.at(chain_id).GetLogPriorAtPrev();
  return mhswapterm;
} 


void MCMC::compute_localAcceptanceRate_chainSwap(unsigned int rankTemp1, unsigned int rankTemp2, unsigned int swapSuccess)
{
  if(rankTemp1 > rankTemp2)
    {
      AccepRate_swapChains_local(rankTemp1,rankTemp2)++;
      if(swapSuccess == 1)
	AccepRate_swapChains_local(rankTemp2,rankTemp1)++;
    }
  else if(rankTemp1 < rankTemp2)
    {
      AccepRate_swapChains_local(rankTemp2,rankTemp1)++;
      if(swapSuccess == 1)
	AccepRate_swapChains_local(rankTemp1,rankTemp2)++;
    }
  else
    {
      std::cout << "***** Error in MCMC::compute_localAcceptanceRate_chainSwap() *****\n";
      std::cout << "\tTwo arguments rankTemp1 and rankTemp2 are the same as " << rankTemp1 <<"but they should be different.\n";
    }
  return;
}

int MCMC::SwapChainHeats_WithinProcess(int chain1, int chain2)
{
  if(chain1>=n_chains || chain2 >=n_chains)
    {
      std::cout << "***** Error in SwapChainHeats_WithinProcess() *****\n";
      std::cout << "\tchain IDs (chain1 = " << chain1 << " or chain2 = " << chain2 <<") are larger or equal to n_chains = " << n_chains << " on process "<< process_id << "\n";
    }

	int swapvar = 0;
		
	double temp1 = 0.0;
	double temp2 = 0.0;
	try {
		temp1 = GetChainTemperature(chain1);
		//temp1 = chains.at(chain1).GetTemperature();
	} catch (std::exception &e) {
		std::cout << "Can't access temperature of " << chain1 << " on process \n";
		return swapvar;
	}
	try {
		temp2 = GetChainTemperature(chain2);
		//	temp2 = chains.at(chain2).GetTemperature();
	} catch (std::exception &e) {
		std::cout << "Can't access temperature of " << chain2 << " on process \n";
		return swapvar;
	}
	//AS: compute swap probability - MH ratio? let's call this mhterm
	// double mhterm = exp(GetMHSwapTerm(chain1) * temp1 - GetMHSwapTerm(chain2) * temp2);
	// Modified by Yujin - 12/9/2014
	double mhterm = exp((temp1 - temp2) *(GetMHSwapTerm(chain2) - GetMHSwapTerm(chain1)));
	// double mhterm = (GetMHSwapTerm(chain1) * temp1)/(GetMHSwapTerm(chain2) * temp2);
	
	if (mhterm >= 1.0 || mhterm > runiform()) 
	  {
	    // Swap temperatures
	    ChangeMCMCState(chain1, temp2);
	    ChangeMCMCState(chain2, temp1);

	    // Recompute the ranks of temperatures
	    chains.at(chain1).compute_rankTemp();
	    chains.at(chain2).compute_rankTemp();
	    
	    // return variable
	    swapvar = 1;
	} //AS: within process swap completed

	// Compute acceptance rates for chain swap
	unsigned int rankTemp1 = chains.at(chain1).getRankTemp();
	unsigned int rankTemp2 = chains.at(chain2).getRankTemp();
	compute_localAcceptanceRate_chainSwap(rankTemp1,rankTemp2,swapvar);

	return swapvar;
}	/* End of function MCMC::SwapChainHeats() */





// Yujin's new version - updated Arun's version --- YC 8/13/2014
/**
 * \function MCMC:SwapChainHeats_BwProcesses Performs swapping of heats between two processes
 * identified by swapA and swapB. This function is essentially called only if MPI is enabled
 * \param n_chains - number of chains on each process
 * \param currentid - current process ID
 * \return int swapvar - indicator of swap being successful or not
 */
int MCMC::SwapChainHeats_BwProcesses(int swapper, int swappee, unsigned int n_chains, int currentid)
{  
  if (swapper == swappee) 
    {
      std::cout << "*** Error *** in MCMC::SwapChainHeats_BwProcesses()\n";
      std::cout << "\t Swapper and Swappee assignment might have bugs. They cannot be the same!\n";
      std::cout << "\t swapper = "<< swapper <<" swappee = "<< swappee <<"\n";
    }
  
  //AS: first need to figure out if swapper and swappee are on the same or diff processes
  int swapAflag = 0;
  int swapBflag = 0;
  int whichElementA = 0;
  int whichElementB = 0;
  // Relaced "floor" by "int" to avoid rounding error -- YC 8/13/2014
  int procIdForA = static_cast<int>(swapper/n_chains);
  int procIdForB = static_cast<int>(swappee/n_chains);

  if (procIdForA == currentid) {
    swapAflag = 1;
  }
  if (procIdForB == currentid) {
    swapBflag = 1;
  }
  whichElementA = swapper - procIdForA * n_chains;
  whichElementB = swappee - procIdForB * n_chains;


  if (whichElementA == whichElementB && procIdForA == procIdForB) 
    {
      std::cout << "*** Error *** in MCMC::SwapChainHeats_BwProcesses()\n";
      std::cout << "Swapper and swappee are: " << swapper << "\t" << swappee << "\n";
      std::cout << "Their process IDs and chain IDs (on the process) are the same and this cannot be happening!\n";
      return 0;
    }

  //AS: if both flags are 0, then the swap attempt doesn't involve current process, so continue
  if (swapAflag == 0 && swapBflag == 0) 
    {
      std::cout << "*** Error *** in MCMC::SwapChainHeats_BwProcesses()\n";
      std::cout << "Processor " <<currentid << " is not swapping in step , so continuing on to next generation \n";
      return 0;
    }	

  
  
  // YC simplied the below. 8/13/2014
  // if both swapflags are 1, then do swapping within the process.
  // if one is 1 and the other is 0, then do swapping between processes.
  int swapvar = 0;
  double abeta, bbeta = 0.0;
  if (swapAflag == 1 && swapBflag == 1) 
    {      
      //AS: if both flags are 1, then the swapping chains are on the same process - so go into within process swap
      //AS: within process swap between swapper[step] and swappee[step]
      std::cout << "Processor " << currentid << " is going into within procesor swap in step \n";
      return SwapChainHeats_WithinProcess(whichElementA, whichElementB);
    }
  else if(swapAflag == 1)
    {
      // Yujin modified this part -- 9/4/2014
      // On the process with swapper, it receive the temperature and MHterm from swappee
      // and compute the MH ratio and determine if it rejects or not.
      MPI::Status status;
      double mhtermA = 0.0;
      double mhtermB = 0.0;

      // Get the rank of temperature of swapper
      unsigned int rankTemp1 = chains.at(whichElementA).getRankTemp();
      unsigned int rankTemp2 = 0;

      // Obtain the MHterm and temperature -- YC 9/4/2014
      mhtermA = GetMHSwapTerm(whichElementA);
      abeta = chains.at(whichElementA).GetTemperature();
      try 
	{
	  // Receive the MHterm and temperature of swappee -- YC 9/4/2014
	  MPI::COMM_WORLD.Recv(&mhtermB, 1, MPI::DOUBLE, procIdForB, 235, status);
	  MPI::COMM_WORLD.Recv(&bbeta, 1, MPI::DOUBLE, procIdForB, 201, status);
	} 
      catch (MPI::Exception e) 
	{
	  std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
	  return 0;
	}
	
      // Compute the rank of temperature of swappee
      rankTemp2 = static_cast<int> (1/bbeta-1)/10;  

      // Compute the MH ratio 
      // double mhterm = exp(mhtermA * abeta - mhtermB * bbeta);
      // Modified by YC - 12/9/2014
      double mhterm = exp((abeta - bbeta) *(mhtermB - mhtermA));
      if (mhterm >= 1.0 || mhterm > runiform()) 
	{
	  ChangeMCMCState(whichElementA, bbeta);
	  chains.at(whichElementA).compute_rankTemp();
	  swapvar = 1;
	} 
      else 
	{
	  swapvar = 0;
	}
      
      // Send the accept/reject result to swappee -- YC 9/4/2014
      try 
	{
	  MPI::COMM_WORLD.Send(&swapvar, 1, MPI::INT, procIdForB, 301);
	  if(swapvar == 1)
	    {
	      // Send the temperature to swappee -- YC 9/4/2014
	      MPI::COMM_WORLD.Send(&abeta, 1, MPI::DOUBLE, procIdForB, 200);
	    }
	} 
      catch (MPI::Exception e) 
	{
	  std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
	  return 0;
	}

      // Compute the local acceptance rate for chain swap
      compute_localAcceptanceRate_chainSwap(rankTemp1,rankTemp2,swapvar);


    } // END of else if(swapAflag == 1)
  else if(swapBflag == 1)
    {
      // On the process with swappee, it send the temperature and mhterm to swapper
      // and then receive the accept/reject result from the swapper and the temperature if accepted.
      MPI::Status status;
      double mhtermB = 0.0;

      // Obtain MHterm and temperature
      mhtermB = GetMHSwapTerm(whichElementB);
      bbeta = chains.at(whichElementB).GetTemperature();
      try 
	{
	  // Send the MHterm and temperature to swapper
	  MPI::COMM_WORLD.Send(&mhtermB, 1, MPI::DOUBLE, procIdForA, 235);
	  MPI::COMM_WORLD.Send(&bbeta, 1, MPI::DOUBLE, procIdForA, 201);
	}
      catch (MPI::Exception e) 
	{
	  std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
	}

      // Receive the accept/reject result from swapper
      try 
	{
	  MPI::COMM_WORLD.Recv(&swapvar, 1, MPI::INT, procIdForA, 301, status);
	} 
      catch (MPI::Exception e) 
	{
	  std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
	  return 0;
	}
      
      if(swapvar == 1)
	{
	  // Receive the temperature of swapper
	  try 
	    {
	      MPI::COMM_WORLD.Recv(&abeta, 1, MPI::DOUBLE, procIdForA, 200, status);
	    } 
	  catch (MPI::Exception e) 
	    {
	      std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
	      return 0;
	    }
	  // Swap the temperature
	  ChangeMCMCState(whichElementB, abeta);
	  chains.at(whichElementB).compute_rankTemp();
	}

    } // END of else if(swapBflag == 1)


  return swapvar;	
}
/* End of function MCMC::SwapChainHeats_BwProcesses() */




/**
 * \function MCMC::FilePrintMCMCState File output function
 * \param fp1 Output stream or pointer to file 1 - all iteration writes
 * \param fp2 Output stream or pointer to file 2 - only cold iteration writes
 * \param swapping Integer flag if swapping is allowed or not
 * \return NULL
 */
void MCMC::FilePrintMCMCState(std::ofstream &fp1, std::ofstream &fp2, int swapping)
{
	double post = 0.0;
	double b = 0.0;
	unsigned int i = 0;
	unsigned int j = 0;
	for ( j = 0; j < chains.size(); j++) {
		try {
			post = chains.at(j).GetPosterior();
		} catch (std::exception &e) {
			std::cout << "Can't access posterior of " << j << " on process " << process_id << "\n";
		}
		try {
			i = chains.at(j).GetCurrentIteration();
		} catch (std::exception &e) {
			std::cout << "Can't access iteration of " << j << " on process " << process_id << "\n";
		}
		try {
			b = chains.at(j).GetTemperature();
		} catch (std::exception &e) {
			std::cout << "Can't access temperature of " << j << " on process " << process_id << "\n";
		}
		if (b == 1.0) {
			fp2 << i << "\t";
			fp2 << post << "\t";
			if (swapping == 1) {
				fp2 << "swapped\n";
			} else {
				fp2 << "unswapped\n";
			}
		}
		fp1 << i << "\t";
		fp1 << post << "\n";
	}
	return;
}
/* End of function MCMC::FilePrintMCMCState() */

