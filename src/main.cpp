/* MIST copyright 2016 by Yujin Chung and Jody Hey */


#include <iostream>
#include <time.h>
#include <vector>
#include <fstream>
#include <valarray>
#include <string>
#include "misc.hpp"
#include "IM.hpp"
#include "coaltree.hpp"
#include "Chain.hpp"
#include "popTree.hpp"
#include "histograms.hpp"
#include <memory>
#include <stdlib.h>
#include <fstream>
#include <unistd.h> // usleep()
#include "optimization.hpp"


using namespace std;

MersenneTwister mt;


/**
 * Main function
 */
int main (int argc, char *argv[])
{
  
  unsigned int numprocesses = 1;
  
  #ifdef MPI_ENABLED
  MPI::Init(argc, argv);
  MPI::COMM_WORLD.Set_errhandler(MPI::ERRORS_THROW_EXCEPTIONS);
  MPI::Status status;
#endif
  int *swapper;
  int *swappee;


  /// YC 5/29/2014
  /// the ID of the current process. Chains on the same process have the same 'currentid'.
  int currentid = 0; 


  // Get the number of processes - YC 9/26/2014
#ifdef MPI_ENABLED
  try 
    {
      numprocesses = MPI::COMM_WORLD.Get_size();
    } catch (MPI::Exception e) {
    std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
    MPI::COMM_WORLD.Abort(-1);
    return 0;
  }

  // Asign process IDs. - YC 9/26/2014
  try 
    {
      currentid = MPI::COMM_WORLD.Get_rank();
    } catch (MPI::Exception e) {
    std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
    MPI::COMM_WORLD.Abort(-1);
    return 0;
  }
  #endif
  
  clock_t clock1, clock2, clock3, clock4;

  // Chains have different random seeds - fixed by YC 9/11/2014
  mt.init_genrand(time(NULL)+time(NULL)*currentid);

  
  /// YC 5/29/2014
  /// Headnode has currentid = 0.
  if(currentid == 0)
    {
      // Welcome message
      cout << "Welcome!\n\n";
    }
  
	    
  IM im;
  // Initialization of 'im' includes reading data set
  // computing the empirical distribution etc.
  
  unsigned int execute=im.initialization(argc, argv, currentid);
  if(execute==1)
    {
      
      
      // popTree* poptree = new popTree;
      // poptree->initialization(im);
      
      unsigned int MLmodes = im.get_MLmodes();	
      
      /*
	if(currentid ==0 && (MLmodes ==0 || MLmodes == 1) && im.get_nLoci() <=100)
	{
	std::cout << "Input data summary:\n";
	im.print_data_pi();
	}
      */
      
      
      // FIXME YC 3/31/2014
      // We may replace 'n_totalgen' by 'n_treeSample' if we
      // want to save the sampled trees after burning and thinning.
      // im.initialize_MCMCsetting(n_chains,n_treeSample, samplingPrior);//,nBurning, thinning);
      unsigned int n_treeSample = im.get_nSampledTrees();
      unsigned int n_chains = im.get_nChains();
      unsigned int numLoci = im.get_nLoci();
      unsigned int samplingPrior = im.get_samplingFromPriorOnly();
      unsigned int lociInParallel = im.get_lociInParallel();
      
      unsigned int totalNumChains=n_chains;
      if(lociInParallel == 0)
	totalNumChains *= numprocesses;
      
      
      /////////////////////
      /// M mode - MCMC ///
      /////////////////////
      
      // Check the starting time
      clock_t clock_start_Mmode = clock();
      
      
      if(MLmodes == 1)
	{
	  
	  MCMC mcmc_run; // An independent run of mcmc chains.
	  // For several run, we can create a vector of 'mcmc'.
	  
	  unsigned int nBurning = im.get_nBurning();
	  unsigned int thinning =im.get_thinning();
	  unsigned int printFreq = im.get_printFreq();
	  
	  mcmc_run.SetUpperKappa(10.0);
	  if(currentid ==0)
	    {
	      std::cout << "#############################\n";
	      std::cout << "#  Markov Chain Monte Carlo #\n";
	      std::cout << "#############################\n\n";
	      std::cout << "Settings:\n\t";
	      if(lociInParallel == 0)
		std::cout << numprocesses << " CPUs and "<< n_chains << " chains per CPU ("<< totalNumChains <<" chains in total)\n\n";	
	      else
		std::cout << numprocesses << " CPUs and " << totalNumChains <<" chains in total\n\n";
	      // std::cout << "Initial trees:\n\n";
	    }

 
    
	  /// added by YC 5/29/2014
	  if(samplingPrior == 0)
	    {
	      mcmc_run.InitializeMCMC(currentid, im, numprocesses);
	    }
	  else
	    {
	      mcmc_run.InitializeMCMC_usePriorOnly(currentid, im, numprocesses);		
	    }
	  if(currentid == 0 && totalNumChains >=2)
	    {
	      std::cout << "The incremental heating values (temperature = 10) for Metropolis-coupling:\n";
	      for(unsigned int i=0; i<totalNumChains; i++)
		{
		  double heat = 1/(1+10.0 * (double) i);
		  std::cout << heat <<", ";
		}
	      std::cout << "\n";
	    }
	  
	  if(currentid ==0)
	    {
	      std::cout << "Starting MCMC simulation...\n\n";
	    }
	  
	  // Added by YC 5/29/2014
#ifdef MPI_ENABLED
	  if (totalNumChains > 1) 
	    {
	      swapper = new int[nBurning];
	      swappee = new int[nBurning];
	      for (int i = 0; i < nBurning; i++)
		{
		  // Modified by YC 8/13/2014
		  // -- YC added the below condition statement "if".
		  //    The swapper and swappee need to be assigned once on the head node, 
		  //    because they are broadcasted.
		  // -- Replace "rand()" by "runiform()" using mt.cc
		  if (currentid == 0) 
		    {
		      swapper[i] = int(runiform() * (totalNumChains));
		      swappee[i] = int(runiform() * (totalNumChains));
		      while (swapper[i] == swappee[i]) 
			{
			  swappee[i] = int(runiform() * (totalNumChains));
			}
		    }
		  
		  MPI::COMM_WORLD.Bcast(&swapper[i], 1, MPI_INT, 0);
		  MPI::COMM_WORLD.Bcast(&swappee[i], 1, MPI_INT, 0);	
		}
	    }	    
#endif
	  
	  //---- Burn-in period ----//
	  if(currentid ==0)
	    {
	      std::cout << "Burn-in: " << nBurning <<" iterations \n";
	    }
	  
	  for(unsigned int i=0; i<nBurning; i++)
	    {
	      // Update each chain
	      if (totalNumChains > 1 ) 
		mcmc_run.UpdateAllChains(0,im.getLoci(),0, swapper[i], swappee[i]);
	      else
		mcmc_run.UpdateAllChains(0,im.getLoci(),0, 0, 0);
	      
	      
	      /*
		
		#ifdef MPI_ENABLED
		MPI::COMM_WORLD.Barrier();
		#endif	   
		
		// print the log-likelihoods of all chains
		if(currentid ==0 && i==0)
		{
		std::cout << "Iter Log-likelihoods of chains\n";
		}
		
		if(i == printFreq* static_cast<int>(i/(printFreq)))
		{
		  if(currentid == 0)
		    {
		      std::cout << i << " ";
		    }
		  mcmc_run.print_allLogLik();
		}
	  
		#ifdef MPI_ENABLED
		MPI::COMM_WORLD.Barrier();
		#endif	  
		
		if(numLoci <=100)
		{
		if(i != 0 && i == 10*printFreq* static_cast<int>(i/(10*printFreq)) )
		{
		// print the cold chain state	  
		if(currentid ==0)
		{
		std::cout << "\niter = " << i <<"\n";
		      std::cout << "temperature log(Priortrees) Locus log-lik Kappa mutScaler Tree\n";
		      }
		      mcmc_run.print_coldChainState();
		      
		      MPI::COMM_WORLD.Barrier();
		      }
		      }
	      */
	    } // END of Burning
	  
	  if(currentid ==0)
	    {
	      std::cout << "Sampling: " << n_treeSample <<" trees \n";
	      std::cout <<"Thinning: "<<thinning <<" iterations.\n";
	    }
	  
	  unsigned int nSample = 0;
	  for(unsigned int i=0, j=0; i<(n_treeSample-1)*thinning+1; i++, j++)
	    {
	      unsigned int swapA = 0;
	      unsigned int swapB = 0;
	      
	      
	      if ( totalNumChains > 1) 
		{
		  if(currentid == 0)
		    {
		      swapA = rand() % (totalNumChains);
		      swapB = rand() % (totalNumChains);
		      while (swapA == swapB) {
			swapB = rand() % (totalNumChains);
		      }
		    }
#ifdef MPI_ENABLED
		  MPI::COMM_WORLD.Barrier();
#endif
		  MPI::COMM_WORLD.Bcast(&swapA, 1, MPI_INT, 0);
		  MPI::COMM_WORLD.Bcast(&swapB, 1, MPI_INT, 0);
#ifdef MPI_ENABLED
		  MPI::COMM_WORLD.Barrier();
#endif
		}
	      
	      
	      if(i==0) // saving
		{
		  if ( totalNumChains > 1) 
		    mcmc_run.UpdateAllChains(0, im.getLoci(), 1, swapA, swapB);
		  else
		    mcmc_run.UpdateAllChains(0, im.getLoci(), 1, 0, 0);
		  
		  mcmc_run.collectUpdates(0);      
		  
		}
	      else if(j==thinning) // saving
		{
		  nSample++;
		  if ( totalNumChains > 1) 
		    mcmc_run.UpdateAllChains(nSample, im.getLoci(),1,swapA, swapB);	    
		  else
		    mcmc_run.UpdateAllChains(nSample, im.getLoci(),1,0, 0);	
		  
		  mcmc_run.collectUpdates(nSample);   
		  
		  j=0;
		  
		}
	      else // thinning
		{
		  if ( totalNumChains > 1) 
		    mcmc_run.UpdateAllChains(nSample, im.getLoci(),0,swapA, swapB);	     
		  else
		    mcmc_run.UpdateAllChains(nSample, im.getLoci(),0,0,0);
		  
		}
	 
	      
	      // print the log-likelihoods of all chains
	      /*
		if(currentid ==0 && i==0)
		{
		std::cout << "Iter Log-likelihoods of chains\n";
		}
		
		if(i == printFreq* static_cast<int>(i/(printFreq)))
		{	      
		if(currentid == 0)
		{
		std::cout << i << " ";
		}	      
		mcmc_run.print_allLogLik();
		}

		#ifdef MPI_ENABLED
		MPI::COMM_WORLD.Barrier();
		#endif	  
		if(numLoci <=100 && ((i != 0 && i == 10*printFreq* static_cast<int>(i/(10*printFreq)) )||i==(n_treeSample-1)*thinning))
		{
		// print the cold chain state	  
		if(currentid ==0)
		{
		  std::cout << "\niter = " << i <<"\n";
		  std::cout << "temperature log(Priortrees) Locus log-lik Kappa mutScaler Tree\n";
		}
		mcmc_run.print_coldChainState();
		
		MPI::COMM_WORLD.Barrier();
		
		// compute and print tree acceptance rate
		MPI::COMM_WORLD.Barrier();
		if(lociInParallel ==0)
		{
		  mcmc_run.computeGlobalTreeAcceptanceRate();
		  MPI::COMM_WORLD.Barrier();
		  mcmc_run.compute_globalAcceptanceRate_chainSwap();
		  }
		  if(lociInParallel ==1)
		  {
		  if(currentid ==0)
		  std::cout << "\nTree acceptance rate (%)\n";	    
		  mcmc_run.print_globalTreeAcceptanceRate();
		  }
		  else
		  {
		  if(currentid ==0)
		  {
		  std::cout << "\nTree acceptance rate (%)\n";	    
		      mcmc_run.print_globalTreeAcceptanceRate();
		      if(totalNumChains > 1 && lociInParallel ==0)
		      {
		      std::cout << "Rate of global chain swap (%) in the upper triangular part\n"
		      << "and the number of swap tries in the lower triangular part\n";
		      mcmc_run.print_AccepRate_swapChains_global();
		      }
		      }
		      std::cout << "\n";		    
		      }
		      MPI::COMM_WORLD.Barrier();
		      }
	      */
	    }// End of sampling

	  // compute and print tree acceptance rate
	  if(lociInParallel ==0)
	    {
	      mcmc_run.computeGlobalTreeAcceptanceRate();
#ifdef MPI_ENABLED
	      MPI::COMM_WORLD.Barrier();
#endif	  
	      mcmc_run.compute_globalAcceptanceRate_chainSwap();
	    }
	  if(lociInParallel ==1)
	    {
	      if(currentid ==0)
		std::cout << "\nTree acceptance rate (%)\n";	    
	      mcmc_run.print_globalTreeAcceptanceRate();
	    }
	  else
	    {
	      if(currentid ==0)
		{
		  std::cout << "\nTree acceptance rate (%)\n";	    
		  mcmc_run.print_globalTreeAcceptanceRate();
		  if(totalNumChains > 1 && lociInParallel ==0)
		    {
		      std::cout << "Rate of global chain swap (%) in the upper triangular part\n"
				<< "and the number of swap tries in the lower triangular part\n";
		      mcmc_run.print_AccepRate_swapChains_global();
		    }
		}
	      std::cout << "\n";		    
	    }
	  MPI::COMM_WORLD.Barrier();
	  if(currentid ==0)
	    std::cout << "\nThe end of MCMC simulation\n\n";
	  
	  
	  // Check the end time of M mode
	  
	  if (currentid == 0) 
	    {
	      mcmc_run.saveCoalTimes2File_chain0();
	      if(im.get_LikelihoodModelType()==2)
		mcmc_run.saveKappa2File_chain0();
	      if(im.get_nLoci() > 1 && im.get_multiLocusSpecific_mutationRate()==2)
		mcmc_run.saveMutationScalers2File_chain0();
	      mcmc_run.saveTipIDs_chain0();
	      mcmc_run.saveTreeIDs_chain0();
	      mcmc_run.saveListTopo2File_chain0();
	      mcmc_run.saveLogEachLikelihoods_chain0();
	      mcmc_run.saveEachLogPrior_chain0();
	      // mcmc_run.saveLogJointLikelihoods_chain0();
	      // mcmc_run.saveLogPriorTrees_chain0();
	      // mcmc_run.saveRNormalDist_chain0();	  
	    }
      
	  clock_t clock_end_Mmode = clock();
	  float runTime1  = (float) (clock_end_Mmode-clock_start_Mmode)/CLOCKS_PER_SEC;
	  if(currentid == 0 )
	    std::cout << "\nThe running time of step 1 (MCMC) is " << runTime1 <<" seconds.\n";
	  MPI::COMM_WORLD.Barrier();
	  
	  /*
	  if(currentid ==0)
	    {
	      coldCh.deleteAll();
	    }
	  */
      
#ifdef MPI_ENABLED
	  if (numprocesses > 1 ) 
	    {
	      delete[] swapper;
	      delete[] swappee;
	    }
	  mcmc_run.deleteTrees();
#endif
	}// END of if(MLmodes ==0 || MLmodes == 1)

  
  
  
      // L mode - Approximation of the posterior density of parameters
    
      if(MLmodes >= 2 ) 
	{	  
	  if(currentid ==0)
	    {
	      std::cout <<"#############################\n";
	      std::cout << "  Step 2: model inference.\n";
	      std::cout <<"#############################\n";
	      std::cout << numprocesses << " CPUs\n\n";	
	    }
	  popTree* poptree = new popTree;
	  Chain coldCh;
	  coldCh.initializeLmode(im, currentid, numprocesses);

	  coldCh.read_LmodeInputFile(im);
	  
	  // std::cout <<"Done with reading input files.\n\n";
	  
	  // Check the start time of L mode
	  clock_t clock_start_Lmode = clock();
	  
	  	  
	  poptree->initialize_popTree(im, currentid); 
	  // coldCh.prepare_Lmode(poptree);
	  clock_t clock_end = clock();
	  
	  /*
	  if(currentid ==0)
	    {
	      double runTime_part = (double) (clock_end-clock_start_Lmode)/CLOCKS_PER_SEC;
	      std::cout <<"Done in preparing L mode. (running time = "<< runTime_part <<" sec)\n";
	    }
	  */
    
	  /// Marginal histograms of demographic parameters
	  clock_t clock_start = clock();
	  
	  // Finding the maximum a poterior probability estimate
	  if( MLmodes >= 2)
	    {
	      MaxPosterior MAPestimate;
	      
	      // if(im.get_migRateMax()!=0)
		coldCh.prepare_Lmode(poptree, im);
	      
	      // std::cout <<"Done with reading input files.\n";

	      if(MLmodes == 3)
		{
		  Eigen::MatrixXd paraVector(1,6);
		  unsigned int nPara_popSizes = 1;
		  // sampling populations
		  paraVector(0,0) = im.get_truePara()(0,0);
		  if(im.get_samePopulationSizes() == 1)
		    paraVector(0,1) = im.get_truePara()(0,0);
		  else
		    {
		      paraVector(0,1) = im.get_truePara()(0,1);
		      nPara_popSizes++;
		    }
		  // ancestral population
		  if(im.get_ancPop() == 1)
		    {
		      if(im.get_samePopulationSizes() == 1)
			paraVector(0,2) = im.get_truePara()(0,0);
		      else
			{
			  paraVector(0,2) = im.get_truePara()(0,2);
			  nPara_popSizes++;
			}
		    }
		  // migration rates
		  paraVector(0,3) = im.get_truePara()(0,nPara_popSizes);
		  if(im.get_sameMigrationRates() == 1)
		    paraVector(0,4) =im.get_truePara()(0,nPara_popSizes);
		  else
		    paraVector(0,4) =im.get_truePara()(0,nPara_popSizes+1);
		  //-- splitting time --//
		  if(im.get_ancPop() == 1)	
		    paraVector(0,5) = im.get_truePara()(0,(im.get_truePara()).cols()-1);
		  else // no ancestral population
		    paraVector(0,5) =im.get_splittingTimeMax();
		  
		  // double truePosterior = MAPestimate.computeJointDensity_MPI_overSubSample(paraVector,im, poptree, coldCh, numprocesses,currentid);
		  double truePosterior = MAPestimate.computeLogJointDensity_MPI_overSubLoci_ESS(paraVector, im, poptree, coldCh, numprocesses, currentid);  
		  if(currentid == 0)
		    {
		      std::cout << "True parameters: "<< im.get_truePara() <<"\n";
		      std::cout << "The posterior density (log) of the true parameters is " << truePosterior <<"\n\n";
		    }
		}// END of if(MLmodes == 3)
	      if(MLmodes == 2) // || MLmodes==4)
		{
		  if(currentid == 0)
		    {
		      std::cout <<"\nPreparing the optimization (differential evolution)\n";
		    }	  
		  // coldCh.prepare_Lmode(poptree);
		  MAPestimate.initiate(im, poptree, coldCh, numprocesses, currentid);	  
		  if(currentid == 0)
		    {
		      std::cout <<"\nStarting the optimization (differential evolution)\n";
		    }
		  MAPestimate.DE(im, poptree, coldCh, numprocesses, currentid);	  
		  if(currentid == 0)
		    {
		      std::cout <<"\n\nDone with the optimization (differential evolution)\n\n"; 
		      
		      // The total computing time for getting eigen values and eigen vectors
		      /*
		      std::cout << "The total computing time for getting eigen values and eigen vectors: "
				<<  MAPestimate.get_totalComputingTime_eigen().count()/1000000 <<"(sec)\n";
		      std::cout << "The function was called " << MAPestimate.get_totalNum_eigenFunctionCalls() << " times\n";
		      std::cout << "The total computing time for computing the coalescent conditional probabilities: "
				<<  MAPestimate.get_totalComputingTime_condiProb().count()/1000000  <<"(sec)\n";
		      std::cout << "The function was called " << MAPestimate.get_totalNum_condiProbFunctionCalls() << " times\n";
		      */
		    }
		}// END of if(MLmodes == 3)
	  
	      if(MLmodes==5)
		{
		  Eigen::MatrixXd paraVector(1,6);
		  unsigned int nPara_popSizes = 1;
		  // sampling populations
		  paraVector(0,0) = im.get_truePara()(0,0);
		  if(im.get_samePopulationSizes() == 1)
		    paraVector(0,1) = im.get_truePara()(0,0);
		  else
		    {
		      paraVector(0,1) = im.get_truePara()(0,1);
		      nPara_popSizes++;
		    }
		  // ancestral population
		  if(im.get_ancPop() == 1)
		    {
		      if(im.get_samePopulationSizes() == 1)
			paraVector(0,2) = im.get_truePara()(0,0);
		      else
			{
			  paraVector(0,2) = im.get_truePara()(0,2);
			  nPara_popSizes++;
			}
		    }
		  // migration rates
		  paraVector(0,3) = im.get_truePara()(0,nPara_popSizes);
		  if(im.get_sameMigrationRates() == 1)
		    paraVector(0,4) =im.get_truePara()(0,nPara_popSizes);
		  else
		    paraVector(0,4) =im.get_truePara()(0,nPara_popSizes+1);
		  //-- splitting time --//
		  if(im.get_ancPop() == 1)	
		    paraVector(0,5) = im.get_truePara()(0,(im.get_truePara()).cols()-1);
		  else // no ancestral population
		    paraVector(0,5) =im.get_splittingTimeMax();
		  double posterior = MAPestimate.computeJointDensity_MPI_overSubSample(paraVector,im, poptree, coldCh, numprocesses,currentid);
		  if(currentid == 0)
		    {
		      std::cout << "Parameters: "<< im.get_truePara() <<"\n";
		      std::cout << "The posterior density of the parameters is " << posterior <<"\n\n";
		    }
		}// END of if(MLmodes==5)
	    }// END of if( MLmodes >= 2)

      
	  // Computing the marginal distribution
	  if( MLmodes == 4) 
	    {
	      
	      marginal marDensity;
	      marDensity.set_gridsize_nJointPoints(im.get_gridsize(),im.get_nJointPoints());
	      if(currentid ==0){
		std::cout << "Preparing computing histograms\n";
	      }
	      marDensity.prepareComputingHistrograms_allProc(im);
	      
	      // new version: parallized over MCMC sample -- by YC 11/18/2014
	      marDensity.prepareComputingHistrograms(im);
	      clock_end = clock();
	      if(currentid ==0)
		{
		  double runTime_part = (double) (clock_end-clock_start)/CLOCKS_PER_SEC;
		  std::cout <<"Done. (running time = "<< runTime_part <<" sec)\n";
		}
	      
	      
	      
	      clock_start = clock();
	      std::cout << "Computing the joint density of demographic parameters on process "<<currentid <<" ...\n";      
	      // new version: parallized over MCMC sample -- by YC 11/18/2014
	      marDensity.computeJointDensity_MPI_overSample(im,poptree,coldCh,numprocesses, currentid);
	      // old version: parallelized over demographic parameter values
	      //marDensity.computeJointDensity_MPI(im,poptree,coldCh,numprocesses, currentid);
	      clock_end = clock();
	      double runTime_part = (double) (clock_end-clock_start)/CLOCKS_PER_SEC;
	      std::cout <<"Done with computing the joint density of demographic parameters on process "<< currentid <<"(running time = "<< runTime_part <<" sec)\n";
	      
#ifdef MPI_ENABLED
	      MPI::COMM_WORLD.Barrier();
#endif
	      
	      
	      
	      clock_start = clock();
	      if(currentid ==0)
		{
		  //marDensity.computeJointDensity(im,poptree,coldCh);
		  std::cout << "Computing the marginal distribution of each demographic parameter...\n";
		}
	      // new version: parallized over MCMC sample -- by YC 11/18/2014
	      marDensity.computeMarginalDensity_MPI_overSample(im,poptree,coldCh,numprocesses, currentid);
	      // old version: parallelized over demographic parameter values
	      // marDensity.computeMarginalDensity_MPI(im,poptree,coldCh,numprocesses, currentid);
	      clock_end = clock();
	      if(currentid ==0)
		{
		  double runTime_part = (double) (clock_end-clock_start)/CLOCKS_PER_SEC;
		  std::cout <<"Done with computing the marginal distribution of each demographic parameter.(running time = "<< runTime_part <<" sec)\n";

		  marDensity.saveHistogram();
		  marDensity.printHistogram();

		  std::cout << "Computing the posterior mean and 95% interval of each demographic parameter...";
		  marDensity.computeMeanIntervals();
		  std::cout <<"Done.\n";
		  marDensity.saveMeanIntervals();
		  //marDensity.printMeanIntervals();
		  
		  
		  //coldCh.deleteAll();
		  
		} // End of if(currentid ==0)
	      
	    } // End of if(MLmodes ==1|| MLmodes ==2)
	
	  // Check the end time of L mode.
	  clock_t clock_end_Lmode = clock();
	  
	  if(currentid ==0)
	    {
	      // std::cout << "\nEnd of L mode.\n\n";
	      float runTime2  = (float) (clock_end_Lmode- clock_start_Lmode)/CLOCKS_PER_SEC;
	      std::cout <<  "The running time of step 2 is " << runTime2 <<" seconds.\n";
	    }
	  
	  coldCh.deleteAll_Lmode();

	  // Delete heap objects
	  poptree->deletePopTree();
	  delete poptree;
	  
	}// END of L mode
    


    }// END of if(execute==1)

  
	
#ifdef MPI_ENABLED
  MPI::Finalize();
  //AS: debug only
  //f1.close();
#endif
  
  
  return 0;
  
}  /// End of main
