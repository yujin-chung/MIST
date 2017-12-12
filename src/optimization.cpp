/* MIST copyright 2016 by Yujin Chung and Jody Hey */


#include <iostream>
#include <math.h>
#include "optimization.hpp"
#include "cppoptlib/solver/lbfgsbsolver.h"



void MaxPosterior::initiate(IM im, popTree* poptree, Chain coldCh, unsigned int nProcs, unsigned int crr_procID)
{
  optimizer = im.get_optimizer();

  //-- Settings for differential evoltuion --//
  nIters = 1000;
  if(optimizer==0) // differential evolution
    {
      nParaVectors = 30;
      F = 0.8;
      CR = 0.9;   
    }
  else if(optimizer==1) // L-BFGS-B
    {
      nParaVectors = 1;
    }

  if(optimizer==0) // DE
    nParaVectors = im.get_nParaVectors();

  //  std::chrono::microseconds initial(0);
  // totalComputingTime_eigen = initial;
  // totalComputingTime_eigen_subMatOfEigenVectors= initial;
  // totalComputingTime_condiProb = initial;

  checkpoint = im.get_checkpoint();
  if(checkpoint ==1|| checkpoint==3)
    {
      howOften_checkpoint = im.get_howOften_checkpoint();
    }

  unsigned int ancPop = im.get_ancPop();
  unsigned int samePopulationSizes = im.get_samePopulationSizes();
  unsigned int sameMigrationRates = im.get_sameMigrationRates();
  double migRateMax = im.get_migRateMax();
  unsigned int nPara_popSizes;
  

  if(poptree->get_age()==0) // single population
    {
      nPara =1;
      nPara_popSizes =1;
    }
  else // more than or equal to 2 populations
    {
      if(ancPop == 0) // island model (no ancestral populations)
	{
	  if(samePopulationSizes == 1)
	    {
	      nPara = 1;
	      nPara_popSizes = 1;
	    }
	  else
	    {
	      nPara = 2;
	      nPara_popSizes = 2;
	    }
	}
      else if(ancPop == 1)
	{
	  if(samePopulationSizes == 1)
	    {
	      nPara = 2; // population size and splitting time
	      nPara_popSizes = 1;
	    }
	  else
	    {
	      nPara = 4; // three population sizes and a splitting time
	      nPara_popSizes = 3;
	    } 
	}
      if(migRateMax !=0)
	{
	  if(sameMigrationRates ==1)
	    nPara += 1;
	  else
	    nPara += 2;     
	}
    }

    
  para_atCrr.resize(nParaVectors,nPara);
  para_atPrev.resize(nParaVectors,nPara);
  posterior_atCrr.resize(nParaVectors);
  posterior_atPrev.resize(nParaVectors);
  
  for(unsigned int i=0; i<nPara; i++)
    {
      if(poptree->get_age() ==0)
	{	  
	  priorsMax.push_back(im.get_popSizeMax());
	}
      else{
	if(i < nPara_popSizes)
	  priorsMax.push_back(im.get_popSizeMax());
	else if(i==nPara-1 && ancPop ==1)
	  priorsMax.push_back(im.get_splittingTimeMax());
	else 
	  priorsMax.push_back(im.get_migRateMax());
      }
    }

  /*
  if(crr_procID == 0)
    marginals.initiate_forDE(im, nPara, priorsMax);
  */
  // coldCh.prepare_Lmode(poptree);


  // Initialize "para_atPrev"
  for(unsigned int i=0; i<nParaVectors; i++)
    {
      if(checkpoint ==2|| checkpoint ==3)
	{
	  // read the checkpoint and save them to "para_atPrev"
	  read_checkpoint();
	}
      else if(checkpoint <= 1)
	{
	  // initialize "para_atPrev" (individuals in DE)
	  for(unsigned int j=0; j<nPara; j++)
	    {
	      double para_temp = 0.0;
	      if(crr_procID == 0)
		para_temp = priorsMax.at(j)*runiform();   
	      
	      MPI::COMM_WORLD.Barrier();
	      MPI::COMM_WORLD.Bcast(&para_temp, 1, MPI_DOUBLE, 0);
	      MPI::COMM_WORLD.Barrier();
	      para_atPrev(i,j)=para_temp;
	    }
	}
      else
	{
	  std::cout << "\n*** Error in void MaxPosterior::initiate() ***\n"
		    << "checkpoint should take a value of 0, 1, 2 or 3, but checkpoint = " 
		    << checkpoint <<".\n\n";
	}
    }

  // Computing the posterior of para_atPrev
  for(unsigned int i=0; i<nParaVectors; i++)
    {      
      Eigen::MatrixXd paraVector(1,6);
      if(poptree->get_age()==0)
	{
	  for(unsigned int j=0; j<6; j++)
	    {
	      if(j==2) // ancestral population
		paraVector(0,j) = para_atPrev(i,0);
	      else
		paraVector(0,j) = 0;
	    }
	}
      else
	{
	  // sampling populations
	  paraVector(0,0) = para_atPrev(i,0);
	  if(samePopulationSizes == 1)
	    paraVector(0,1) = para_atPrev(i,0);
	  else
	    paraVector(0,1) = para_atPrev(i,1);
	  // ancestral population
	  if(ancPop == 1)
	    {
	      if(samePopulationSizes == 1)
		paraVector(0,2) = para_atPrev(i,0);
	      else
		paraVector(0,2) = para_atPrev(i,2);
	    }
	  // migration rates
	  if(migRateMax !=0)
	    {
	      paraVector(0,3) = para_atPrev(i,nPara_popSizes);
	      if(sameMigrationRates == 1)
		paraVector(0,4) = para_atPrev(i,nPara_popSizes);
	      else
		paraVector(0,4) = para_atPrev(i,nPara_popSizes+1);
	    }
	  else
	    {
	      paraVector(0,4) = 0; paraVector(0,5) =0;
	    }
	  //-- splitting time --//
	  if(ancPop == 1)	
	    paraVector(0,5) = para_atPrev(i,nPara-1);	    
	  else // no ancestral population
	    paraVector(0,5) =im.get_splittingTimeMax();
	}

     
      unsigned int lociInParallel = im.get_lociInParallel();

      if(lociInParallel ==1)
	{
	  if(coldCh.get_multiLocusSpecific_mutationRate() ==1) // variable mutation rate scalars
	    {
	      posterior_atPrev.at(i) = computeLogJointDensity_mutationScalars_MPI_overSubLoci(paraVector, im, poptree, coldCh, nProcs, crr_procID);  	      
	    }
	  else // constant mutation rate scalars
	    {
	      posterior_atPrev.at(i) = computeLogJointDensity_MPI_overSubLoci(paraVector, im, poptree, coldCh, nProcs, crr_procID); 
	    }
	}
      else if(lociInParallel ==0)
	{
	  posterior_atPrev.at(i) =  log(computeJointDensity_MPI_overSubSample(paraVector,im, poptree,coldCh, nProcs, crr_procID));
	}
      else
	{
	  std::cout << "lociInParallel should be 1 or 0, but lociInParallel = " << lociInParallel <<"\n";
	}
    }
  
 
  return;
}



/*
template <> double computeLogJointDensity<double, 6>::computeLogJointDensity_mutationScalars_MPI_overSubLoci(TVector paraVec, MaxPosterior MAPestimate, IM im, popTree* poptree, Chain coldCh, unsigned int nProcs, unsigned int crr_procID)
{
  double logJointDensity = 0;

  Eigen::MatrixXd para;
  para.resize(1, MAPestimate.get_nPara());
  for(unsigned int i=0; i< paraVec.size(); i++)
    para(0,i) = paraVec(i);

  logJointDensity =  -1 * (double) MAPestimate.computeLogJointDensity_mutationScalars_MPI_overSubLoci(para, im, poptree, coldCh, nProcs, crr_procID);  	

  // std::cout <<" paraVect = " << paraVec <<"\n";
  std::cout <<"para = "<< para <<"\n";
  std::cout <<" logJointDensity = "<< logJointDensity <<"\n\n";
  

  return logJointDensity;
}
*/

void MaxPosterior::read_checkpoint()
{
  ifstream inFile;
  inFile.open("checkpoint_DE.txt");
  double para_temp;
  for(unsigned int i=0; i<nParaVectors; i++)
    {
      for(unsigned int j=0; j<nPara; j++)
	{
	  inFile >> para_temp;	  
	  para_atPrev(i,j)=para_temp;
	}
    }
  
  inFile.close();

  return;
}

void MaxPosterior::write_checkpoint()
{
  ofstream outFile;
  outFile.open("checkpoint_DE.txt");
  
  double para_temp;
  for(unsigned int i=0; i<nParaVectors; i++)
    {
      for(unsigned int j=0; j<nPara; j++)
	{
	  para_temp = para_atPrev(i,j);	  
	  if(j==nPara-1)
	    outFile << para_temp <<"\n";
	  else
	    outFile << para_temp <<"\t";
	}
    }

  outFile.close();

  return;
}





void MaxPosterior::DE_eachIter(IM im, popTree* poptree, Chain coldCh, unsigned int nProcs, unsigned int crr_procID)
{

  Eigen::MatrixXd newPara;
  newPara.resize(1,nPara);
  for(unsigned int i=0; i<nParaVectors; i++)
    {
      unsigned int recombinationID = runiform_discrete(nPara);
      newPara.setZero();
      double newPara_each = 0.0;
      // select base vectors
      unsigned int b1=i;
      unsigned int b2=i;
      unsigned int b3=i;
      if(crr_procID == 0)
	{
	  while(b1==i)
	    b1 = runiform_discrete(nParaVectors);
	  while(b2 ==i || b2 == b1)
	    b2 = runiform_discrete(nParaVectors);
	  while(b3==i || b3== b2|| b3== b1)
	    b3 = runiform_discrete(nParaVectors);
	}
     
      for(unsigned int j=0; j<nPara; j++)
	{
	  if(crr_procID == 0)
	    {
	      newPara_each = para_atPrev(i,j);
	      if(runiform() < CR ||  j== recombinationID)
		{
		  newPara_each = para_atPrev(b1,j)+F*(para_atPrev(b2,j)-para_atPrev(b3,j)); // differential mutation
		  while(newPara_each <0 || newPara_each >priorsMax.at(j))
		    {
		      if(newPara_each <0)
			newPara_each = para_atPrev(b1,j)-runiform()*(para_atPrev(b1,j));
		      if(newPara_each > priorsMax.at(j))
			newPara_each = para_atPrev(b1,j)+runiform()*(priorsMax.at(j)-para_atPrev(b1,j));
		    }
		  
		  if(newPara_each <0 || newPara_each > priorsMax.at(j))
		    {
		      std::cout << "\n*** Error *** in map::DE_eachIter()\n";
		    }
		}
	    }
	  MPI::COMM_WORLD.Barrier();
	  MPI::COMM_WORLD.Bcast(&newPara_each, 1, MPI_DOUBLE, 0);
	  MPI::COMM_WORLD.Barrier();
	  newPara(0,j) = newPara_each;
	} //END of for(unsigned int j=0; j<nPara; j++)

      Eigen::MatrixXd paraVector(1,6);
      unsigned int nPara_popSizes =1;
      if(poptree->get_age() ==0)// single population
	{
	  for(unsigned int ii=0; ii<6; ii++)
	    {
	      if(ii==2)
		paraVector(0,ii) = newPara(0,0);
	      else
		paraVector(0,ii) = 0;
	    }
	}
      else
	{
	  //-- sampling populations --//
	  paraVector(0,0) = newPara(0,0);
	  if(im.get_samePopulationSizes() == 1)
	    paraVector(0,1) = newPara(0,0);
	  else
	    {
	      paraVector(0,1) = newPara(0,1);
	      nPara_popSizes++;
	    }
	  //-- ancestral population --//
	  if(im.get_ancPop() == 1)
	    {
	      if(im.get_samePopulationSizes() == 1)
		paraVector(0,2) = newPara(0,0);
	      else
		{
		  paraVector(0,2) = newPara(0,2);
		  nPara_popSizes++;
		}
	    }
	  else // no ancestral population
	    paraVector(0,2) = 0;
	  //-- migration rates --//
	  paraVector(0,3) = newPara(0,nPara_popSizes);
	  if(im.get_sameMigrationRates() == 1)
	    paraVector(0,4) = newPara(0,nPara_popSizes);
	  else
	    paraVector(0,4) = newPara(0,nPara_popSizes+1);
	  //-- splitting time --//
	  if(im.get_ancPop() == 1)	
	    paraVector(0,5) = newPara(0,nPara-1);
	  else // no ancestral population
	    paraVector(0,5) =im.get_splittingTimeMax();
	}

      long double logPosterior_newPara = 0;

      unsigned int lociInParallel = im.get_lociInParallel();
      if(lociInParallel ==1)
	{
	  if(coldCh.get_multiLocusSpecific_mutationRate() ==1) // variable mutation rate scalars
	    {
	      logPosterior_newPara = computeLogJointDensity_mutationScalars_MPI_overSubLoci(paraVector, im, poptree, coldCh, nProcs, crr_procID);  
	      
	    }
	  else // constant mutation rate scalars
	    {
	      logPosterior_newPara = computeLogJointDensity_MPI_overSubLoci(paraVector, im, poptree, coldCh, nProcs, crr_procID);  
	    }
	}
      else if(lociInParallel ==0)
	{
	  logPosterior_newPara = (long double) log(computeJointDensity_MPI_overSubSample(paraVector, im, poptree, coldCh, nProcs, crr_procID));
	}
      else
	{
	  std::cout << "lociInParallel should be 1 or 0, but lociInParallel = " << lociInParallel <<"\n";
	}
      
      MPI::COMM_WORLD.Barrier();

      // YC 2/27/2014
      // All the processes share the same posterior densities and parameter values.
      if(logPosterior_newPara > posterior_atPrev.at(i))
	{
	  para_atCrr.row(i) = newPara;
	  posterior_atCrr.at(i) = logPosterior_newPara;
	  
	}
      else
	{
	  para_atCrr.row(i) = para_atPrev.row(i);
	  posterior_atCrr.at(i) = posterior_atPrev.at(i);
	} 
      /*
      if(crr_procID == 0)
	{	  
	  marginals.add(newPara, posterior_newPara);
	}
      */
    }

  MPI::COMM_WORLD.Barrier();
  para_atPrev = para_atCrr;
  posterior_atPrev = posterior_atCrr;
	
  if(crr_procID == 0)
    {	  
      std::vector<long double>::const_iterator iter_max, iter_min;
      iter_max = max_element(posterior_atCrr.begin(), posterior_atCrr.end());
      iter_min = min_element(posterior_atCrr.begin(), posterior_atCrr.end());
      logPosteriorMax = *iter_max;
      logPosteriorMin = *iter_min;
      
      /*
      unsigned int found_max = 0;
      ID_max = 0;
      while(ID_max < nParaVectors && found_max==0 )
	{
	  if(*iter_max == posterior_atCrr.at(ID_max))
	    found_max =1;
	  else
	    ID_max++;
	}    
      */
      // std::cout <<"logPosteriorMax = " << logPosteriorMax <<" logPosteriorMin = " << logPosteriorMin <<"\n";
    }
  MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Bcast(&logPosteriorMax, 1, MPI_LONG_DOUBLE, 0);
  MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Bcast(&logPosteriorMin, 1, MPI_LONG_DOUBLE, 0);
  MPI::COMM_WORLD.Barrier();

  
  return;
}


void MaxPosterior::DE(IM im, popTree* poptree, Chain coldCh, unsigned int nProcs, unsigned int crr_procID)
{
  /*
  if(crr_procID == 0)
    std::cout << "In MaxPosterior::DE()\n";
  */

  unsigned int conv = 0;
  unsigned int iter=0;
  while(conv == 0 && iter < nIters)
    {
      //if(crr_procID==0)
	//std::cout << "Starting DE_eachIter\n";
      DE_eachIter(im, poptree, coldCh, nProcs, crr_procID);
      //if(crr_procID==0)
      //std::cout << "Ending DE_eachIter\n";

      // YC 1/12/2015
      // There are two criteria used to determine the convergence
      // the first criterion is same as that in IMa2
      // and the 2nd criterion is added. It is very difficult
      // to satisfy the 1st criterion if the maximum posterior value is very large, 
      // but the 2nd criterion would determine the convergence earlier.
      maxDist = 0;
      if(logPosteriorMax-logPosteriorMin < pow(10,-4) && abs((logPosteriorMax-logPosteriorMin)/logPosteriorMax) < pow(10,-4)) 
	{
	  /*
	  if(crr_procID ==0)
	    {
	      // std::vector<long double>::const_iterator iter_max;
	      // iter_max = max_element(posterior_atCrr.begin(), posterior_atCrr.end());
	      Eigen::MatrixXd maxPara = para_atCrr.row(ID_max);
	      for(unsigned int ii=0; ii<nParaVectors; ii++)
		{
		  if(ii!=ID_max)
		    {
		      double dd = max((para_atCrr.row(ii)-maxPara).maxCoeff(),(maxPara-para_atCrr.row(ii)).maxCoeff()  );
		      if(maxDist <dd)
			maxDist = dd;
		    }
		}
	    }
	  MPI::COMM_WORLD.Barrier();
	  MPI::COMM_WORLD.Bcast(&maxDist, 1, MPI_DOUBLE, 0);
	  MPI::COMM_WORLD.Barrier();
	  
	  if(maxDist < pow(10,-4))
	    {
	      conv = 1;
	    }
	  */
	  conv =1;
	  // std::cout <<"iter= "<<iter <<" logPosteriorMax-logPosteriorMin="<<logPosteriorMax-logPosteriorMin <<" (logPosteriorMax-logPosteriorMin)/logPosteriorMax="<<(logPosteriorMax-logPosteriorMin)/logPosteriorMax <<"\n";
	}

      if(crr_procID == 0)
	{
	  if(iter - 100* static_cast<unsigned int>(iter/100) == 0)
	    {
	      /*
	      unsigned int found_min = 0;	      
	      unsigned int count_min = 0;
	      while( count_min < nParaVectors &&  found_min==0)
		{
		  if(logPosteriorMin==posterior_atCrr.at(count_min))
		    found_min = 1;
		  else
		    count_min++;
		}
	      */
	      std::cout << "\n\niter = " << iter 
			<<": the largest log(posterior) = "<< logPosteriorMax << ", the smallest log(posterior) = " << logPosteriorMin <<"\n";
	      /*
	      if(maxDist !=0)
		std::cout << "maxDist = " << maxDist <<"\n";
	      */
	      std::cout << "Estimates with the largest posterior: " << para_atCrr.row(ID_max) <<"\n";
	      // std::cout << "Miminum a posterior estimates: " << para_atCrr.row(count_min) <<"\n";   
	      /*
	      std::cout << "The total computing time for getting eigen values and eigen vectors: "
			<<  totalComputingTime_eigen.count()/1000000 <<"(sec)\n";
	      std::cout << "The function was called " << totalNum_eigenFunctionCalls << " times\n";
	      std::cout << "The total computing time for computing the coalescent conditional probabilities: "
			<<  totalComputingTime_condiProb.count()/1000000  <<"(sec)\n";
	      std::cout << "The function was called " << totalNum_condiProbFunctionCalls << " times\n";
	      */
	    }
	  else if(iter - 10* static_cast<unsigned int>(iter/10) == 0)
	    {
	      std::cout <<".";
	    }
	}

      // write a checkpoint
      if(crr_procID ==0)
	{
	  if(checkpoint == 1|| checkpoint==3)
	    {
	      if(iter - howOften_checkpoint * static_cast<unsigned int>(iter/howOften_checkpoint) ==0)
		{
		  write_checkpoint();
		}
	    }
	}


      iter++;

    }
  if(conv==0 && crr_procID==0)
    {
      std::cout << "\nWarning: the optimization did not converge.\n";
    }

 
  if(crr_procID == 0)
    {
      std::cout << "\nMaximum a posterior estimates: " << para_atCrr.row(ID_max)
		<< "\nlog(posterior density) = " << logPosteriorMax
		<<" \n at iteration " << iter
		<<"\n";
      /*
      std::cout << "Computing the posterior means...\n";
      marginals.computeMeanIntervals_forDE();
      std::cout << "Done.\n";
      std::cout << "Saving the posterior means..\n";
      marginals.saveMeanIntervals();
       std::cout << "Done.\n";
      //std::cout << "Saving the histograms..\n";
      marginals.saveHistogram_forDE();
      //std::cout << "Done\n";
      */
     
    }
 
  
  return;
}


// 'demographicPara' is a 1x6 matrix.
double MaxPosterior::computeJointDensity_MPI_overSubSample(Eigen::MatrixXd demographicPara, IM im, popTree* poptree, Chain coldCh, unsigned int nProcs, unsigned int crr_procID)
{
  // REMOVE
  //if(crr_procID == 0)
  //  std::cout << "\n In MaxPosterior::computeJointDensity_MPI_overSubSample\n";

  poptree->replacePara(demographicPara);

 
  coldCh.compute_partialJointPosteriorDensity_overSubSample(poptree, im, crr_procID, nProcs);
  std::vector<long double> expectationOfEachCoalProb = coldCh.get_expectationOfCoalProb();
  unsigned int nloci = coldCh.GetNumLoci();
  double posterior = 1;
  if(nloci != expectationOfEachCoalProb.size())
    {
      std::cout << "** Error In MaxPosterior::computeJointDensity_MPI_overSubSample() **\n";
      std::cout << "nloci = " << nloci 
		<< " expectationOfEachCoalProb.size() = " << expectationOfEachCoalProb.size() <<"\n\n";
    }
  for(unsigned int lc=0; lc< nloci; lc++)
    {
      double each_local = expectationOfEachCoalProb.at(lc);
      double each_global =0;
      MPI::COMM_WORLD.Barrier();
      MPI::COMM_WORLD.Allreduce(&each_local, &each_global, 1, MPI_DOUBLE, MPI_SUM);
      MPI::COMM_WORLD.Barrier();
      if(crr_procID ==0)
	{
	  posterior *=each_global;
	}
      MPI::COMM_WORLD.Barrier();
      
    }

  if(crr_procID ==0)
    {	
      Eigen::Vector3d paraMax = im.get_paraMax();
      double priorPopTree = poptree->computeJointPrior(paraMax);
      posterior *=priorPopTree;
    }

  if(crr_procID ==0)
    {
      totalComputingTime_eigen += coldCh.get_eachComputingTime_eigen();
      totalComputingTime_condiProb += coldCh.get_eachComputingTime_condiProb();
      totalNum_eigenFunctionCalls++;
      totalNum_condiProbFunctionCalls += coldCh.get_countCondiProbFunctionCalls();
      // std::cout << "totalComputingTime_eigen = " << totalComputingTime_eigen.count()/1000 << "milliseconds\n";
    }
 

  return posterior;
}

// 3/16/2016
// 'demographicPara' is a 1x6 matrix.
long double MaxPosterior::computeLogJointDensity_MPI_overSubLoci_ESS(Eigen::MatrixXd demographicPara, IM im, popTree* poptree, Chain coldCh, unsigned int nProcs, unsigned int crr_procID)
{
  // REMOVE
  // if(crr_procID == 0)
  //  std::cout << "\n In MaxPosterior::computeJointDensity_MPI_overSubSample\n";

  poptree->replacePara(demographicPara);

  coldCh.compute_partialJointPosteriorDensity_overSubLoci_ESS(poptree, im, crr_procID, nProcs);
  std::vector<long double> logExpectationOfEachCoalProb = coldCh.get_logExpectationOfCoalProb();
  std::vector<long double> logExpectationOfEachCoalProbSquared = coldCh.get_logExpectationOfCoalProbSquared();
  unsigned int numSubLoci = coldCh.getNumSubLoci();

  long double posterior = 1;
  long double posteriorSquared = 1;
  if(numSubLoci != logExpectationOfEachCoalProb.size())
    {
      std::cout << "** Error In MaxPosterior::computeJointDensity_MPI_overSubSample() **\n";
      std::cout << "numSubLoci = " << numSubLoci 
		<< " expectationOfEachCoalProb.size() = " << logExpectationOfEachCoalProb.size() <<"\n\n";
    }
  long double local_logPosterior_subLoci =0;
  // long double local_logPosteriorSquared_subLoci =0;
  long double local_sumESS =0;
  long double global_logPosterior = 0;
  // long double global_logPosteriorSquared = 0;
  long double global_sumESS =0;
  for(unsigned int lc=0; lc< numSubLoci; lc++)
    {
      local_logPosterior_subLoci += logExpectationOfEachCoalProb.at(lc);
      local_sumESS += coldCh.GetNumIterations()*exp(2*logExpectationOfEachCoalProb.at(lc)-logExpectationOfEachCoalProbSquared.at(lc)); 
      // local_logPosteriorSquared_subLoci += logExpectationOfEachCoalProbSquared.at(lc);
      // std::cout<<"crrProcID="<<crr_procID<< "lc = " << lc <<" local_logPosterior_subLoci = "<< local_logPosterior_subLoci 
      //	       <<" local_logPosteriorSquared_subLoci = "<< local_logPosteriorSquared_subLoci <<"\n";
      //std::cout <<"crrProcID="<<crr_procID<< " lc = " << lc <<" nMCMC="<<coldCh.GetNumIterations()
      //		<< " each ESS= "<< coldCh.GetNumIterations()*exp(2*logExpectationOfEachCoalProb.at(lc)-logExpectationOfEachCoalProbSquared.at(lc)) <<"\n";
    }
  
  MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Allreduce(&local_logPosterior_subLoci, &global_logPosterior, 1, MPI_LONG_DOUBLE, MPI_SUM);
  MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Allreduce(&local_sumESS, &global_sumESS, 1, MPI_LONG_DOUBLE, MPI_SUM);
  // MPI::COMM_WORLD.Allreduce(&local_logPosteriorSquared_subLoci, &global_logPosteriorSquared, 1, MPI_LONG_DOUBLE, MPI_SUM);
  MPI::COMM_WORLD.Barrier();
  
  long double logPosterior = 0;
  long double ESS = 0;
  if(crr_procID ==0)
    {	
      Eigen::Vector3d paraMax = im.get_paraMax();
      double priorPopTree = poptree->computeJointPrior(paraMax);
      //posterior = exp(global_logPosterior+log(priorPopTree));
      logPosterior = global_logPosterior +log(priorPopTree);
      ESS = global_sumESS/coldCh.GetNumLoci();
      // ESS = coldCh.GetNumIterations()*exp(2*global_logPosterior - global_logPosteriorSquared);
	// ESS = n(\bar w)^2 / \bar w^2
      /*
      std::cout << "coldCh.GetNumIterations()=" <<coldCh.GetNumIterations()
	       <<" global_logPosterior = " << global_logPosterior
	       <<"  global_logPosteriorSquared="<< global_logPosteriorSquared<<"\n";
      */
      std::cout << "Average effective sample size (ESS) of the importance sampling is "<< ESS 
		<< "(<"<<coldCh.GetNumIterations() << ", "<<ESS/coldCh.GetNumIterations()*100 <<"%).\n";
      // std::cout << "Effective sample size of the importance sampling (product of sum, log) is "<< coldCh.GetNumLoci()*log(coldCh.GetNumIterations())+2*global_logPosterior - global_logPosteriorSquared <<".\n";
    }
  
  if(crr_procID ==0)
    {
      totalComputingTime_eigen += coldCh.get_eachComputingTime_eigen();
      totalComputingTime_condiProb += coldCh.get_eachComputingTime_condiProb();
      totalNum_eigenFunctionCalls++;
      totalNum_condiProbFunctionCalls += coldCh.get_countCondiProbFunctionCalls();
      // std::cout << "totalComputingTime_eigen = " << totalComputingTime_eigen.count()/1000 << "milliseconds\n";
    }
  
  // std::cout << "crr_procID=" << crr_procID << " logPosterior = " << logPosterior <<"\n";   
  
  return logPosterior;
}



// 1/4/2016
// 'demographicPara' is a 1x6 matrix.
long double MaxPosterior::computeLogJointDensity_MPI_overSubLoci(Eigen::MatrixXd demographicPara, IM im, popTree* poptree, Chain coldCh, unsigned int nProcs, unsigned int crr_procID)
{

  int valid = poptree->replacePara(demographicPara);
  if(valid == -1)
    return 1*numeric_limits<double>::infinity();
  else if(valid == 1)
    {

      coldCh.compute_partialJointPosteriorDensity_overSubLoci(poptree, im, crr_procID, nProcs);
      std::vector<long double> logExpectationOfEachCoalProb = coldCh.get_logExpectationOfCoalProb();
      unsigned int numSubLoci = coldCh.getNumSubLoci();
      
      long double posterior = 1;
      if(numSubLoci != logExpectationOfEachCoalProb.size())
	{
	  std::cout << "** Error In MaxPosterior::computeJointDensity_MPI_overSubSample() **\n";
	  std::cout << "numSubLoci = " << numSubLoci 
		    << " expectationOfEachCoalProb.size() = " << logExpectationOfEachCoalProb.size() <<"\n\n";
	}
      long double local_logPosterior_subLoci =0;
      long double global_logPosterior = 0;
      for(unsigned int lc=0; lc< numSubLoci; lc++)
	{
	  local_logPosterior_subLoci += logExpectationOfEachCoalProb.at(lc);
	  // std::cout<< "lc = " << lc <<" local_logPosterior_subLoci = "<< local_logPosterior_subLoci <<"\n";
	}
      // std::cout << "crr_procID=" << crr_procID << " local_logPosterior_subLoci = " << local_logPosterior_subLoci <<"\n";
      
      MPI::COMM_WORLD.Barrier();
      MPI::COMM_WORLD.Allreduce(&local_logPosterior_subLoci, &global_logPosterior, 1, MPI_LONG_DOUBLE, MPI_SUM);
      MPI::COMM_WORLD.Barrier();
      
      long double logPosterior = 0;
      if(crr_procID ==0)
	{	
	  Eigen::Vector3d paraMax = im.get_paraMax();
	  double priorPopTree = poptree->computeJointPrior(paraMax);
	  //posterior = exp(global_logPosterior+log(priorPopTree));
	  logPosterior = global_logPosterior; //  +log(priorPopTree);
	  // std::cout <<"log(priorPopTree) = " << log(priorPopTree) <<"\n";
	} 
      return logPosterior;
    }
}



// 10/27/2015
// 'demographicPara' is a 1x6 matrix.
long double MaxPosterior::computeLogJointDensity_mutationScalars_MPI_overSubLoci(Eigen::MatrixXd demographicPara, IM im, popTree* poptree, Chain coldCh, unsigned int nProcs, unsigned int crr_procID)
{

  int valid = poptree->replacePara(demographicPara);
  if(valid == -1)
    return 1*numeric_limits<double>::infinity();
  else if(valid == 1)
    {
      poptree->replacePara(demographicPara);
      
      coldCh.compute_partialJointPosteriorDensity_mutationScalars_overSubLoci(poptree, im, crr_procID, nProcs);
      
      unsigned int nSample = coldCh.get_nSubSample();
      std::vector<long double> fullLogJointPosterior_local = coldCh.get_partialLogJointCoalProb();
      std::vector<long double> fullLogJointPosterior_global;
      fullLogJointPosterior_global.resize(nSample);
      
      // std::vector<double> expectationOfEachCoalProb = coldCh.get_expectationOfCoalProb();
      // unsigned int numSubLoci = coldCh.getNumSubLoci();
      
      // double posterior = 1;
      if(nSample != fullLogJointPosterior_local.size())
	{
	  std::cout << "** Error In MaxPosterior::computeLogJointDensity_mutationScalars_MPI_overSubLoci() **\n";
	  std::cout << "nSample = " << nSample 
		    << " fullLogJointPosterior_local.size() = " << fullLogJointPosterior_local.size() <<"\n\n";
	}
      
      long double local_logPosterior =0;
      long double global_logPosterior = 0;
      long double sumJointPosterior = 0.0;
      long double logSumJointPosterior = 0.0;
      long double max_logPosterior = -1*numeric_limits<long double>::infinity();
      for(unsigned int ii=0; ii< nSample; ii++)
	{
	  local_logPosterior = fullLogJointPosterior_local.at(ii);
	  global_logPosterior = 0;
	  MPI::COMM_WORLD.Barrier();
	  MPI::COMM_WORLD.Allreduce(&local_logPosterior, &global_logPosterior, 1, MPI_LONG_DOUBLE, MPI_SUM);
	  MPI::COMM_WORLD.Barrier();
	  
	  if(crr_procID==0)
	    {
	      Eigen::Vector3d paraMax = im.get_paraMax();
	      double priorPopTree = poptree->computeJointPrior(paraMax);
	      fullLogJointPosterior_global.at(ii) = global_logPosterior;
	      
	      sumJointPosterior += exp(global_logPosterior);
	      
	      if(max_logPosterior < global_logPosterior)
		max_logPosterior = global_logPosterior;
	    }
	}

      
      if(crr_procID==0)
	{
	  long double ratios = 0.0;
	  for(unsigned int ii=0; ii< nSample; ii++)
	    {
	      ratios += exp(fullLogJointPosterior_global.at(ii)-max_logPosterior);
	    }      
	  logSumJointPosterior = max_logPosterior + log(ratios);
	}
      
      long double logPosterior = 0;
      if(crr_procID ==0)
	{	
	  Eigen::Vector3d paraMax = im.get_paraMax();
	  double priorPopTree = poptree->computeJointPrior(paraMax);
	  //posterior = exp(global_logPosterior+log(priorPopTree));
	  logPosterior = logSumJointPosterior + log(priorPopTree)-log(nSample);
	}
      
      if(crr_procID ==0)
	{
	  totalComputingTime_eigen += coldCh.get_eachComputingTime_eigen();
	  totalComputingTime_condiProb += coldCh.get_eachComputingTime_condiProb();
	  totalNum_eigenFunctionCalls++;
	  totalNum_condiProbFunctionCalls += coldCh.get_countCondiProbFunctionCalls();
	  // std::cout << "totalComputingTime_eigen = " << totalComputingTime_eigen.count()/1000 << "milliseconds\n";
	}
      
      
      return logPosterior;
    }
}


