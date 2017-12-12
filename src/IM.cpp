/*MIST 2016 Yujin Chung and Jody Hey */

/*
 * IM.cpp
 *
 */


#include <iostream>
#include <fstream>
#include "IM.hpp"




void IM::initialize_MCMCsetting(int nChains, int nTreeSample, int samplingFromPriors)
{
	// FIXME
	// YC 3/4/2014
	// Later this function should be modified
	// so that the variables below will be
	// initialized by the user's input.
	n_chains = nChains;
	//n_Burning = nBurning;
	//thinning = nthinning;
	//n_MCMCgen = n_Burning+(n_sampledTrees-1)*thinning+1; //nMCMCgen;
	n_sampledTrees = nTreeSample;
	samplingFromPriorOnly = samplingFromPriors;
}

 // If the returning value is 0, then stop the software. If 1, execute MCMC.
unsigned int IM::initialization(int argc, char *argv[], unsigned int processID)
{   
  // Defaults
  splittingTimeMax = pow(10,6);
  popSizeMax = pow(10,6);
  migRateMax = 0;
  nPairsMut =0;
  nParaVectors = 100;
  checkpoint = 0;
  n_chains = 1;
  n_Burning =100;
  thinning = 10;
  n_sampledTrees=100;
  priorType=0; // improper prior
  lociInParallel = 1;
  samplingFromPriorOnly=0;
  popSizeMax = 10;
  multiLocusSpecific_mutationRate= 1; // constant per-site mutation rate for M mode
  nPairsMut =0;
  printFreq =10;
  TreesWithMaxP =0;
  newickTreeFormat = 0;
  samePopulationSizes = 0;
  sameMigrationRates = 0;
  checkpoint = 0;
  ancPop=1;

  locus locus_tmp;
  locus_tmp.initialize_multiLocusSpecific_mutationRate(multiLocusSpecific_mutationRate);
  loci.push_back(locus_tmp);

  if(argc <= 4){
    if(processID==0)
      {
	std::cout <<"This program is run through a command line interface.\n"
		  <<"To execute the program, please read the manual.\n";
      }
    return 0;

  }else{  
    if(processID == 0)
      {
	std::cout << "Command line: " ;
	for(unsigned int i=0; i< argc; i++)
	  {
	    std::cout << argv[i] <<" ";
	  }
	std::cout << "\n\n";
      }    
    
    unsigned int counter = 1;
    while(counter < argc)
      {
	char *pstr = argv[counter];
	char flag = toupper(pstr[1]); // convert lowercase letter to uppercase
	switch(flag)
	  {
	  case 'R':
	    MLmodes = atoi(argv[counter+1]);
	    if(MLmodes>=2)
	      multiLocusSpecific_mutationRate = 0;
	    break;
	  case 'B':
	    n_Burning = atoi(argv[counter+1]);
	    break;
	  case 'N':
	    if(MLmodes ==1)
	      thinning = atoi(argv[counter+1]);
	    else
	      {
		n_loci = atoi(argv[counter+1]);
		counter++;
		unsigned int nGeneCopies = atoi(argv[counter+1]);	
		loci.at(0).set_nGeneCopies(nGeneCopies);
	      }
	    break;
	  case 'L':
	    if(MLmodes <= 1) // MCMC
	      n_sampledTrees = atoi(argv[counter+1]);
	    else if(MLmodes >=2) // L mode
	      {
		n_sampledTrees = atoi(argv[counter+1]);
		nTrees2compute = n_sampledTrees;
		/*
		counter++;
		nTrees2compute = atoi(argv[counter+1]);
		if(nTrees2compute > n_sampledTrees)
		  {
		    std::cout << "\n *** Error in IM::initialization() ***\n";
		    std::cout << "n_sampledTrees(=" << n_sampledTrees << ")"
			      << " is smaller than nTrees2compute(="<< nTrees2compute <<")\n\n";
		  }
		counter++;
		TreesWithMaxP = atoi(argv[counter+1]);
		counter++;
		lociInParallel = atoi(argv[counter+1]);
		*/
	      }
	    break;
	  case 'D':
	    samplingFromPriorOnly = atoi(argv[counter+1]);
	    break;
	      case 'C':
		if(MLmodes == 1) // MCMC
		  n_chains = atoi(argv[counter+1]);
		else if(MLmodes >= 2 || MLmodes <= 5) // reading the true parameters
		  {
		    samePopulationSizes = atoi(argv[counter+1]);
		    counter++;
		    sameMigrationRates = atoi(argv[counter+1]);
		  }
		break;
	      case 'G':
		if(MLmodes == 1) // MCMC
		  {
		    priorType = atoi(argv[counter+1]);
		    if(priorType==0)
		      changeOfRatePoint =0;
		    else if(priorType ==1)
		      changeOfRatePoint =1;
		    else if(priorType == 2)
		      {
			counter++;
			changeOfRatePoint = atof(argv[counter+1]);
		      }
		  }
		else if(MLmodes >= 2 )
		  priorType = atoi(argv[counter+1]);		
		else if(MLmodes == 5)
		  gridsize = atoi(argv[counter+1]);
		break;
	  case 'I':
		if(MLmodes==1) // MCMC
		  {
		    ifstream inFile;
		    string word;
		    unsigned int line_counter = 0;
		    unsigned int locus_counter = 0;
		    loci.resize(0);
		    
		    inFile.open(argv[counter+1]);
		    while(!inFile.eof() && (locus_counter == 0 || locus_counter < n_loci))
		      {
			locus locus_tmp;
			locus_tmp.initialize_multiLocusSpecific_mutationRate(multiLocusSpecific_mutationRate);
			if(line_counter == 0)
			  {
			    // inFile >> poptree_string;
			    unsigned int nCopies;
			    inFile >> nCopies;
			    locus_tmp.set_nGeneCopies(nCopies);
			    inFile >> n_loci;
			    if(processID == 0)
			      {
				// std::cout << "Population tree is" << poptree_string <<".\n";
				std::cout << "Reading the data file, '" << argv[counter+1]<<"'...\n";
				std::cout << nCopies <<" gene copies and " << n_loci <<" loci.\n";
			      }
			    line_counter +=2;
			  }
			else
			  {
			    unsigned int nSites;
			    inFile >> nSites;
			    locus_tmp.set_nSites(nSites);
			    inFile >> word;
			    locus_tmp.initializeLikelihoodModel(word);
			    line_counter++;
			    vector<vector<char> > seq; /// The sequence of a locus, containing
			    /// the sequences of "n_genes" genes.
			    /// 2-d array: "n_genes" by "n_sites"
			    unsigned int nGeneCopies = locus_tmp.get_nGeneCopies();
			    for(unsigned int i=0; i<nGeneCopies; i++)
			      {		
				inFile >> word;
				locus_tmp.popNames.push_back(word);
				inFile >> word;
				locus_tmp.tipNames.push_back(word);
				inFile >> word;
				vector<char> row(word.c_str(), word.c_str() + word.size() + 1u);
				seq.push_back(row);
			      }
			    line_counter +=  nGeneCopies; //locus_tmp.n_geneCopies;
			    locus_tmp.compute_pi_uniqSeq(seq);
			    loci.push_back(locus_tmp);
			    locus_counter++;
			  }
		      }
		    inFile.close();
		  }
		else
		  nParaVectors = atoi(argv[counter+1]);
		break;
	      case 'J':
		nJointPoints = atoi(argv[counter+1]);
		break;
	      case 'M': 
	      if(MLmodes == 1) // MCMC  Default: 1	
		  lociInParallel=atof(argv[counter+1]);
	      else
		migRateMax  = atof(argv[counter+1]); // upper bound of migration rate
	      break;
	  case 'O':
	    optimizer = atoi(argv[counter+1]);
	    break;
	  case 'Q':
	      popSizeMax = atof(argv[counter+1]); // upper bound of population size
	      break;
	    case 'T':
	      splittingTimeMax = atof(argv[counter+1]); // upper bound of splitting time
	      break;
	    case 'S':
	      if(MLmodes == 1) // MCMC
		{
		  nPairsMut = atof(argv[counter+1]); 
		}
	      else if(MLmodes >=2) // L mode
		{
		  newickTreeFormat = atoi(argv[counter+1]); // 1 if the input gene trees for L mode are in newick format
		  if(newickTreeFormat == 1)
		    {
		      counter++;
		      newickTreeFileName = argv[counter+1];
		    }
		}
	      break;
	    case 'A':
	      ancPop = atof(argv[counter+1]); // 0 if island model; 1 if isolation model
	      break;
	    case 'P':
	      if(MLmodes == 1) // MCMC
		{
		  printFreq = atof(argv[counter+1]); // print frequency  
		}
	      else if(MLmodes >=2 && MLmodes <= 5) // reading the true parameters
		{
		  unsigned int npara = atoi(argv[counter+1]);
		  truePara.resize(1,npara);
		  for(unsigned int i=0; i<npara; i++)
		    {
		      double para = atof(argv[counter+2+i]);
		      truePara(0,i) = para;
		    }
		  counter += npara;
		}
	      break;
	    case 'U':
	      multiLocusSpecific_mutationRate= atof(argv[counter+1]);
	      //  loci.at(0).initialize_multiLocusSpecific_mutationRate(multiLocusSpecific_mutationRate);
	      break;
	    case 'H':
	      poptree_string = std::string(argv[counter+1]);
	      break;	 
	    case 'W':
	      checkpoint = atoi(argv[counter+1]);
	      if(checkpoint ==1 || checkpoint ==3)
		{
		  counter++;
		  howOften_checkpoint = atoi(argv[counter+1]);
		}
	      break;
	    default:
	      std::cout << "Error in command line: flag "<< flag <<" not found\n\n";
	    }
	  
	  counter+=2;
	}

  loci.at(0).initialize_multiLocusSpecific_mutationRate(multiLocusSpecific_mutationRate);

  if(ancPop == 0) // if island model (no ancestral populations)
    splittingTimeMax = std::numeric_limits<double>::infinity();

  Eigen::Vector3d paraMax(popSizeMax, splittingTimeMax, migRateMax);
  
  
  return 1;
  }
}

void IM::print_data_pi()
{
  //vector<locus>::iterator iter;
  for(unsigned int i = 0; i < n_loci; i++)
    {
      cout << "Locus " << i <<":\n";
      loci.at(i).print_likelihoodModel();
      loci.at(i).print_data();
      if(loci.at(i).getLikelihoodModel() == 2) // HKY
	{	  
	  loci.at(i).print_pi();
	}
    }
}


int IM::get_nChains()
{
	return n_chains;
}

int IM::get_nMCMCgen()
{
	return n_MCMCgen;
}

Eigen::Vector3d IM::get_paraMax()
{
	Eigen::Vector3d paraMax(popSizeMax,splittingTimeMax,migRateMax);
	return paraMax;
}
