/*MIST 2016 Yujin Chung and Jody Hey */


#include <iostream>
#include <vector>
#include <fstream>
#include "Eigen/Dense"
#include "IM.hpp"
#include "popTree.hpp"
#include "misc.hpp"


void locus::print_likelihoodModel()
{
  if(modelType == 0)
    std::cout << "Infinite-site model\n";
  else if(modelType ==1)
    std::cout << "JC model\n";
  else if(modelType ==2)
    std::cout << "HKY model\n";
  return;
}

void locus::initializeLikelihoodModel(string model)
{
  char flag = toupper(model[0]);
  switch(flag)
    {
    case 'I': // inifinite-site model
      modelType = 0;
      done_computingSumOfMuationsOnBranches = 0;
      break;
    case 'J':
      modelType = 1;
      break;
    case 'H':
      modelType = 2;
    }
  return;
}

unsigned int locus::get_n_sites_uniq()
{
	return n_sites_uniq;
}


/*
 * Explore 'seq_uniq', find the same site with the argument 'site'
 * and return the index of the site in 'seq_uniq'.
 */
int locus::IDsameSeq(vector<unsigned int> site)
{
  int ID = -1;
  int i=0;
  while(ID ==-1 && i<n_sites_uniq)
    {
      if(seq_uniq[i] == site)
	ID = i;
      else
	i++;
    }
  return ID;
}


/*
 * Computing the empirical distribution of nucleotides (pi)
 * and save the unique sequences with their frequencies.
 *
 * @param data The nucleotides read from the data file.
 */
void locus::compute_pi_uniqSeq_HKY(vector<vector<char> > data)
{
  int i,ID_sameSeq;
  char base;

  std::valarray<double> pi_temp (0.0,4);
	
	for(i=0; i<n_sites; i++)
	{
		vector<unsigned int> site;
		site.resize(n_geneCopies);
		for(unsigned int j=0; j<n_geneCopies; j++)
		{
			base = data[j][i];
	        if (base == 'a' || base == 'A')
	        {
	            pi_temp[0]++;
	            site[j] = 0;
	        }
	        else if (base == 'c' || base == 'C')
	        {
	            pi_temp[1]++;
	            site[j] = 1;
	        }
	        else if (base == 'g' || base == 'G')
	        {
	            pi_temp[2]++;
	            site[j] = 2;
	        }
	        else if (base == 't' || base == 'u' || base == 'T' || base == 'U')
	        {
	            pi_temp[3]++;
	            site[j] = 3;
	        }
	        // FIXME
	        /*
	        else if (c == 'n' || c == '-' || c == 'N' || c == '.')
	                  {
	                    L[li].seq[i][j] = -1;
	                  }
	                  else
	                  {
	                    IM_err(IMERR_DATAERROR,"BAD BASE in locus %d species %i base %i: %c",li, i + 1, j + 1, c);
	                  }
	                  */

		}
		if(i==0)
		{
			seq_uniq.push_back(site);
			n_sites_uniq = 1;
			freq_uniqueSeq.push_back(1);
		}
		else
		{
			ID_sameSeq = IDsameSeq(site);
			if(ID_sameSeq < 0)
			{
				seq_uniq.push_back(site);
				n_sites_uniq++;
				freq_uniqueSeq.push_back(1);
			}
			else
				freq_uniqueSeq.at(ID_sameSeq)++;
		}
	}

	double normalizer = pi_temp.sum();
	pi_temp = pi_temp/normalizer;
	pi.resize(4);
	pi = pi_temp;
}

/*
 * Computing the empirical distribution of nucleotides (pi)
 * and save the unique sequences with their frequencies.
 *
 * @param data The nucleotides read from the data file.
 */
void locus::compute_uniqSegSites_IS(vector<vector<char> > data)
{
  // std::cout << "In locus::compute_uniqSegSites_IS\n";

  int ID_sameSeq; // can be negative
  char ancestralBase, derivedBase;
  unsigned int flag_violateIS = 0;
  unsigned int found_derivedBase = 0;
  nSNPs = 0;
  n_sites_uniq = 0;
  vector<unsigned int> site;

  // std::cout << "n_sites = " << n_sites <<"\n";

  for(unsigned int i=0; i<n_sites; i++)
    {
      found_derivedBase =0;
      flag_violateIS = 0;
      site.resize(n_geneCopies);
      unsigned int j=0;
      while(j<n_geneCopies && flag_violateIS ==0)
	{
	  if(j==0)
	    {
	      ancestralBase = data[j][i];
	      site[j] = 0;
	    }
	  else 
	    {
	      // if(strcmp(data[j][i], ancestralBase) == 0) // same base
	      if(data[j][i] == ancestralBase) // same base
		site[j] = 0;
	      else
		{
		  if(found_derivedBase == 0)
		    {
		      derivedBase = data[j][i];
		      site[j] = 1;
		      found_derivedBase = 1;
		    }
		  else if(data[j][i] == derivedBase) // same base
		    {
		      site[j] = 1;
		    }
		  else // different from ancestral and derived bases.
		    {
		      flag_violateIS = 1;
		    }
		}
	      
	    }
	  /*
	  std::cout << "site= " <<i << "genecopy = " << j <<"\n"
		    << "data[j][i]= " <<data[j][i] 
		    << " site[j] = " << site[j]  <<"\n";
	  */
	  j++;
	} // End of while(j<n_geneCopies && flag_violateIS ==0) 

      if(flag_violateIS == 0)
	{
	  if(found_derivedBase == 1)
	    {
	      nSNPs++;
	      if(seq_uniq.size() ==0)
		n_sites_uniq = 1;
	      else
		n_sites_uniq++;
	      seq_uniq.push_back(site);
	      freq_uniqueSeq.push_back(1);
	    }
	}
    }// END of for(unsigned int i=0; i<n_sites; i++)

  if(n_sites_uniq==0)
    {
      n_sites_uniq = 1;
      seq_uniq.push_back(site);
      freq_uniqueSeq.push_back(1);
    }
  return;
}

void locus::compute_pi_uniqSeq(vector<vector<char> > data)
{
  if(modelType == 0)
    compute_uniqSegSites_IS(data); // Actually not uniq site. All segregating sites.
  else if(modelType == 2)
    compute_pi_uniqSeq_HKY(data);
  return;
}


void locus::print_pi()
{
	char bases[]={'A','C','G','T'};
	for(int i=0;i<4; i++)
	{
		cout << bases[i] <<": " << pi[i] <<", ";
	}
	cout <<"\n";
}

void locus::print_data()
{
  if(modelType == 0) // IS
    cout << n_geneCopies <<" genes and " << n_sites <<" sites (" << nSNPs<< " SNPs)\n";
  else 
    cout << n_geneCopies <<" genes and " << n_sites <<" sites (" <<n_sites_uniq <<" unique sites)\n";
  cout << "Population\t Individual\n";
  for(unsigned int i=0;i<n_geneCopies;i++)
    {
      cout << popNames.at(i) <<"\t " <<tipNames.at(i) << "\n";
    }
  /*
    cout << "Population\t Individual\t Sequence\n";
    for(unsigned int i=0;i<n_geneCopies;i++)
    {
    cout << popNames.at(i) <<"\t " <<tipNames.at(i) << "\t ";
    for(int j=0;j< n_sites_uniq;j++)
    cout << (seq_uniq.at(j)).at(i) <<"";
    cout <<"\n";
    }
    cout << "frequencies:\n";
    for(int i=0; i<n_sites_uniq; i++)
    cout << "\t site " << i+1 << ": " << freq_uniqueSeq.at(i) <<"\n";
	*/
  cout << "\n";
}

void IM::initialize_nloci(int nloci)
{
	n_loci = nloci;
}


vector<int> popTree::find_popIDs2migrate_moveDown(vector<int> popIDs2migrate, double lowerB, double upperB)
{
	// REMOVE
	// std::cout << "In popTree::find_popIDs2migrate_moveDown().\n";
	// std::cout << "popTree is \n";
	// print_poptree();
	// std::cout << "popIDs2migrate is..\n\t";
	// print_vectorInt(popIDs2migrate);


	if(age < upperB && par->age > lowerB)
		popIDs2migrate.push_back(popID);
	if(age > lowerB && !isTip)
	{
		vector<int> tmp1 = desc[0]->find_popIDs2migrate_moveDown(popIDs2migrate, lowerB, upperB);
		popIDs2migrate.resize(tmp1.size());
		popIDs2migrate = tmp1;

		vector<int> tmp2 = desc[1]->find_popIDs2migrate_moveDown(popIDs2migrate, lowerB, upperB);
		popIDs2migrate.resize(tmp2.size());
		popIDs2migrate = tmp2;
	}

	return popIDs2migrate;
}

vector<int> popTree::find_popIDs2migrate_moveUpDown(vector<int> popIDs2migrate, double lowerB, double upperB)
{
	if(par->desc[0]->popID != popID)
	{
		vector<int> tmp = desc[0]->find_popIDs2migrate_moveDown(popIDs2migrate,lowerB,upperB);
		popIDs2migrate.resize(tmp.size());
		popIDs2migrate = tmp;
	}
	else if(par->desc[1]->popID != popID)
	{
		vector<int> tmp = desc[1]->find_popIDs2migrate_moveDown(popIDs2migrate,lowerB,upperB);
		popIDs2migrate.resize(tmp.size());
		popIDs2migrate = tmp;
	}

	if(!(par->isRoot))
	{
		vector<int> tmp =par->find_popIDs2migrate_moveUpDown(popIDs2migrate,lowerB,upperB);
		popIDs2migrate.resize(tmp.size());
		popIDs2migrate = tmp;
	}

	return popIDs2migrate;
}

/***
 * Note that the given population tree should NOT be the root node
 */
vector<int> popTree::find_popIDs2migrate()
{
	// REMOVE
	// std::cout << "in popTree::find_popIDs2migrate()\n";
	// print_poptree();

	vector<int> popIDs2migrate;

	if(par->desc[0]->popID != popID)
	{
		popIDs2migrate = par->desc[0]->find_popIDs2migrate_moveDown(popIDs2migrate,age,par->age);
	}
	else if(par->desc[1]->popID != popID)
	{
		popIDs2migrate = par->desc[1]->find_popIDs2migrate_moveDown(popIDs2migrate,age,par->age);
		// REMOVE
		//std::cout << "popIDs2migrate is..\n\t";
		//print_vectorInt(popIDs2migrate);
	}

	if(!(par->isRoot))
	{
		popIDs2migrate = par->find_popIDs2migrate_moveUpDown(popIDs2migrate,age,par->age);
	}

	return popIDs2migrate;
}

