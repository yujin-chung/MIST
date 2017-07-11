/* MIST copyright 2016 by Yujin Chung and Jody Hey */

#include <stdlib.h>
#include <iostream>
#include <vector>
#include <valarray>
#include <cmath>
#include "misc.hpp"
#include "coaltree.hpp"

using namespace std;



/**
 *  LOCAL FUNCTIONS
 */



/**
 * YC 1/17/2014
 *
 * Find the row and column indices of the minimum value of the given matrix 'mat'
 * @param mat Double array. The numbers of rows and columns should be the same.
 * @param size The number of rows (columns) of 'mat'
 */
vector<int> Matrix::get_nodeIDs_minMat ()
{
	vector<int> ids; // to save the indices (row and column) of the minimum value of the given matrix.
	ids.resize(2);

	if(nrow==1 && ncol==1)
		return ids;
	else
	{
	int random_n = 0; // random number
	int idx_min = 0;

	valarray<double> sub_mat = (mat[mat >= 0.0]);

	// return the IDs of the nodes to join
	valarray<int> index(mat.size());
	for(std::size_t i=0;i < mat.size();i++)
		index[i] = i;
	valarray<int> sub_index(index[mat==sub_mat.min()]);

	if(sub_index.size()>1) // if there are multiple minimum values
	{
		random_n = rand() % sub_index.size(); // Randomly pick a value from {0,...,mat.size()-1}
		idx_min = sub_index[random_n];
	}
	else // only one minimum value
		idx_min = sub_index[0];


	ids.at(0) = div(idx_min , ncol).quot;
	ids.at(1) = div(idx_min , ncol).rem+1;

	return ids;
	}
}

/*
 * private member 'mat' is reduced as a (nrow-1)x(ncol-1) matrix
 * The rows and columns (whose IDs are in 'nodeIDs') are removed
 * and new column and row are added.
 * Distances in the new column are the weighted average of the distances
 * in the removed rows and columns.
 * @param nodeIDs IDs of rows and columns to remove. nodeIDs[0] < nodeIDs[1].
 */
void Matrix::update_distMat_IS(vector<int> nodeIDs, vector<int> sizes)
{
  Matrix newMat;
  newMat.replace(mat,nrow,ncol);
  
  // Remove rows.
  if(nodeIDs.at(1) == nrow)
    newMat.del_row(nodeIDs.at(1)-1);
  else
    newMat.del_row(nodeIDs.at(1));
  newMat.del_row(nodeIDs.at(0));
  // Remove columns
  newMat.del_col(nodeIDs.at(1)-1);
  if(nodeIDs.at(0)==0)
    newMat.del_col(0);
  else
    newMat.del_col(nodeIDs.at(0)-1);
  
  // Declare new row and column to add
  valarray<double> row(-1.0,ncol-1),col(0.0,nrow-1);
  // Add a new row
  newMat.add_row(row);
  
  int w1, w2;
  double d1, d2;
  // double minDist=0.0;
  // Compute the weighted distances for the new column
  for(int i=0;i<nrow-1;i++)
    {      
      if(i<nodeIDs.at(0))
	{
	  // w1 = sizes[i]*sizes[nodeIDs.at(0)];
	  // w2 = sizes[i]*sizes[nodeIDs.at(1)];
	  d1 = val(i,nodeIDs.at(0)-1);
	  d2=  val(i,nodeIDs.at(1)-1);
	}
      else if(i+1<nodeIDs.at(1))
	{
	  // w1 = sizes[nodeIDs.at(0)]*sizes[i+1];
	  // w2 = sizes[i+1]*sizes[nodeIDs.at(1)];
	  d1 = val(nodeIDs.at(0),i);
	  d2 = val(i+1,nodeIDs.at(1)-1);
	}
      else
	{
	  // w1 =sizes[nodeIDs.at(0)]*sizes[i+2];
	  // w2 = sizes[nodeIDs.at(1)]*sizes[i+2];
	  d1 = val(nodeIDs.at(0),i+1);
	  d2 = val(nodeIDs.at(1),i+1);
	}     
      // col[i] = (w1*d1+w2*d2)/(w1+w2);
      col[i] = (d1+d2)/2;
      /*
      if(i==0)
	minDist = col[i];
      if(col[i] < minDist)
	minDist=col[i];
      */
    }
  //for(int i=0;i<nrow-1;i++)
  //  col[i] = col[i]- minDist;

  // Add a new column
  newMat.add_col(col);

  // Replace 'mat' by the new distance matrix
  replace(newMat);
}

/*
 * private member 'mat' is reduced as a (nrow-1)x(ncol-1) matrix
 * The rows and columns (whose IDs are in 'nodeIDs') are removed
 * and new column and row are added.
 * Distances in the new column are the weighted average of the distances
 * in the removed rows and columns.
 * @param nodeIDs IDs of rows and columns to remove. nodeIDs[0] < nodeIDs[1].
 */
void Matrix::update_distMat(vector<int> nodeIDs, vector<int> sizes)
{
  Matrix newMat;
  newMat.replace(mat,nrow,ncol);
  
  // Remove rows.
  if(nodeIDs.at(1) == nrow)
    newMat.del_row(nodeIDs.at(1)-1);
  else
    newMat.del_row(nodeIDs.at(1));
  newMat.del_row(nodeIDs.at(0));
  // Remove columns
  newMat.del_col(nodeIDs.at(1)-1);
  if(nodeIDs.at(0)==0)
    newMat.del_col(0);
  else
    newMat.del_col(nodeIDs.at(0)-1);
  
  // Declare new row and column to add
  valarray<double> row(-1.0,ncol-1),col(0.0,nrow-1);
  // Add a new row
  newMat.add_row(row);
  
  int w1, w2;
  double d1, d2;
  // Compute the weighted distances for the new column
  for(int i=0;i<nrow-1;i++)
    {
      if(i<nodeIDs.at(0))
	{
	  w1 = sizes[i]*sizes[nodeIDs.at(0)];
	  w2 = sizes[i]*sizes[nodeIDs.at(1)];
	  d1 = val(i,nodeIDs.at(0)-1);
	  d2=  val(i,nodeIDs.at(1)-1);
	}
      else if(i+1<nodeIDs.at(1))
	{
	  w1 = sizes[nodeIDs.at(0)]*sizes[i+1];
	  w2 = sizes[i+1]*sizes[nodeIDs.at(1)];
	  d1 = val(nodeIDs.at(0),i);
	  d2 = val(i+1,nodeIDs.at(1)-1);
	}
      else
	{
	  w1 =sizes[nodeIDs.at(0)]*sizes[i+2];
	  w2 = sizes[nodeIDs.at(1)]*sizes[i+2];
	  d1 = val(nodeIDs.at(0),i+1);
	  d2 = val(nodeIDs.at(1),i+1);
	}
      col[i] = (w1*d1+w2*d2)/(w1+w2);
      
    }
  // Add a new column
  newMat.add_col(col);

  // Replace 'mat' by the new distance matrix
  replace(newMat);
}

/**
 * YC 1/17/2014
 *
 * Computing the pairwise distance matrix for a locus. This distance
 * matrix is used when UPGMA tree (as an initial tree in MCMC) is
 * constructed. The distance here is the number of sites with
 * different nucleotides of two sequences.
 *
 * Note that the distance matrix 'distmat' is (#genes-1) by (#genes-1) matrix,
 * where 'distmat[i][j]' is the distance of (i+1)th and (j+2)th sequences,
 * for i,j=0,...,(#genes-2). Then,
 * (1) distmat[i][i-1] = 0 for i>=1;
 * (2) ddistmat[i][j] = distmat[j+1][i-1] for i>=1.
 * Therefore, we need the upper triangular part of the matrix 'distmat',
 * i.e., d[i][j] for j is greater or equal to i.
 *
 * @param index_locus The index of the locus for which distance matrix is computed.
 */
Matrix locus::get_distMat()
{
  Matrix distMat;
  distMat.initialize(n_geneCopies-1,n_geneCopies-1);
  double dist = 0.0;
  
  // Computing pairwise distances of sequences
  for (unsigned int i = 0; i < n_geneCopies -1; i++)
    {
      for (unsigned int j = 0; j < n_geneCopies-1; j++)
	{
	  if(j>=i)
	    {
	      dist = 0.0;
	      // Computing the distance between (i+1)th and (j+2)th sequences
	      for (int k = 0; k < n_sites_uniq; k++)
		dist += (seq_uniq.at(k).at(j+1) != seq_uniq.at(k).at(i)) * freq_uniqueSeq.at(k);
	      distMat.replace(dist,i,j);
	    }
	  else
	    {
	      distMat.replace(-1,i,j);
	    }
	  
	}
    }
  
  return distMat;
}




/**
 * YC 1/17/2014
 *
 * Build the UPGMA tree for a locus. The UPGMA tree is used
 * as an initial tree in MCMC.
 */

node* get_UPGMA (locus lc)
{
	// REMOVE
  // std::cout << "\nIn get_UPGMA()\n";

  vector<node*> list_nodes; // Nodes to join
  vector<node*>::iterator iter_node; // iterator for 'list_nodes';
  vector<int> list_sizes; // sizes of nodes (clade sizes) in 'list_nodes'
  vector<int>::iterator iter_size; // iterator for 'list_nodes';
  node *new_node = 0; // newly joined node
  vector<int> nodeIndex; // node ID's to merge.
  nodeIndex.resize(2);
  unsigned int nGeneCopies = lc.get_nGeneCopies();
  unsigned int n_nodes = nGeneCopies; // the number of nodes to join

  Matrix distmat;
  distmat.initialize(nGeneCopies-1,nGeneCopies-1); // Declare a distance matrix

  // Computing a distance matrix
  distmat.replace(lc.get_distMat());
  int i = 0; // indices for loop statements
  
  // Initialize the list of node and the list of node-sizes.
  for(i=0; i<n_nodes; i++)
    {
      //FIXME
      new_node = new node;
      new_node->initialization();
      //node tr;
      //new_node = &tr;
      new_node->isTip = 1; // The initial nodes to join are tips.
      new_node->set_isRoot(0);
      new_node->age = 0.0;
      new_node->tipID = i+1;
      new_node->isLikelihoodNULL =1;
      new_node->assign_siblingOrder(2);
      list_nodes.push_back(new_node);
      list_sizes.push_back(new_node->size_tree());
      
      //REMOVE
      //std::cout << "size of new node is " << new_node->size_tree() <<"\n";
      
    }

  // REMOVE
  for(unsigned int i=0; i< (unsigned) n_nodes; i++)
    list_nodes.at(i)->print_coaltree();

  // ---- Build UPGMA tree ---- //
  while(n_nodes>1)
    {
      if(n_nodes==2)
	{
	  nodeIndex[0] = 0; nodeIndex[1] = 1;
	}
      else
	{
	  // Get the node ID's with the minimum value of 'distmat'
	  nodeIndex = distmat.get_nodeIDs_minMat();
	  
	  // REMOVE
	  /*
	  for(unsigned int i=0; i< (unsigned) n_nodes; i++)
	    list_nodes.at(i)->print_coaltree();
	  distmat.print();
	  std::cout << "nodeIndex[0]=" << nodeIndex[0] << " nodeIndex[1]="<< nodeIndex[1] <<"\n";
	  */
	  // update the distance matrix 'distmat'
	  if(n_nodes>=3)
	    {
	      distmat.update_distMat(nodeIndex,list_sizes);
	    }
	}
      // Join the nodes with the minimum distance
      //FIXME
      new_node = new node;
      new_node->initialization();
      //node tmp;
      //new_node = &tmp;
      new_node->isTip = 0;
      new_node->set_isRoot(0);
      new_node->isLikelihoodNULL = 1;
      new_node->assign_siblingOrder(2);
      list_nodes[nodeIndex[0]]->assign_siblingOrder(0);
      list_nodes[nodeIndex[1]]->assign_siblingOrder(1);
      list_nodes[nodeIndex[0]]->par = new_node;
      list_nodes[nodeIndex[1]]->par = new_node;
      
      new_node->desc[0]=list_nodes[nodeIndex[0]];
      new_node->desc[1]=list_nodes[nodeIndex[1]];
      
      // Branch length is randomly generated from Exponential(1)
      // FIXME: YC 1/23/2014
      // the rate for Exponential distribution can be replaced
      // by the initial coalescent rate as done in IMa2.
      new_node->age = rexpdist(1)
	+ max(list_nodes[nodeIndex[0]]->age,list_nodes[nodeIndex[1]]->age);
      
      
      // update the list of nodes to join
      iter_node = list_nodes.begin()+nodeIndex[1];
      list_nodes.erase(iter_node);
      iter_node = list_nodes.begin()+nodeIndex[0];
      list_nodes.erase(iter_node);
      list_nodes.push_back(new_node);
      
      iter_size = list_sizes.begin()+nodeIndex[1];
      list_sizes.erase(iter_size);
      iter_size = list_sizes.begin()+nodeIndex[0];
      list_sizes.erase(iter_size);
      list_sizes.push_back(new_node->size_tree());
      
      //REMOVE
      //std::cout << "size of new node is " << new_node->size_tree() <<"\n";
      
      n_nodes = list_nodes.size();
      
      
    }
  
  list_nodes[0]->set_isRoot(1);
  return list_nodes[0];
}

/**
 * YC 1/17/2014
 *
 * Build the UPGMA tree for a locus. The UPGMA tree is used
 * as an initial tree in MCMC.
 */
/*
node* get_initialTree_IS (locus lc)
{
	// REMOVE
  // std::cout << "\nIn get_UPGMA_IS()\n";


  locus lc_redu = lc;
  vector<node*> list_nodes; // Nodes to join
  std::vector< std::vector<unsigned int>> list_nodeIDs;
  vector<node*>::iterator iter_node; // iterator for 'list_nodes';
  vector<int> list_sizes; // sizes of nodes (clade sizes) in 'list_nodes'
  vector<int>::iterator iter_size; // iterator for 'list_nodes';
  node *new_node = 0; // newly joined node
  vector<int> nodeIndex; // node ID's to merge.
  nodeIndex.resize(2);
  unsigned int nGeneCopies = lc.get_nGeneCopies();
  unsigned int n_nodes = nGeneCopies; // the number of nodes to join
  unsigned int nSites= lc.get_n_sites_uniq();
  Matrix distmat;
  distmat.initialize(nGeneCopies-1,nGeneCopies-1); // Declare a distance matrix

  // Computing a distance matrix
  distmat.replace(lc.get_distMat());
  
  // Initialize the list of node and the list of node-sizes.
  for(unsigned int i=0; i<n_nodes; i++)
    {
      //FIXME
      new_node = new node;
      new_node->initialization();
      //node tr;
      //new_node = &tr;
      new_node->isTip = 1; // The initial nodes to join are tips.
      new_node->isRoot = 0;
      new_node->age = 0.0;
      new_node->tipID = i+1;
      new_node->isLikelihoodNULL =1;
      new_node->assign_siblingOrder(2);
      list_nodes.push_back(new_node);
      list_sizes.push_back(new_node->size_tree());
      std::vector<unsigned int> clade; clade.push_back(i);
      list_nodeIDs.push_back(clade);
      //REMOVE
      //std::cout << "size of new node is " << new_node->size_tree() <<"\n";
      
    }

  // REMOVE
  for(unsigned int i=0; i< (unsigned) n_nodes; i++)
    list_nodes.at(i)->print_coaltree();

  while(n_nodes>1)
    {
      vector<unsigned int> site;
      for(unsigned int s=0; s< nSites & n_nodes>1 ; s++)
	{
	  site = lc.getSite_uniqSeq(s);
	}
    }


  // ---- Build UPGMA tree ---- //
  while(n_nodes>1)
    {
      if(n_nodes==2)
	{
	  nodeIndex[0] = 0; nodeIndex[1] = 1;
	}
      else
	{
	  // Get the node ID's with the minimum value of 'distmat'
	  nodeIndex = distmat.get_nodeIDs_minMat();
	  
	  // REMOVE
	  for(unsigned int i=0; i< (unsigned) n_nodes; i++)
	    list_nodes.at(i)->print_coaltree();
	  distmat.print();
	  std::cout << "nodeIndex[0]=" << nodeIndex[0] << " nodeIndex[1]="<< nodeIndex[1] <<"\n";

	  // update the distance matrix 'distmat'
	  if(n_nodes>=3)
	    {
	      // locus lc_new = lc_redu;
	      unsigned int nSites= lc_redu.get_n_sites_uniq();
	      // unsigned int nSeq = lc_redu.get_nGeneCopies();
	      
	      std::vector<std::vector<unsigned int> > seq;
	      std::vector<unsigned int> site;
	      std::vector<unsigned int> new_site;
	      
	      for(unsigned int i=0; i<nSites; i++)
		{
		  site = lc_redu.getSite_uniqSeq(i);
		  // std::cout << "site.size() = " << site.size() <<"\n";
		  unsigned int targetAllele = site.at(nodeIndex[0]);
		  if(site.at(nodeIndex[1]) != targetAllele)
		    {
		      std::cout <<"\n Error in node* get_initialTree_IS()\n";
		      std::cout << "site.at(nodeIndex[0]) = " << site.at(nodeIndex[0])
				<<" site.at(nodeIndex[1]) = " << site.at(nodeIndex[1]) <<"\n";
		    }	
		  unsigned int noMut =0;
		  for(unsigned int j=0; j<site.size() & noMut>1 ; j++)
		    {
		      if(j != nodeIndex[0] & j!= nodeIndex[1])
			{
			  new_site.push_back(site.at(j));
			  if(site.at(j) == targetAllele)
			    noMut++;
			}
		    }
		  if(noMut > 0)
		    {
		      new_site.push_back(targetAllele);
		      seq.push_back(new_site);
		    }
		}
	      if(seq.size() ==0)
		seq.push_back(site);
	      lc_redu.set_nSites(seq.size());
	      lc_redu.set_n_sites_uniq(seq.size());
	      lc_redu.set_seq_uniq(seq);  
	      
	      distmat.replace(lc_redu.get_distMat());
	      distmat.update_distMat_IS(nodeIndex,list_sizes);
	    }
	}
      // Join the nodes with the minimum distance
      //FIXME
      new_node = new node;
      new_node->initialization();
      //node tmp;
      //new_node = &tmp;
      new_node->isTip = 0;
      new_node->isRoot = 0;
      new_node->isLikelihoodNULL = 1;
      new_node->assign_siblingOrder(2);
      list_nodes[nodeIndex[0]]->assign_siblingOrder(0);
      list_nodes[nodeIndex[1]]->assign_siblingOrder(1);
      list_nodes[nodeIndex[0]]->par = new_node;
      list_nodes[nodeIndex[1]]->par = new_node;
      
      new_node->desc[0]=list_nodes[nodeIndex[0]];
      new_node->desc[1]=list_nodes[nodeIndex[1]];
      
      // Branch length is randomly generated from Exponential(1)
      // FIXME: YC 1/23/2014
      // the rate for Exponential distribution can be replaced
      // by the initial coalescent rate as done in IMa2.
      new_node->age = rexpdist(1)
	+ max(list_nodes[nodeIndex[0]]->age,list_nodes[nodeIndex[1]]->age);
      
      
      // update the list of nodes to join
      iter_node = list_nodes.begin()+nodeIndex[1];
      list_nodes.erase(iter_node);
      iter_node = list_nodes.begin()+nodeIndex[0];
      list_nodes.erase(iter_node);
      list_nodes.push_back(new_node);
      
      iter_size = list_sizes.begin()+nodeIndex[1];
      list_sizes.erase(iter_size);
      iter_size = list_sizes.begin()+nodeIndex[0];
      list_sizes.erase(iter_size);
      list_sizes.push_back(new_node->size_tree());
      
      //REMOVE
      //std::cout << "size of new node is " << new_node->size_tree() <<"\n";
      
      n_nodes = list_nodes.size();
      
      
    }
  
  list_nodes[0]->isRoot = 1;
  return list_nodes[0];
}
*/

/**
 * YC 1/17/2014
 *
 * Build the UPGMA tree for a locus. The UPGMA tree is used
 * as an initial tree in MCMC.
 */

node* get_UPGMA_IS (locus lc)
{
	// REMOVE
  // std::cout << "\nIn get_UPGMA_IS()\n";

  vector<node*> list_nodes; // Nodes to join
  vector<node*>::iterator iter_node; // iterator for 'list_nodes';
  vector<int> list_sizes; // sizes of nodes (clade sizes) in 'list_nodes'
  vector<int>::iterator iter_size; // iterator for 'list_nodes';
  node *new_node = 0; // newly joined node
  vector<int> nodeIndex; // node ID's to merge.
  nodeIndex.resize(2);
  unsigned int nGeneCopies = lc.get_nGeneCopies();
  unsigned int n_nodes = nGeneCopies; // the number of nodes to join
  locus lc_redu = lc;
  lc_redu.set_n_sites_uniq(lc.get_n_sites_uniq());
  lc_redu.set_seq_uniq(lc.get_seq_uniq());  
  lc_redu.set_nGeneCopies(lc.get_nGeneCopies());
  lc_redu.set_freq_uniqueSeq(lc.get_freq_uniqueSeq());

  Matrix distmat;
  distmat.initialize(nGeneCopies-1,nGeneCopies-1); // Declare a distance matrix

  // Computing a distance matrix
  distmat.replace(lc.get_distMat());
  
  // Initialize the list of node and the list of node-sizes.
  for(unsigned int i=0; i<n_nodes; i++)
    {
      //FIXME
      new_node = new node;
      new_node->initialization();
      //node tr;
      //new_node = &tr;
      new_node->isTip = 1; // The initial nodes to join are tips.
      new_node->set_isRoot(0);
      new_node->age = 0.0;
      new_node->tipID = i+1;
      new_node->isLikelihoodNULL =1;
      new_node->assign_siblingOrder(2);
      list_nodes.push_back(new_node);
      list_sizes.push_back(new_node->size_tree());
      
      //REMOVE
      //std::cout << "size of new node is " << new_node->size_tree() <<"\n";
      
    }

  // REMOVE
  /*
  for(unsigned int i=0; i< (unsigned) n_nodes; i++)
    list_nodes.at(i)->print_coaltree();
  */

  // ---- Build UPGMA tree ---- //
  while(n_nodes>1)
    {
      if(n_nodes==2)
	{
	  nodeIndex[0] = 0; nodeIndex[1] = 1;
	}
      else
	{
	  // Get the node ID's with the minimum value of 'distmat'
	  nodeIndex = distmat.get_nodeIDs_minMat();
	  
	  // REMOVE
	  /*
	  for(unsigned int i=0; i< (unsigned) n_nodes; i++)
	    list_nodes.at(i)->print_coaltree();
	  distmat.print();
	  std::cout << "nodeIndex[0]=" << nodeIndex[0] << " nodeIndex[1]="<< nodeIndex[1] <<"\n";
	  */

	  // update the distance matrix 'distmat'
	  if(n_nodes>=3)
	    {
	      // locus lc_new = lc_redu;
	      unsigned int nSites= lc_redu.get_n_sites_uniq();
	      // unsigned int nSeq = lc_redu.get_nGeneCopies();
	      
	      std::vector<std::vector<unsigned int> > seq;
	      std::vector<unsigned int> site;
	      std::vector<unsigned int> new_site; new_site.resize(0);
	      unsigned int targetAllele =0;
	      // std::cout << "nSites = " << nSites <<"\n";
	      for(unsigned int s=0; s<nSites; s++)
		{
		  new_site.resize(0);
		  site = lc_redu.getSite_uniqSeq(s);
		  // std::cout << "site.size() = " << site.size() <<"\n";
		  targetAllele = site.at(nodeIndex[0]);
		  // std::cout << "targetAllele = " << targetAllele <<"\n";
		  if(site.at(nodeIndex[1]) != targetAllele)
		    {
		      std::cout <<"\n Error in node* get_initialTree_IS()\n";
		      std::cout << "site.at(nodeIndex[0]) = " << site.at(nodeIndex[0])
				<<" site.at(nodeIndex[1]) = " << site.at(nodeIndex[1]) <<"\n";
		    }	
		  unsigned int noMut =0;
		  for(unsigned int g=0; g< site.size(); g++)
		    {
		      // std::cout << "g= " << g << " site.at(g) = " << site.at(g) <<"\n";
		      if(g != nodeIndex[0] & g!= nodeIndex[1])
			{
			  new_site.push_back(site.at(g));
			  if(site.at(g) == targetAllele)
			    noMut++;
			}
		    }
		  // std::cout << "noMut = " << noMut <<"\n";
		  if(noMut > 0)
		    {
		      new_site.push_back(targetAllele);
		      seq.push_back(new_site);
		    }
		}
	      if(seq.size() ==0)
		{
		  new_site.push_back(targetAllele);
		  seq.push_back(new_site);
		}
	      lc_redu.set_nSites(seq.size());
	      lc_redu.set_n_sites_uniq(seq.size());
	      lc_redu.set_seq_uniq(seq);  
	      lc_redu.set_nGeneCopies(new_site.size());
	      // lc_redu.set_freq_uniqueSeq(lc.get_freq_uniqueSeq());

	      // std::cout << "new nsites = " << seq.size() << " newNGeneCopies = " << new_site.size() <<"\n";

	      distmat.replace(lc_redu.get_distMat());
	      // distmat.update_distMat_IS(nodeIndex,list_sizes);
	    }
	}
      // Join the nodes with the minimum distance
      //FIXME
      new_node = new node;
      new_node->initialization();
      //node tmp;
      //new_node = &tmp;
      new_node->isTip = 0;
      new_node->set_isRoot(0);
      new_node->isLikelihoodNULL = 1;
      new_node->assign_siblingOrder(2);
      list_nodes[nodeIndex[0]]->assign_siblingOrder(0);
      list_nodes[nodeIndex[1]]->assign_siblingOrder(1);
      list_nodes[nodeIndex[0]]->par = new_node;
      list_nodes[nodeIndex[1]]->par = new_node;
      
      new_node->desc[0]=list_nodes[nodeIndex[0]];
      new_node->desc[1]=list_nodes[nodeIndex[1]];
      
      // Branch length is randomly generated from Exponential(1)
      // FIXME: YC 1/23/2014
      // the rate for Exponential distribution can be replaced
      // by the initial coalescent rate as done in IMa2.
      if(lc.get_multiLocusSpecific_mutationRate() ==1)
	{	  
	  double nSegSites = (double) lc.get_n_sites_uniq();
	  double nSites = (double) lc.get_nSites();
	  double nSeq = (double) lc.get_nGeneCopies();
	  double invRate = std::min(20.0, std::max(0.01/(2*nSeq-2), nSegSites/nSites/(2*nSeq-2) ));
	  double coalT = abs(rNormalDistribution(0.0, invRate));
	  new_node->age = coalT + max(list_nodes[nodeIndex[0]]->age,list_nodes[nodeIndex[1]]->age);
	  // new_node->age = rexpdist(1/invRate)
	  //  + max(list_nodes[nodeIndex[0]]->age,list_nodes[nodeIndex[1]]->age);
	}
      else
	new_node->age = rexpdist(1)
	  + max(list_nodes[nodeIndex[0]]->age,list_nodes[nodeIndex[1]]->age);

      
      
      // update the list of nodes to join
      iter_node = list_nodes.begin()+nodeIndex[1];
      list_nodes.erase(iter_node);
      iter_node = list_nodes.begin()+nodeIndex[0];
      list_nodes.erase(iter_node);
      list_nodes.push_back(new_node);
      
      iter_size = list_sizes.begin()+nodeIndex[1];
      list_sizes.erase(iter_size);
      iter_size = list_sizes.begin()+nodeIndex[0];
      list_sizes.erase(iter_size);
      list_sizes.push_back(new_node->size_tree());
      
      //REMOVE
      //std::cout << "size of new node is " << new_node->size_tree() <<"\n";
      
      n_nodes = list_nodes.size();
      
      
    }
  
  list_nodes[0]->set_isRoot(1);
  return list_nodes[0];
}


/**
 *  GLOBAL FUNCTIONS
 */


/**
 * Initialize coalescent trees for a loci when HKY model (or
 * substitution model) is considered.
 * The initial tree is UPGMA tree
 */

void initial_coalTree_HKY(vector<node*> tr, vector<locus> loci)
{
	//tr.resize(loci.size());
	//std::cout << "loci size: "<<loci.size()<<"\n";
	//std::cout << "tree size: "<<tr.size()<<"\n";
	for(unsigned int i=0; i < loci.size(); i++)
	{
		// tr.at(i) = get_UPGMA(loci.at(i));
		//tr.at(i)// = new node;
		tr.at(i) = get_UPGMA(loci.at(i));
		//tr.at(i)->desc[0]->par = tr.at(i);
		//tr.at(i)->desc[1]->par = tr.at(i);
		tr.at(i)->print_coaltree();
	}
}


/**
 * Initialize coalescent trees for a loci when HKY model (or
 * substitution model) is considered.
 * The initial tree is UPGMA tree
 */
node* initial_coalTree_HKY(locus lc)
{
	//REMOVE
	//std::cout << "In initial_coalTree_HKY()\n";

	// node *new_node;
	//vector<node*> trees;
	//vector<node*>::size_type i;
	//for(i=0; i < loci.size(); i++)
	//{
		node *new_node;
		new_node = get_UPGMA(lc);
		new_node->assignPopulations2Tips(lc);
		//trees.push_back(new_node);
		//}

	//REMOVE
	//std::cout << "End of initial_coalTree_HKY()\n";

	//return trees;
		return new_node;

}

/// OLD version - YC 11/4/2014
/**
 * Initialize coalescent trees for a loci when HKY model (or
 * substitution model) is considered.
 * The initial tree is UPGMA tree
 */
vector<node*> initial_coalTree_HKY(vector<locus> loci)
{
	//REMOVE
	//std::cout << "In initial_coalTree_HKY()\n";

	// node *new_node;
	vector<node*> trees;
	vector<node*>::size_type i;
	for(i=0; i < loci.size(); i++)
	{
		node *new_node;
		new_node = get_UPGMA(loci[i]);
		new_node->assignPopulations2Tips(loci[i]);
		trees.push_back(new_node);
	}

	//REMOVE
	//std::cout << "End of initial_coalTree_HKY()\n";

	return trees;

}





// YC 3/13/2015
// Remove singletons from the sequences and get UPGMA tree.
/**
 * Initialize coalescent trees when the infinite
 * site model is considered.
 * The initial tree is ...
 */
node* initial_coalTree_IS(locus lc)
// void makeIS (int ci, int li, int nosimmigration)
{
  // std::cout << "\nIn initial_coalTree_IS()\n";

  locus lc_noSingleton = lc;
  unsigned int nSites= lc_noSingleton.get_n_sites_uniq();
  unsigned int nSeq = lc_noSingleton.get_nGeneCopies();
  //unsigned int nUniqSites= lc_noSingleton.get_n_sites_uniq();
  // Under IS model, n_sites = n_sites_uniq

  vector<vector<unsigned int> > seq;
  vector<unsigned int> site;

  for(unsigned int i=0; i<nSites; i++)
    {
      site = lc_noSingleton.getSite_uniqSeq(i);
      // std::cout << "site.size() = " << site.size() <<"\n";
      unsigned int nZero = 0;      
      for(unsigned int j=0; j<site.size(); j++)
	{
	  if(site.at(j) ==0)
	    nZero++;
	}
      if(nZero != 1 && nZero != nSeq-1)
	seq.push_back(site);
    }
  if(seq.size() ==0)
    seq.push_back(site);
  lc_noSingleton.set_nSites(seq.size());
  lc_noSingleton.set_n_sites_uniq(seq.size());
  lc_noSingleton.set_seq_uniq(seq);


  node *new_node;
  new_node = get_UPGMA_IS(lc_noSingleton);
  new_node->assignPopulations2Tips(lc_noSingleton);
  
  new_node->computeTotalLengthsOfLineages();

  //std::cout << "nSites = " << nSites <<"\n";
  //std::cout << "seq.size(nsites without singleton) = " << seq.size() <<"\n";
  // new_node->print_coaltree();

  return new_node;
 
}                               /* makeIS */





// FIXME - YC
// I didn't fully understand this function, but just follow how IMa2 
// creates a new tree.
/**
 * Initialize coalescent trees when the infinite
 * site model is considered.
 * The initial tree is ...
 */
node* initial_coalTree_IS_fromIMa2(locus lc)
// void makeIS (int ci, int li, int nosimmigration)
{
  //int c[2], newedge, *curid;
  double ptime;

  unsigned int nGenes = lc.get_nGeneCopies();
  unsigned int curgenes = nGenes;
  unsigned int cursites = lc.get_n_sites_uniq();
  node *gtree; // struct edge *gtree = C[ci]->G[li].gtree;
  Eigen::MatrixXi distMat;
  distMat.resize(nGenes,nGenes); // distmat = alloc2Dint (L[li].numgenes, L[li].numgenes);
  vector<int> singletons;
  singletons.resize(cursites);

  vector<vector<unsigned int> > Seqs;
  Seqs.resize(cursites);
  for(unsigned int i=0; i<cursites; i++)
    {
      Seqs.at(i) =lc.getSite_uniqSeq(i);
    }

  // tempseq = alloc2Dint (L[li].numgenes, L[li].numsites);
  //curid = static_cast<int *> (malloc ((L[li].numgenes) * (sizeof (int))));
  // singletons = static_cast<int *> (malloc ((L[li].numsites) * (sizeof (int))));
  // curgenes = L[li].numgenes;
  // cursites = L[li].numsites;
  //for (i = 0; i < L[li].numgenes; i++)
  //  curid[i] = i;

  node *new_node = 0; // newly joined node
  vector<node*> list_nodes; // Nodes to join
  vector<node*>::iterator iter_node; // iterator for 'list_nodes';
  // Initialize the list of node and the list of node-sizes.
  for(unsigned int i=0; i<nGenes; i++)
    {
      new_node = new node;
      new_node->initialization();
      //node tr;
      //new_node = &tr;
      new_node->isTip = 1; // The initial nodes to join are tips.
      new_node->set_isRoot(0);
      new_node->age = 0.0;
      new_node->tipID = i+1;
      new_node->isLikelihoodNULL =1;
      new_node->assign_siblingOrder(2);
      list_nodes.push_back(new_node);
      // list_sizes.push_back(new_node->size_tree());
      
      //REMOVE
      //std::cout << "size of new node is " << new_node->size_tree() <<"\n";
      
    }

  // C[ci]->G[li].mignum = 0;
  // for (i = 0; i < 2 * L[li].numgenes - 1; i++)
  // {
  //  gtree[i].up[0] = gtree[i].up[1] = gtree[i].down = -1; // Tips
  //}

  /*
  for (i = 0; i < L[li].numgenes; i++)
  {
    for (j = 0; j < L[li].numsites; j++)
      tempseq[i][j] = L[li].seq[i][j];
  }
  */


  vector<unsigned int> coal;
  coal.resize(2);

  ptime = 0.0;
  distMat.setZero();
  while (curgenes > 1)
    {
      //std::cout << "curgenes = " << curgenes <<"\n";

      unsigned int nPairsSameSeq = 0;
      // Build a distance matrix of the full data.
      for(unsigned int i = 0; i < curgenes; i++) 
	{
	  for (unsigned int j = i + 1; j < curgenes; j++) 
	    {
	      // distMat(i,j) = distmat[j][i] = 0;
	      for (unsigned int site = 0; site < cursites; site++)
		{
		  // if (tempseq[i][site] != tempseq[j][site])
		  if(Seqs.at(site).at(i) != Seqs.at(site).at(j))
		    {
		      distMat(i,j)++; // distmat[i][j]++;
		      distMat(i,j)++; // distmat[j][i]++;
		    }
		}
	      // if (distmat[i][j] == 0)
	      if (distMat(i,j) == 0)
		nPairsSameSeq++; // Count the number of pairs of identical sequences
	    }
	} // END of for(unsigned int i = 0; i < curgenes; i++) 

      //std::cout << "nPairsSameSeq = " << nPairsSameSeq <<"\n";

      // All the sequences are different
      if (nPairsSameSeq == 0) 
	{
	  // Identifying singletons
	  for(unsigned int i = 0; i < cursites; i++)
	    {
	      unsigned int nDerivedAlles = 0;
	      for (unsigned int j = 0; j < curgenes; j++)
		nDerivedAlles += Seqs.at(i).at(j);// tempseq[j][i];
	      if (nDerivedAlles == 1 || nDerivedAlles == (curgenes - 1))
		singletons.at(i) = 1;	      
	      else
		singletons.at(i) = 0;
	    }
	  // Re-build a distance matrix
	  nPairsSameSeq = 0;	  
	  for(unsigned int i = 0; i < curgenes; i++)
	    {
	      for (unsigned int j = i + 1; j < curgenes; j++)
		{
		  distMat.setZero();
		  // distmat[i][j] = distmat[j][i] = 0;
		  for (unsigned int site = 0; site < cursites; site++)
		    {
		      //if (tempseq[i][site] != tempseq[j][site] && singletons[site] == 0)
		      // If the site is not singleton, we count the number of differences
		      if(Seqs.at(site).at(i) != Seqs.at(site).at(j) && singletons.at(site) == 0)
			{
			  distMat(i,j)++; // distmat[i][j]++;
			  distMat(i,j)++; // distmat[j][i]++;
			}
		    }
		  if (distMat(i,j) == 0)
		    nPairsSameSeq++; //num++;
		}
	    }
	} // END of if (nPairsSameSeq == 0) 

      // std::cout << "nPairsSameSeq = " << nPairsSameSeq <<"\n";

      // Still all the sequences are different
      if (nPairsSameSeq == 0)
	{
	  std::cout << "\n*** Error in initial_coalTree_IS() ***\n";
	  printf ("Data not compatible with infinite sites model\n");
	  /*
	  for (unsigned int i = 0; i < curgenes; i++)
	    {
	      printf ("seq %i: ", curid[i]);
	      for (j = 0; j < cursites; j++)
		printf ("%i", tempseq[i][j]);
	      printf ("\n");
	    }
	  */
	  // IM_err(IMERR_INFINITESITESFAIL,"locus #: %d, name: %s,  sequence#: %d", li,L[li].name, curid[i]);
	  
	}

      // Determine which nodes will join
      unsigned int coalnum = runiform_discrete(nPairsSameSeq); // (int) (uniform () * nParisSampleSeq);
      unsigned int num = 0;
      for (unsigned int i = 0; i < curgenes; i++)
	{
	  for (unsigned int j = i + 1; j < curgenes; j++)
	    {
	      if (distMat(i,j) == 0)
		{
		  num++;
		  if (num > coalnum)
		    {
		      if (i < j)
			{
			  coal.at(0) = i;
			  coal.at(1) = j;
			}
		      else
			{
			  coal[0] = j;
			  coal[1] = i;
			}
		      i = curgenes;
		      j = curgenes;
		    }
		}
	    }
	}
      // REMOVE
      //std::cout << "coal.at(0)  = " << coal.at(0)<< " coal.at(1) = "<< coal.at(1) << "\n";
      if(cursites >= 2)
	{
	  for (unsigned int i = 0; i < cursites; i++)
	    {
	      //std::cout << "Seqs.at(i).at(coal.at(0)) = " << Seqs.at(i).at(coal.at(0))
	      //		<<  " Seqs.at(i).at(coal.at(1)) = " <<Seqs.at(i).at(coal.at(1)) <<"\n";
	      if (Seqs.at(i).at(coal.at(0)) != Seqs.at(i).at(coal.at(1)) )
		{
		  for (unsigned int j = 0; j < curgenes; j++)
		    {
		      for (unsigned int k = i; k < cursites - 1; k++)
			Seqs.at(k).at(j) = Seqs.at(k+1).at(j);
		    }
		  i--;
		  cursites--;
		}
	    }
	}

      // Join the nodes 
      new_node = new node;
      new_node->initialization();
      new_node->isTip = 0;
      new_node->set_isRoot(0);
      new_node->isLikelihoodNULL = 1;
      new_node->assign_siblingOrder(2);
      list_nodes[coal.at(0)]->assign_siblingOrder(0);
      list_nodes[coal.at(1)]->assign_siblingOrder(1);
      list_nodes[coal.at(0)]->par = new_node;
      list_nodes[coal.at(1)]->par = new_node;
      new_node->desc[0]=list_nodes[coal.at(0)];
      new_node->desc[1]=list_nodes[coal.at(1)];
      new_node->age = rexpdist(1)+ max(list_nodes[coal.at(0)]->age,list_nodes[coal.at(1)]->age);

      // update the list of nodes to join          
      iter_node = list_nodes.begin()+coal.at(1);
      list_nodes.erase(iter_node);
      iter_node = list_nodes.begin()+coal.at(0);
      list_nodes.erase(iter_node);
      list_nodes.push_back(new_node);


      //newedge = 2 * L[li].numgenes - curgenes;
      //c[0] = curid[coal[0]];
      //c[1] = curid[coal[1]];
      //addedge (ci, li, newedge, c, curgenes, &ptime,nosimmigration);
      if(coal.at(1) < curgenes-1)
	{
	  for (unsigned int i = coal.at(1); i < curgenes - 1; i++)
	    {
	      for (unsigned int j = 0; j < cursites; j++)
		Seqs.at(j).at(i) = Seqs.at(j).at(i + 1);
	      // curid[i] = curid[i + 1];
	    }
	}
      vector<int> newSeq;
      for(unsigned int i=0; i< cursites; i++)
	newSeq.push_back(Seqs.at(i).at(coal.at(0)));
      for(unsigned int i=coal.at(0); i<curgenes-1;i++)
	{
	  for (unsigned int j = 0; j < cursites; j++)
	    Seqs.at(j).at(i) = Seqs.at(j).at(i+1);
	}
      for(unsigned int i=0; i< cursites; i++)
	Seqs.at(i).at(curgenes-2) = newSeq.at(i);
      // curid[coal[0]] = newedge;
      curgenes--;
    }
  list_nodes[0]->set_isRoot(1);

  
  list_nodes[0]->assignPopulations2Tips(lc);

  return list_nodes[0];
 
  /* 
  C[ci]->G[li].roottime = ptime;
  C[ci]->G[li].root = L[li].numlines - 1;
  gtree[C[ci]->G[li].root].mig[0].mt = -1;
  gtree[C[ci]->G[li].root].time = TIMEMAX;
  orig2d_free2D ((void **) tempseq, L[li].numgenes);
  orig2d_free2D ((void **) distmat, L[li].numgenes);
  XFREE (curid);
  XFREE (singletons);
  C[ci]->G[li].hilike = -1e20;
  */
}                               /* makeIS */



